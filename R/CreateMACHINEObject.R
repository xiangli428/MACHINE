#' @title Create a MACHINE object.
#' 
#' @description 
#' Create a MACHINE object from GWAS summary statistics and LD matrices of
#' multiple ancestries. GWAS z-scores, LD matrices, MACHINE hyper-parameters, 
#' MCMC samples, and outputs are stored.
#' 
#' @usage 
#' MACHINE = CreateMACHINEObject(Z,
#'                               R,
#'                               N,
#'                               SNP_ID = NULL,
#'                               POP_ID = NULL,
#'                               in_sample_LD = FALSE,
#'                               lambda = NULL,
#'                               R_eig = NULL,
#'                               a = 0.005,
#'                               b = NULL,
#'                               c = 0.2,
#'                               coverage = 0.95,
#'                               purity = 0.5,
#'                               tol = 1e-8)
#'
#' @param Z Z-score matrix with each row for a variant and each column for an
#' ancestry. Missing values are set as NA.
#' @param R A list of P correlation matrices that can be coerced to 
#' \code{\link[Matrix]{dsCMatrix-class}} format. The size of each matrix matches 
#' the number of variants present in the corresponding ancestry. The diagonal 
#' elements of each matrix will be coerced to zeros.
#' @param N GWAS sample size. For binary traits, 'N' is the sum of the
#' number of cases and the number of controls.
#' @param SNP_ID Identifiers of SNPs. The default is c("SNP_1", ...).
#' @param in_sample_LD Logical. Setting in_sample_LD = TRUE if in-sample LD
#' matrix is used.
#' @param lambda Hyper-parameter quantifying the discrepancy between z-scores
#' and R. If in_sample_LD == TRUE, lambda will be set as 0. Otherwise, if
#' lambda is not provided, lambda will be estimated.
#' @param R_eig An output of "eigen" function. A list with two components
#' "values" and "vector". If in_sample_LD == TRUE, don't need to provide it.
#' If provided, it will be used to create h2D2 object directly. Otherwise, 
#' it will be computed.
#' @param a Shape parameters for sigma^2. Either a positive real number or an
#' M-vector of positive real numbers.
#' @param b Shape parameter for 1-h^2. A positive real number. If not specified,
#' it will be estimated by an empirical Bayesian approach.
#' @param c Shape parameters for eta. Either a positive real number or an
#' M-by-P matrix of positive real numbers.
#' @param coverage A number between 0 and 1 specifying the required coverage
#' of the credible sets.
#' @param purity A number between 0 and 1 specifying the minimum 
#' absolute correlation allowed in a credible set.
#' @param tol Eigenvalues less than tol will be coerced to zeros.
#' 
#' @return An h2D2 object. See \code{\link{h2D2-class}}.
#' The input z-scores will be modified to z / sqrt(1 + z^2/N - 1/N).
#' 
#' @export

CreateMACHINEObject <- function(Z,
                                R,
                                N,
                                SNP_ID = NULL,
                                POP_ID = NULL,
                                in_sample_LD = F,
                                lambda = NULL,
                                R_eig = NULL,
                                a = 0.005,
                                b = NULL,
                                c = 0.2,
                                coverage = 0.95,
                                purity = 0.5,
                                tol = 1e-8)
{
  # Check z.
  Z = as.matrix(Z)
  P = ncol(Z)
  M = nrow(Z)
  if(M <= 1)
  {
    stop("Should contain at least 2 SNPs.")
  }
  
  # Check N.
  if(length(N) != P)
  {
    stop("The length of 'N' should be equal to the number of ancestries.")
  }
  if(any(!is.numeric(N)) | any(N < 0))
  {
    stop("'N' must be a vector of positive numbers.")
  }
  
  for(k in 1:P)
  {
    Z[,k] = Z[,k] / sqrt(1 + Z[,k]^2 / N[k] - 1 / N[k])
  }
  
  # Check SNP ID.
  if(is.null(SNP_ID))
  {
    SNP_ID = paste("SNP", 1:M, sep = "_")
  }
  if(length(SNP_ID) != M)
  {
    warning(sprintf("%s %s", "The length of 'SNP_ID' is not equal to the",
                    "number of SNPs. Rename SNPs."))
    SNP_ID = paste("SNP", 1:M, sep = "_")
  }
  
  # Check POP ID.
  if(is.null(POP_ID))
  {
    POP_ID = paste("POP", 1:P, sep = "_")
  }
  if(length(POP_ID) != P)
  {
    warning(sprintf("%s %s", "The length of 'POP_ID' is not equal to the",
                    "number of ancestries. Rename ancestries."))
    POP_ID = paste("POP", 1:P, sep = "_")
  }
  
  rownames(Z) = SNP_ID
  colnames(Z) = POP_ID
  names(N) = POP_ID
  
  # Check R.
  for(k in 1:P)
  {
    R[[k]] %<>% as.matrix()
    SNPs = SNP_ID[!is.na(Z[,k])]
    rownames(R[[k]]) = colnames(R[[k]]) = SNPs
    
    if(nrow(R[[k]]) != sum(!is.na(Z[,k])) | ncol(R[[k]]) != sum(!is.na(Z[,k])))
    {
      stop("Dimensions of LD matrix are not equal to the length of effect sizes.")
    }

    if(any(abs(R[[k]]) > 1))
    {
      # warning("There cannot be a value greater than 1 in the LD matrix.")
      # warning("Values larger than 1 or smaller than -1 are coerced to be 1 or -1.")
      R[[k]][R[[k]] > 1] = 1
      R[[k]][R[[k]] < -1] = -1
    }

    if(!isSymmetric(R[[k]], check.attributes = F))
    {
      warning("LD matrix must be symmetric. Coerce to be symmetric.")
      R[[k]] = (R[[k]] + t(R[[k]])) / 2
    }

    if(any(diag(R[[k]]) != 0))
    {
      diag(R[[k]]) = 0
      warning("The diagonal elements of LD matrix are coerced to zeros.")
    }
    
    R_full = matrix(0, M, M)
    rownames(R_full) = colnames(R_full) = SNP_ID
    R_full[SNPs, SNPs] = R[[k]]
    R[[k]] = R_full
    # R[[k]] = as(as(as(R_full, "dMatrix"), "symmetricMatrix"), "CsparseMatrix")
  }
  names(R) = POP_ID
  
  #Check lambda
  if(!in_sample_LD)
  {
    if(is.null(R_eig))
    {
      R_eig = foreach(k = 1:P) %do%
      {
        eigen(R[[k]][!is.na(Z[,k]), !is.na(Z[,k])], symmetric = T)
      }
      for(k in 1:P)
      {
        R_eig[[k]]$values = R_eig[[k]]$values + 1
      }
    }
    
    values = foreach(k = 1:P) %do%
    {
      R_eig[[k]]$values[R_eig[[k]]$values < tol] = 0
      R_eig[[k]]$values
    }
    tvectors = foreach(k = 1:P) %do%
    {
      t(R_eig[[k]]$vectors)
    }
    
    if(is.null(lambda))
    {
      lambda = foreach(k = 1:P, .combine = "c") %do%
      {
        idx = which(!is.na(Z[,k]))
        z = Z[idx,k]
        j = which.max(abs(z))
        dz = z - R[[k]][idx,idx[j]] * z[j]
        dz[j] = 0
        Udz2 = (tvectors[[k]] %*% dz)^2
        
        gr = function(lambda)
        {
          denom = lambda + values[[k]]
          sum(1 / denom - Udz2 / denom^2)
        }
        
        if(gr(1e-100) > 0)
        {
          opt = list(root = 0)
        } else {
          opt = uniroot(gr, interval = c(1e-100, 1e10), tol = 1e-10)
        }
        
        if(opt$root <= 1e-100)
        {
          max(mean(dz^2) - 1, 0)
        } else {
          max(mean(dz^2) - 1, opt$root)
        }
      }
    }
  } else {
    values = list()
    tvectors = list()
    lambda = rep(0,P)
  }
  
  names(lambda) = POP_ID
  
  for(k in 1:P)
  {
    R[[k]] = as(as(as(R[[k]], "dMatrix"), "symmetricMatrix"), "CsparseMatrix")
  }
  
  # Hyper-parameters
  if(is.null(a))
  {
    a = rep(0.005, M)
  } else if(length(a) == 1){
    a = rep(a, M)
  }
  if((length(a) != M) | any(a < 0))
  {
    stop("Shape parameters 'a' should be an M-vector of positive numbers.")
  }
  
  if(!is.null(b))
  {
    if(b <= 1)
    {
      warning("Setting b<=1 may lead to divergent result.")
    }
  }
  
  if(is.null(c))
  {
    c = matrix(0.2, M, P)
  } else if(length(c) == 1) {
    c = matrix(c, M, P)
  } else if(any(dim(c) != c(M,P)) | any(c < 0)) {
    stop("Shape parameters 'c' should be an M-by-P matrix of positive numbers.")
  }
  rownames(c) = SNP_ID
  colnames(c) = POP_ID
  
  # Check coverage and purity
  if(coverage < 0 | coverage > 1)
  {
    stop("Coverage must be in a range between 0 and 1.")
  }
  if(purity < 0 | purity > 1)
  {
    stop("Purity must be in a range between 0 and 1.")
  }
  
  return(new("MACHINE",
             M = M,
             P = P,
             N = N,
             Z = Z,
             R = R,
             SNP_ID = SNP_ID,
             POP_ID = POP_ID,
             values = values,
             tvectors = tvectors,
             lambda = lambda,
             a = a,
             b = b,
             c = c,
             mcmc_samples = list("n_samples" = 0, "n_burnin" = 0),
             mcmc_mean = list(),
             mcmc_sd = list(),
             mcmc_quantile = list(),
             CL = matrix(0, M, P),
             coverage = coverage,
             purity = purity,
             CS = list()))
}