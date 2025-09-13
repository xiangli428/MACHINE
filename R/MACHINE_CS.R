#' @title Identify credible sets for MACHINE object.
#' 
#' @description 
#' A greedy search algorithm to identify all credible sets achieving target
#' "coverage" with minimum absolute correlation not less than "purity".
#' 
#' @usage
#' CS = MACHINE_CS(MACHINE, coverage = 0.95, purity = 0.5)
#' 
#' @param MACHINE A MACHINE object with MCMC samples.
#' @param coverage A number between 0 and 1 specifying the required coverage
#' of the credible sets.
#' @param purity A number between 0 and 1 specifying the minimum 
#' absolute correlation allowed in a credible set.
#' 
#' @return A list of lists. Each list is for an ancestry and comprises the 
#' following elements:
#' \describe{
#'   \item{sets}{A list in which each element is a vector containing the
#'   indices of the variables in the CS.}
#'   \item{CL}{The credible level of each CS.}
#'   \item{purity}{The purity of each CS.}
#' }
#' 
#' @export

MACHINE_CS = function(MACHINE, coverage = 0.95, purity = 0.5)
{
  P = MACHINE@P
  n = MACHINE@mcmc_samples$n_samples
  
  set_CL = function(beta, R, idx, n)
  {
    if(length(idx) >= 3)
    {
      eig = tryCatch(eigs_sym(R[idx,idx], k = 1),
                     error = function(e){eigen(R[idx,idx], symmetric = T)})
      v = beta[,idx] %*% eig$vectors[,1]
    } else {
      v = beta[,idx[1]] + beta[,idx[2]] * sign(R[idx[1],idx[2]])
    }
    abs(sum(v > 0) - sum(v < 0)) / n
  }
  
  CS = list()
  
  for(p in 1:P)
  {
    id = which(!is.na(MACHINE@Z[,p]))
    M = length(id)
    
    beta = t(MACHINE@mcmc_samples$beta[[p]][id,])
    
    R = as.matrix(MACHINE@R[[p]])[id,id]
    
    LD_pairs = foreach(j = 1:M) %do%
    {
      idx = which(abs(R[,j]) >= purity)
      idx[order(abs(R[idx,j]), decreasing = T)]
    }
    diag(R) = 1
    
    CL = data.frame("index" = 1:M, "CL" = MACHINE@CL[id,p])
    CL = arrange(CL, -CL, index)
    
    CS[[p]] = list("sets" = list(), "CL" = c())
    n1 = sum(CL$CL >= coverage)
    if(n1 > 0)
    {
      for(r in 1:n1)
      {
        CS[[p]]$sets[[r]] = CL$index[r]
        CS[[p]]$CL[r] = CL$CL[r]
      }
    }
    r = n1 + 1
    
    candidate = CL$index[(n1+1):M]
    candidate %<>% setdiff(which(apply(beta == 0, 2, all)))
    allSNPs = candidate
    CL = arrange(CL, index)
    
    while(length(candidate) > 0)
    {
      j = candidate[1]
      candidate = candidate[-1]
      tmp_CL = CL$CL[j]
      tests = intersect(LD_pairs[[j]], allSNPs)
      set = c(j)
      
      while(length(tests) > 0 & tmp_CL < coverage)
      {
        while(tmp_CL < coverage)
        {
          if(length(tests) > 0)
          {
            k = tests[which.max(CL$CL[tests])]
            tmp_CL = tmp_CL + CL$CL[k]
            tests = intersect(tests, LD_pairs[[k]])
            set %<>% append(k)
          } else {
            break
          }
        }
        tmp_CL = set_CL(beta, R, set, n)
      }
      
      if(tmp_CL >= coverage)
      {
        set_ori = set = set[order(CL$CL[set])]
        while(T)
        {
          for(l in set_ori)
          {
            if(length(set) > 2)
            {
              k = which(set == l)
              test_CL = set_CL(beta, R, set[-k], n)
              if(test_CL >= coverage)
              {
                set = set[-k]
                tmp_CL = test_CL
              }
            } else {
              break
            }
          }
          
          if(length(set_ori) == length(set) | length(set) == 2)
          {
            break
          } else {
            set_ori = set
          }
        }
        
        CS[[p]]$sets[[r]] = sort(set)
        CS[[p]]$CL[r] = tmp_CL
        r = r + 1
        candidate = setdiff(candidate, set)
        allSNPs = setdiff(allSNPs, set)
      }
    }
    
    CS[[p]]$purity = foreach(set = CS[[p]]$sets, .combine = "rbind") %do%
    {
      if(length(set) > 1)
      {
        R_sub = abs(R[set, set])
        R_sub = R_sub[upper.tri(R_sub)]
        data.frame("min.abs.corr" = min(R_sub),
                   "mean.abs.corr" = mean(R_sub),
                   "median.abs.corr" = median(R_sub))
      } else {
        data.frame("min.abs.corr" = 1,
                   "mean.abs.corr" = 1,
                   "median.abs.corr" = 1)
      }
    }
    
    CS[[p]]$sets = foreach(set = CS[[p]]$sets) %do%
    {
      set = id[set]
      names(set) = NULL
      set
    }
  }
  names(CS) = MACHINE@POP_ID
  return(CS)
}