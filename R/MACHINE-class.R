#' @title The MACHINE class
#' 
#' @description 
#' The MACHINE object is used to store GWAS summary data, LD matrices, variables
#' useful for MCMC sampling, MACHINE hyper-parameters, posterior samples, and
#' credible sets. Objects can be created by CreateMACHINEObject.
#' 
#' @slot M Number of SNPs.
#' @slot P Number of ancestries.
#' @slot N A P-vector of GWAS sample sizes for all ancestries. For binary traits, 
#' 'N' is the sum of the number of cases and the number of controls.
#' @slot Z Z-score matrix with M rows and P columns. Missing values are set as NA.
#' @slot R A list of P correlation matrices in \code{\link[Matrix]{dsCMatrix-class}} 
#' format. The size of each matrix matches the number of variants present in the
#' corresponding ancestry. The diagonal elements of each matrix are set as 0.
#' @slot SNP_ID IDs of SNPs. A character vector of length M.
#' @slot POP_ID IDs of ancestries. A character vector of length P.
#' @slot values A list of P vectors. Each vector comprises eigenvalues of 
#' corresponding R (with diagonal elements as 1).
#' @slot tvectors A list of P matrices. For each matrix, each column of its 
#' transpose is an eigenvector of corresponding R.
#' @slot lambda A P-vector of hyper-parameters quantifying the discrepancy 
#' between z-scores and R for each ancestry.
#' @slot a An M-vector of shape parameters for xi.
#' @slot b Shape parameter for 1-h^2.
#' @slot c An M-by-P matrix of shape parameters for eta.
#' @slot mcmc_samples A list storing MCMC samples.
#' @slot mcmc_mean A list storing posterior sample means.
#' @slot mcmc_sd A list storing posterior sample standard errors.
#' @slot mcmc_quantile A list storing posterior sample 2.5\%, 25\%, 50\%, 75\%, 
#' and 97.5\% quantiles.
#' @slot CL An M-by-P-matrix of credible levels of variables for each ancestry.
#' @slot coverage A number between 0 and 1 specifying the "coverage" of the 
#' credible sets.
#' @slot purity A number between 0 and 1 specifying the minimum absolute 
#' correlation allowed in a credible set.
#' @slot CS A list of identified credible sets in each ancestry.
#' See \code{\link{MACHINE_CS}}.
#' 
#' @export

setClass("MACHINE", 
         slots = c(
           M = "numeric", #Number of SNPs
           P = "numeric", #Number of populations
           N = "numeric", #Number of individuals
           Z = "matrix", #P*M matrix of z-scores
           R = "list", #LD matrix
           SNP_ID = "character",
           POP_ID = "character",
           values = "list",
           tvectors = "list",
           lambda = "numeric", #P inflation factors
           a = "numeric", #Shape parameters for averaged sigma2
           b = "ANY", #Shape parameter for 1-h2
           c = "matrix", #Shape parameters for ancestry-specific sigma2
           mcmc_samples = "list",
           mcmc_mean = "list",
           mcmc_sd = "list",
           mcmc_quantile = "list",
           CL = "matrix",
           coverage = "numeric",
           purity = "numeric",
           CS = "list"
         ))