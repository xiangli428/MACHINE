#' @title MACHINE package
#' @name MACHINE
#' @description MACHINE: Multi-AnCestry Heritability INducEd Dirichlet decomposition
#' 
#' @import Rcpp
#' @import RcppEigen
#' @import Matrix
#' @import foreach
#' @import dplyr
#' @import magrittr
#' @import RSpectra
#' 
#' @details 
#' This package implements a multi-ancestry fine-mapping method based on the 
#' "Multi-AnCestry Heritability INducEd Dirichlet decomposition" (MACHINE) prior.
#' It is a novel multi-ancestry fine-mapping method that utilizes a continuous 
#' global-local shrinkage prior.
#' This method can be applied to both quantitative and binary traits.
#' An MCMC algorithm is employed to obtain samples from the posterior 
#' distribution.
#' This method can also provide "credible sets" of candidate causal variants, 
#' which are generally provided by fine-mapping methods based on 
#' discrete-mixture priors.
#' 
#' @keywords internal
#' 
#' @references https://doi.org/10.1016/j.ajhg.2023.12.007

"_PACKAGE"

