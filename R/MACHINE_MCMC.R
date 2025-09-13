#' @title MCMC sampling under MACHINE prior.
#' 
#' @description
#' Using MCMC sampler to fit posterior distribution under MACHINE prior.
#' 
#' @usage
#' MACHINE = MACHINE_MCMC(MACHINE, mcmc_n = 5500, burn_in = 500, thin = 1, 
#'                        stepsize = 2, seed = 42, get_CS = TRUE, 
#'                        n_chain = 3, fold = 1.1,
#'                        pre_mcmc_n = 400, pre_burn_in = 200, pre_p = 0.05,
#'                        pre_maxiter = 10, pre_miniter = 2)
#' 
#' @param MACHINE A MACHINE object.
#' @param mcmc_n Number of posterior samples to generate.
#' @param burn_in Number of early samples to discard.
#' @param thin Thinning parameter of the MCMC chain. The default is 1, 
#' indicating no thinning.
#' @param stepsize Step size for proposal.
#' @param seed Random seed for MCMC sampling.
#' @param get_CS Logical. If get_CS = TRUE, detect CS after MCMC. Otherwise, 
#' CS will not be detected.
#' @param n_chain Number of parallel chains with different temperature.
#' @param fold Fold change of temperature between two neighbored parallel chains.
#' For instance, if \code{n_chain = 3} and \code{fold = 1.1}, the temperature
#' for the three parallel chains will be 1, 1.1, 1.21.
#' @param pre_mcmc_n Number of posterior samples to generate in each iteration
#' of pretraining. If b is given, the pretrain will not run.
#' @param pre_burn_in Number of early samples to discard in each iteration of
#' pretraining.
#' @param pre_p p-value threshold for the pretraining.
#' @param pre_maxiter Maximum number of pretraining iterations allowed.
#' @param pre_miniter Minimum number of pretraining iterations required.
#' 
#' @return An MACHINE object with MCMC samples. See \code{\link{MACHINE-class}}.
#' 
#' @import Matrix
#' @export

MACHINE_MCMC = function(MACHINE, mcmc_n = 5500, burn_in = 500, thin = 1,
                        stepsize = 1, seed = 42, get_CS = T, 
                        n_chain = 3, fold = 1.1,
                        pre_mcmc_n = 400, pre_burn_in = 200, pre_p = 0.05,
                        pre_maxiter = 10, pre_miniter = 2)
{
  if(!is.numeric(mcmc_n) | mcmc_n < 0)
  {
    stop("'mcmc_n' should be a non-negative integer.")
  }
  mcmc_n %<>% round()

  if(!is.numeric(burn_in) | burn_in < 0)
  {
    stop("'burn_in' should be a non-negative integer.")
  }
  burn_in %<>% round()

  if(!is.numeric(thin) | thin <= 0)
  {
    stop("'thin' should be a positive integer.")
  }
  thin %<>% round()
  thin %<>% max(1)

  if(!is.numeric(stepsize) | stepsize <= 0)
  {
    stop("'stepsize' should be a positive number.")
  }

  if(!is.numeric(seed) | seed < 0)
  {
    stop("'seed' should be a non-negative integer.")
  }
  seed %<>% round()

  if(!is.numeric(pre_mcmc_n) | mcmc_n < 0)
  {
    stop("'pre_mcmc_n' should be a non-negative integer.")
  }
  pre_mcmc_n %<>% round()

  if(!is.numeric(pre_burn_in) | pre_burn_in < 0)
  {
    stop("'pre_burn_in' should be a non-negative integer.")
  }
  pre_burn_in %<>% round()

  if(pre_mcmc_n <= pre_burn_in)
  {
    stop("'pre_mcmc_n' should be larger than 'pre_burn_in'.")
  }
  
  P = MACHINE@P
  M = MACHINE@M
  N = MACHINE@N
  lambda = MACHINE@lambda
  values = MACHINE@values
  tvectors = MACHINE@tvectors
  Z = MACHINE@Z
  
  obs = !is.na(Z)
  
  R = foreach(k = 1:P, .combine = "cbind") %do%
  {
    MACHINE@R[[k]]
  }
  R %<>% as("dgCMatrix")
  
  dW = matrix(0,M,P)
  
  W = foreach(k = 1:P, .combine = "cbind") %do%
  {
    if(lambda[k] == 0)
    {
      dW[,k] = N[k]
      MACHINE@R[[k]] * N[k]
    } else {
      w = t(tvectors[[k]]) %*% (
        tvectors[[k]] * values[[k]]^2 / (values[[k]] + lambda[k])) * N[k]
      w_full = matrix(0,M,M)
      w_full[obs[,k],obs[,k]] = (w+t(w)) / 2
      dW[,k] = diag(w_full)
      diag(w_full) = 0
      as(as(as(w_full, "dMatrix"), "generalMatrix"), "CsparseMatrix")
    }
  }
  
  mu = foreach(k = 1:P, .combine = "rbind") %do%
  {
    if(lambda[k] == 0)
    {
      sparseMatrix(i = rep(1, sum(obs[,k])), j = which(obs[,k]), 
                   x = Z[obs[,k],k] * sqrt(N[k]), dims = c(1,M))
    } else {
      z = Z[obs[,k],k]
      sparseMatrix(i = rep(1, sum(obs[,k])), j = which(obs[,k]), x = (
        t(tvectors[[k]]) %*% ((tvectors[[k]] %*% z) * values[[k]] / (
          values[[k]] + lambda[k])))[,1] * sqrt(N[k]), dims = c(1,M))
    }
  }
  
  LD_pairs_pop = foreach(k = 1:P, .combine = "cbind") %do%
  {
    LD = as.matrix(MACHINE@R[[k]]^2)
    LD[LD < MACHINE@purity^2] = 0
    LD = t(LD / pmax(rowSums(LD), 5e-324))
    as(LD, "dgCMatrix")
  }
  rm(LD)
  
  LD_pairs = foreach(k = 1:P, .combine = "pmax") %do%
  {
    as.matrix(MACHINE@R[[k]]^2)
  }
  LD_pairs[LD_pairs < MACHINE@purity^2] = 0
  LD_pairs = LD_pairs * (1 - as.matrix(dist(obs, method = "maximum")))
  LD_pairs = t(LD_pairs / (rowSums(LD_pairs) + 5e-324))
  LD_pairs %<>% as("dgCMatrix")
  
  sample = list()
  
  if(MACHINE@mcmc_samples$n_samples == 0)
  {
    sample[[1]] = matrix(0, M*P, n_chain) #beta
    sample[[2]] = foreach(n = 1:n_chain, .combine = "cbind") %do%
    {
      foreach(k = 1:P, .combine = "c") %do%
      {
        log(1e-4 / M / rowSums(obs) * obs[,k])
      }
    }
    sample[[3]] = matrix(1, M*P, n_chain) #psi
  } else {
    sample = MACHINE@mcmc_samples$last_sample
  }
  
  temp = fold^c(1:n_chain - 1)
  
  if(is.null(MACHINE@b))
  {
    asum_pop = foreach(k = 1:P, .combine = "+") %do%
    {
      id = which(obs[,k])
      sum(MACHINE@a[id] * rowSums(obs[id,]) * MACHINE@c[id,k] /
            rowSums(MACHINE@c[id,] * obs[id,], na.rm = T))
    } / P
    sample = MACHINE_pretrain(MACHINE, sample, obs, R, dW, W, mu, LD_pairs,
                              LD_pairs_pop, asum_pop, n_chain, temp,
                              pre_mcmc_n, pre_burn_in, pre_p,
                              pre_maxiter, pre_miniter, stepsize, seed)
  }
  
  if(mcmc_n > 0)
  {
    samples = MACHINE_sampling(MACHINE, sample, obs, R, dW, W, mu, 
                               LD_pairs, LD_pairs_pop, n_chain, temp, 
                               mcmc_n, thin, stepsize, seed)
    
    # Add new MCMC samples
    
    rownames(samples[[1]]) = rownames(samples[[2]]) = rownames(samples[[3]]) =
      rep(MACHINE@SNP_ID, P)
    rownames(samples[[4]]) = MACHINE@POP_ID
    
    colnames(samples[[1]]) = colnames(samples[[2]]) = colnames(samples[[3]]) =
      colnames(samples[[4]]) = names(samples[[5]]) = paste(
        "mcmc", (MACHINE@mcmc_samples$n_samples+MACHINE@mcmc_samples$n_burnin+1):
          (MACHINE@mcmc_samples$n_samples+MACHINE@mcmc_samples$n_burnin+mcmc_n),
        sep = '_')
    
    if(MACHINE@mcmc_samples$n_samples == 0)
    {
      for(k in 1:P)
      {
        MACHINE@mcmc_samples$beta[[k]] = samples[[1]][((k-1)*M+1):(k*M),]
        MACHINE@mcmc_samples$tau[[k]] = samples[[2]][((k-1)*M+1):(k*M),]
        MACHINE@mcmc_samples$psi[[k]] = samples[[3]][((k-1)*M+1):(k*M),]
      }
      names(MACHINE@mcmc_samples$beta) = names(MACHINE@mcmc_samples$tau) =
        names(MACHINE@mcmc_samples$psi) = MACHINE@POP_ID
      MACHINE@mcmc_samples$h2_beta = samples[[4]]
      MACHINE@mcmc_samples$h2 = samples[[5]]
    } else {
      for(k in 1:P)
      {
        MACHINE@mcmc_samples$beta[[k]] %<>% cbind(samples[[1]][((k-1)*M+1):(k*M),])
        MACHINE@mcmc_samples$tau[[k]] %<>% cbind(samples[[2]][((k-1)*M+1):(k*M),])
        MACHINE@mcmc_samples$psi[[k]] %<>% cbind(samples[[3]][((k-1)*M+1):(k*M),])
      }
      MACHINE@mcmc_samples$h2_beta %<>% cbind(samples[[4]])
      MACHINE@mcmc_samples$h2 %<>% c(samples[[5]])
    }
    
    MACHINE@mcmc_samples$last_sample = samples[[6]]
  }
  
  # burn in
  
  MACHINE@mcmc_samples$n_burnin = MACHINE@mcmc_samples$n_burnin + burn_in
  
  if(burn_in > 0)
  {
    for(l in c("beta", "tau", "psi"))
    {
      for(k in 1:P)
      {
        MACHINE@mcmc_samples[[l]][[k]] = 
          MACHINE@mcmc_samples[[l]][[k]][,-c(1:burn_in)]
      }
    }
    MACHINE@mcmc_samples$h2_beta = MACHINE@mcmc_samples$h2_beta[,-c(1:burn_in)]
    MACHINE@mcmc_samples$h2 = MACHINE@mcmc_samples$h2[-c(1:burn_in)]
  }
  
  MACHINE@mcmc_samples$n_samples = ncol(MACHINE@mcmc_samples$beta[[1]])
  
  #summary statistics
  MACHINE@mcmc_mean$beta = sapply(MACHINE@mcmc_samples$beta, rowMeans)
  MACHINE@mcmc_sd$beta = sapply(MACHINE@mcmc_samples$beta, function(x){
    return(apply(x, 1, sd))})
  MACHINE@mcmc_quantile$beta = foreach(k = 1:P) %do%
  {
    t(apply(MACHINE@mcmc_samples$beta[[k]], 1, quantile, 
            probs = c(0.025,0.25,0.5,0.75,0.975)))
  }
  names(MACHINE@mcmc_quantile$beta) = MACHINE@POP_ID
  
  sigma2 = lapply(MACHINE@mcmc_samples$tau, exp)
  MACHINE@mcmc_mean$sigma2 = sapply(sigma2, rowMeans)
  MACHINE@mcmc_sd$sigma2 = sapply(sigma2, function(x){return(apply(x, 1, sd))})
  MACHINE@mcmc_quantile$sigma2 = foreach(k = 1:P) %do%
  {
    t(apply(sigma2[[k]], 1, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))
  }
  names(MACHINE@mcmc_quantile$sigma2) = MACHINE@POP_ID
  
  MACHINE@mcmc_mean$h2_beta = rowMeans(MACHINE@mcmc_samples$h2_beta)
  MACHINE@mcmc_sd$h2_beta = apply(MACHINE@mcmc_samples$h2_beta, 1, sd)
  MACHINE@mcmc_quantile$h2_beta = t(apply(
    MACHINE@mcmc_samples$h2_beta, 1, quantile, 
    probs = c(0.025,0.25,0.5,0.75,0.975)))
  
  MACHINE@mcmc_mean$h2 = mean(MACHINE@mcmc_samples$h2)
  MACHINE@mcmc_sd$h2 = sd(MACHINE@mcmc_samples$h2)
  MACHINE@mcmc_quantile$h2 = quantile(MACHINE@mcmc_samples$h2, 
                                      probs = c(0.025,0.25,0.5,0.75,0.975))
  
  # CL and CS
  MACHINE@CL = sapply(MACHINE@mcmc_samples$beta, function(y){
    apply(y, 1, function(x){abs(sum(x > 0) - sum(x < 0)) / length(x)})})
  
  if(get_CS)
  {
    MACHINE@CS = MACHINE_CS(MACHINE, MACHINE@coverage, MACHINE@purity)
  }
  
  # Convergence
  n = MACHINE@mcmc_samples$n_samples
  n_2 = floor(n/2)
  MACHINE@mcmc_samples$PSRF_beta = foreach(k = 1:P, .combine = "cbind") %do%
  {
    foreach(j = 1:M, .combine = "c") %do%
    {
      s1 = MACHINE@mcmc_samples$beta[[k]][j,1:n_2]
      s2 = MACHINE@mcmc_samples$beta[[k]][j,(n-n_2+1):n]
      B = n_2 * (mean(s1) - mean(s2))^2 / 2
      W = (var(s1) + var(s2))/2
      sqrt(1 - 1/n_2 + B/(W+5e-324) / n_2)
    }
  }
  
  return(MACHINE)
}