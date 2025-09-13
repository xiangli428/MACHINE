# MACHINE v1.0

The `MACHINE` package implements a multi-ancestry fine-mapping method based on 
the "Multi-AnCestry Heritability INducEd Dirichlet decomposition" (MACHINE) prior.
It is a novel multi-ancestry fine-mapping method that utilizes a continuous 
global-local shrinkage prior.
This method can be applied to both quantitative and binary traits.
An MCMC algorithm is employed to obtain samples from the posterior 
distribution.
This method can also provide "credible sets" of candidate
causal variants, which are generally provided by fine-mapping methods
based on discrete-mixture priors.
We implement an approach for robust fine mapping with out-of-sample LD matrix.

## Quick Start

### 1. Install `MACHINE` from GitHub
```
devtools::install_github("https://github.com/xiangli428/MACHINE")
library(MACHINE)
```

### 2. load test dataset
```
data(MACHINE_test_data)
```
When out-of-sample LD matrix is used, eigendecomposition of each matrix of `R` 
will be performed during the creation of `MACHINE` object. By default, set 
`R_eig = NULL`. It can also be precomputed by `eigen`.
```
R_eig = foreach(k = 1:length(R)) %do%
{
  eigen(R[[k]], symmetric = TRUE)
}
# If the diagonal elements of each matrix of R are already set as 0, run
for(k in 1:length(R))
{
  R_eig[[k]]$values = R_eig[[k]]$values + 1
}
```

### 3. Create an MACHINE object with summary data

```
MACHINE = CreateMACHINEObject(Z,
                              R,
                              N,
                              SNP_ID,
                              POP_ID,
                              in_sample_LD = FALSE,
                              R_eig = R_eig,
                              a = 0.005,
                              c = 0.2)
```
where `a` and `c` specify the hyper-parameters of the prior.
By default, when both `in_sample = TRUE` and `in_sample = FALSE`, we 
recommend setting `b = NULL` and `b` will be estimated by a pre-training 
process before MCMC. 

### 4. MCMC sampling

```
MACHINE = MACHINE_MCMC(MACHINE, mcmc_n = 5500, burn_in = 500, thin = 1, 
                       stepsize = 2, seed = 42, get_CS = TRUE, 
                       n_chain = 3, fold = 1.1,
                       pre_mcmc_n = 400, pre_burn_in = 200, pre_p = 0.05,
                       pre_maxiter = 10, pre_miniter = 2)
```
We introduce parallel tempering to improve mixing of MCMC. Set
`n_chain = 1` can turn off parallel tempering.

### 5. Results

Credible levels of SNPs:
```
MACHINE@CL
```
Credible sets:
```
MACHINE@CS
```

### 6. Citation