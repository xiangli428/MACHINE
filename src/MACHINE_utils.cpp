// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>

#include <omp.h>
//[[Rcpp::plugins(openmp)]]

using namespace std;
using namespace Rcpp;
using namespace Eigen;

//' The inverse Gaussian distribution
//' 
//' @description Random generation for the inverse Gaussian distribution with 
//' parameters \\code{mu} and \\code{lambda}.
//'
//' @usage
//' rinvGauss(mu, lambda)
//' 
//' @param mu The mean parameter. Must be positive.
//' @param lambda The shape parameter. Must be positive.
//' 
//' @details
//' The inverse Gaussian distribution with parameters $mean = \\mu$ and 
//' $shape = \\lambda$ has density
//' $$
//' f(x) = \\sqrt\{\\frac\{\\lambda\}\{2 \\pi x^3\}\} 
//' \\exp \{ -\\frac\{\\lambda (x-\\mu)^2\}\{2 \\mu^2 x\} \}
//' $$
//' 
//' @return A random number from inverse Gaussian distribution.
// [[Rcpp::export]]
double rinvGauss(double nu, double lambda)
{
  double y = R::rchisq(1);
  double x;
  double u = 4 * lambda / (nu * y);
  if(u > 1e-11)
  {
    x = nu * (1 - 2 / (1 + sqrt(1 + u)));
  } else {
    x = lambda / y;
  }
  double z = R::runif(0,1);
  if(z <= nu / (nu + x))
    return(x);
  else
    return(pow(nu, 2)/x);
}

double propose_lognormal(double x, double stepsize = 2)
{
  x *= exp(R::rnorm(0, stepsize));
  return(x);
}

double lpexp(double x, double y)
{
  if(x == R_NegInf && y == R_NegInf)
  {
    return(R_NegInf);
  } else if(x > y)
  {
    return(x + log1p(exp(y-x)));
  } else {
    return(y + log1p(exp(x-y)));
  }
}

double lmexp(double x, double y)
{
  if(y == R_NegInf)
  {
    return(x);
  } else if (exp(y-x) == 1) {
    return(R_NegInf);
  } else {
    return(x + log1p(-exp(y-x)));
  }
}

// [[Rcpp::export]]
List MACHINE_pretrain(S4 MACHINE,
                      List sample,
                      const LogicalMatrix obs,
                      const Eigen::SparseMatrix<double> R,
                      const Eigen::MatrixXd dW,
                      const Eigen::SparseMatrix<double> W,
                      const Eigen::SparseMatrix<double> mu,
                      const Eigen::SparseMatrix<double> LD_pairs,
                      const Eigen::SparseMatrix<double> LD_pairs_pop,
                      const double asum_pop,
                      const unsigned int n_chain,
                      const NumericVector temp,
                      const unsigned int pre_mcmc_n = 400,
                      const unsigned int pre_burn_in = 200,
                      const double pre_p = 0.05,
                      const unsigned int pre_maxiter = 10,
                      const unsigned int pre_miniter = 2,
                      const double stepsize = 2,
                      const unsigned int seed = 428)
{
  Function set_seed("set.seed");
  set_seed(seed);

  // Data
  const int P = MACHINE.slot("P");
  const int M = MACHINE.slot("M");

  // Hyper parameters
  const NumericVector a = MACHINE.slot("a");
  const double asum = sum(a);
  const MatrixXd C = MACHINE.slot("c");
  const double sqrt2 = sqrt(2);
  double b = asum * (1 - 1e-4) * 1e4; // initial b
  double h2_est = 1e-4;
  NumericVector a_C(M);

  // Initial samples
  int i = 0; // index for mcmc iteration
  int ii; // index for thin
  int j, l, j1;//, j1; // index for SNP
  int j2 = 0;
  int k, k2; // index for pop
  int n; // index for chain

  for(j = 0; j < M; ++j)
  {
    a_C(j) = a(j);
    for(SparseMatrix<double>::InnerIterator itz(mu,j); itz; ++itz)
    {
      a_C(j) -= C(j,itz.row());
    }
  }

  MatrixXd beta(M*P, n_chain), tau(M*P, n_chain), psi(M*P, n_chain);
  VectorXd log_sigma2(M);
  NumericVector h2_beta(P);
  VectorXd log_res(n_chain);
  VectorXd ch_v(M*P);

  NumericVector K(M);
  for(j = 0; j < M; ++j)
  {
    K(j) = mu.col(j).nonZeros();
  }

  beta = sample[0];
  tau = sample[1];
  psi = sample[2];

  for(n = 0; n < n_chain; ++n)
  {
    log_res(n) = 0;
    for(j = 0; j < M; ++j)
    {
      for(SparseMatrix<double>::InnerIterator itz(mu,j); itz; ++itz)
      {
        k = itz.row();
        log_res(n) = lmexp(log_res(n), tau(k*M+j, n));
      }
    }
  }

  NumericVector keep_h2_beta(pre_mcmc_n - pre_burn_in);

  // mcmc
  double sigma2_jk, tau_jk_new, sigma2_jk_new, tsigma2_jk, tsigma2_jk_new;
  double u_jk, acc_prob;
  // double c_j, c_j1, c_j2;
  double log_res_j, nu_jk, rss, ch, r, u, q, W_12, p_21;
  double p_12 = 0;
  NumericVector order_sigma2;
  NumericVector s(P);

  NumericMatrix Z = MACHINE.slot("Z");
  NumericVector zmax(M);
  for(j = 0; j < M; ++j)
  {
    for(SparseMatrix<double>::InnerIterator itz(mu,j); itz; ++itz)
    {
      k = itz.row();
      if(abs(Z(j,k)) > zmax(j))
      {
        zmax(j) = abs(Z(j,k));
      }
    }
  }
  Function order("order");
  NumericVector order_SNP = order(zmax, _["decreasing"] = true);
  order_SNP = order_SNP - 1;

  Environment stats = Environment::namespace_env("stats");
  Function t_test = stats["t.test"];
  List test;
  double pval, estimate;

  bool dir = true;

  while(i < pre_maxiter)
  {
    ii = 0;
    while(ii < pre_mcmc_n)
    {
      for(l = 0; l < M; ++l)
      {
        j = order_SNP(l);

        for(SparseMatrix<double>::InnerIterator itz(mu,j); itz; ++itz)
        {
          k = itz.row();

          for(n = 0; n < n_chain; ++n)
          {
            //1. sigma2

            sigma2_jk = exp(tau(k*M+j,n));
            log_res(n) = lpexp(log_res(n), tau(k*M+j,n));

            tau_jk_new = propose_lognormal(
              tau(k*M+j,n) - log_res(n), R::runif(0,stepsize)) + log_res(n);
            acc_prob = (tau_jk_new - log_res(n)) / (tau(k*M+j,n) - log_res(n));
            sigma2_jk_new = exp(tau_jk_new);

            u_jk = itz.value();
            for(SparseMatrix<double>::InnerIterator it(W, k*M+j); it; ++it)
            {
              u_jk -= it.value() * beta(k*M+it.row(),n);
            }
            u_jk /= temp(n);

            tsigma2_jk = 1 / (dW(j,k) / temp(n) + 1 /
              (K(j) * sigma2_jk * psi(k*M+j,n)));
            tsigma2_jk_new = 1 / (dW(j,k) / temp(n) + 1 /
              (K(j) * sigma2_jk_new * psi(k*M+j,n)));

            acc_prob *= sqrt((K(j) * dW(j,k) * sigma2_jk * psi(k*M+j,n) + temp(n)) /
              (K(j) * dW(j,k) * sigma2_jk_new * psi(k*M+j,n) + temp(n))) *
                exp(pow(u_jk, 2) * (tsigma2_jk_new - tsigma2_jk) / 2 +
                (b-1) * (lmexp(log_res(n), tau_jk_new) -
                lmexp(log_res(n), tau(k*M+j,n))));

            if(K(j) == 1)
            {
              acc_prob *= exp(a(j) * (tau_jk_new - tau(k*M+j,n)));
            } else {
              log_res_j = R_NegInf;
              for(SparseMatrix<double>::InnerIterator itz2(mu,j);
                  itz2; ++itz2)
              {
                k2 = itz2.row();
                if(k2 != k)
                {
                  log_res_j = lpexp(log_res_j, tau(k2*M+j,n));
                }
              }
              acc_prob *= exp(C(j,k) * (tau_jk_new - tau(k*M+j,n)) +
                a_C(j) * (lpexp(log_res_j, tau_jk_new) -
                lpexp(log_res_j, tau(k*M+j,n))));
            }

            if((acc_prob > 1) || (R::runif(0,1) < acc_prob))
            {
              tau(k*M+j,n) = tau_jk_new;
              tsigma2_jk = tsigma2_jk_new;
            }
            log_res(n) = lmexp(log_res(n), tau(k*M+j,n));

            //2. beta
            beta(k*M+j,n) = R::rnorm(tsigma2_jk * u_jk, sqrt(tsigma2_jk));

            //3. psi
            if(beta(k*M+j,n) == 0)
            {
              psi(k*M+j,n) = 1;
            } else {
              nu_jk = sqrt2 * exp(tau(k*M+j,n)/2) / abs(beta(k*M+j,n));
              psi(k*M+j,n) = 1 / rinvGauss(nu_jk, 2);
            }
          }
        }
      }

      for(k = 0; k < P; ++k)
      {
        h2_beta(k) = 0;
        for(j = 0; j < M; ++j)
        {
          for(SparseMatrix<double>::InnerIterator it(R,k*M+j); it; ++it)
          {
            l = it.row();
            if(l > j)
            {
              h2_beta(k) += 2 * it.value() * beta(k*M+j,0) * beta(k*M+l,0);
            }
          }
          h2_beta(k) += pow(beta(k*M+j,0), 2);
        }
      }

      if(ii >= pre_burn_in)
      {
        keep_h2_beta(ii - pre_burn_in) = mean(h2_beta);
      }

      ++ii;

      if(((ii + 1) % 5 == 0) && (ii < pre_mcmc_n))
      {
        if(((ii + 1) % 15 == 5))
        {
          for(n = 0; n < n_chain; ++n)
          {
            for(j = 0; j < M; ++j)
            {
              log_sigma2(j) = R_NegInf;
              for(SparseMatrix<double>::InnerIterator itz(mu,j); itz; ++itz)
              {
                log_sigma2(j) = lpexp(log_sigma2(j), tau(itz.row()*M+j,n));
              }
            }

            //Choose the index for switching.
            order_sigma2 = order(log_sigma2, _["decreasing"] = true);
            order_sigma2 = order_sigma2 - 1;

            l = 10 + R::rgeom(0.3);
            if(l > M) l = M - 1;

            for(j = 0; j < l; ++j)
            {
              j1 = order_sigma2(j); //index of SNP 1
              if(LD_pairs.col(j1).nonZeros() >= 1)
              {
                r = R::runif(0,1);
                u = 0;
                for(SparseMatrix<double>::InnerIterator it(LD_pairs,j1); it; ++it)
                {
                  u += it.value();
                  if(u > r)
                  {
                    j2 = it.row();
                    p_12 = it.value();
                    break;
                  }
                }
                p_21 = LD_pairs.coeff(j1,j2);

                rss = 0;

                for(SparseMatrix<double>::InnerIterator itz(mu,j1); itz; ++itz)
                {
                  k = itz.row();
                  q = 0;
                  W_12 = W.coeff(j2,k*M+j1);
                  if(W_12 > 0)
                  {
                    s(k) = 1;
                    for(SparseMatrix<double>::InnerIterator it(W,k*M+j1); it; ++it)
                    {
                      q += it.value() * beta(k*M+it.row(),n);
                    }
                    for(SparseMatrix<double>::InnerIterator it(W,k*M+j2); it; ++it)
                    {
                      q -= it.value() * beta(k*M+it.row(),n);
                    }
                    q += W_12 * (beta(k*M+j1,n) - beta(k*M+j2,n));
                    q *= (beta(k*M+j1,n) - beta(k*M+j2,n));

                    q += (itz.value() - mu.coeff(k,j2)) *
                      (beta(k*M+j2,n) - beta(k*M+j1,n)) + (dW(j1,k) - dW(j2,k)) *
                      (pow(beta(k*M+j1,n),2) - pow(beta(k*M+j2,n),2)) / 2;
                  } else {
                    s(k) = -1;
                    for(SparseMatrix<double>::InnerIterator it(W,k*M+j1); it; ++it)
                    {
                      q += it.value() * beta(k*M+it.row(),n);
                    }
                    for(SparseMatrix<double>::InnerIterator it(W,k*M+j2); it; ++it)
                    {
                      q += it.value() * beta(k*M+it.row(),n);
                    }
                    q -= W_12 * (beta(k*M+j1,n) + beta(k*M+j2,n));
                    q *= (beta(k*M+j1,n) + beta(k*M+j2,n));

                    q += -(itz.value() + mu.coeff(k,j2)) *
                      (beta(k*M+j2,n) + beta(k*M+j1,n)) + (dW(j1,k) - dW(j2,k)) *
                      (pow(beta(k*M+j1,n),2) - pow(beta(k*M+j2,n),2)) / 2;
                  }

                  rss += q / temp(n);

                  if(K(j1) > 1)
                  {
                    rss += (C(j1,k) - C(j2,k)) * (tau(k*M+j2,n) - tau(k*M+j1,n));
                  }
                }

                if(K(j1) > 1)
                {
                  rss += (a_C(j1) - a_C(j2)) * (log_sigma2(j2) - log_sigma2(j1));
                } else {
                  rss += (a(j1) - a(j2)) * (log_sigma2(j2) - log_sigma2(j1));
                }

                acc_prob = exp(rss) * p_21 / p_12;

                if((acc_prob > 1) || (R::runif(0,1) < acc_prob))
                {
                  for(SparseMatrix<double>::InnerIterator itz(mu,j1); itz; ++itz)
                  {
                    k = itz.row();

                    ch = beta(k*M+j1,n);
                    beta(k*M+j1,n) = beta(k*M+j2,n) * s(k);
                    beta(k*M+j2,n) = ch * s(k);

                    ch = tau(k*M+j1,n);
                    tau(k*M+j1,n) = tau(k*M+j2,n);
                    tau(k*M+j2,n) = ch;

                    ch = psi(k*M+j1,n);
                    psi(k*M+j1,n) = psi(k*M+j2,n);
                    psi(k*M+j2,n) = ch;
                  }

                  ch = log_sigma2(j1);
                  log_sigma2(j1) = log_sigma2(j2);
                  log_sigma2(j2) = ch;
                }
              }
            }
          }
        } else if(((ii + 1) % 15 == 10)) {
          for(n = 0; n < n_chain; ++n)
          {
            for(k = 0; k < P; ++k)
            {
              for(j = 0; j < M; ++j)
              {
                log_sigma2(j) = tau(k*M+j,n) + log(K(j));
              }

              //Choose the index for switching.
              order_sigma2 = order(log_sigma2, _["decreasing"] = true);
              order_sigma2 = order_sigma2 - 1;

              l = 10 + R::rgeom(0.3);
              if(l > M) l = M - 1;

              for(j = 0; j < l; ++j)
              {
                j1 = order_sigma2(j); //index of SNP 1
                if(LD_pairs_pop.col(k*M+j1).nonZeros() >= 1)
                {
                  r = R::runif(0,1);
                  u = 0;
                  for(SparseMatrix<double>::InnerIterator it(LD_pairs_pop,k*M+j1);
                      it; ++it)
                  {
                    u += it.value();
                    if(u > r)
                    {
                      j2 = it.row();
                      p_12 = it.value();
                      break;
                    }
                  }

                  p_21 = LD_pairs_pop.coeff(j1,k*M+j2);

                  rss = 0;
                  W_12 = W.coeff(j2,k*M+j1);

                  if(W_12 > 0)
                  {
                    for(SparseMatrix<double>::InnerIterator it(W,k*M+j1); it; ++it)
                    {
                      rss += it.value() * beta(k*M+it.row(),n);
                    }
                    for(SparseMatrix<double>::InnerIterator it(W,k*M+j2); it; ++it)
                    {
                      rss -= it.value() * beta(k*M+it.row(),n);
                    }
                    rss += W_12 * (beta(k*M+j1,n) - beta(k*M+j2,n));
                    rss *= (beta(k*M+j1,n) - beta(k*M+j2,n));

                    rss += (mu.coeff(k,j1) - mu.coeff(k,j2)) *
                      (beta(k*M+j2,n) - beta(k*M+j1,n)) + (dW(j1,k) - dW(j2,k)) *
                      (pow(beta(k*M+j1,n),2) - pow(beta(k*M+j2,n),2)) / 2;
                  } else {
                    for(SparseMatrix<double>::InnerIterator it(W,k*M+j1); it; ++it)
                    {
                      rss += it.value() * beta(k*M+it.row(),n);
                    }
                    for(SparseMatrix<double>::InnerIterator it(W,k*M+j2); it; ++it)
                    {
                      rss += it.value() * beta(k*M+it.row(),n);
                    }
                    rss -= W_12 * (beta(k*M+j1,n) + beta(k*M+j2,n));
                    rss *= (beta(k*M+j1,n) + beta(k*M+j2,n));

                    rss += -(mu.coeff(k,j1) + mu.coeff(k,j2)) *
                      (beta(k*M+j2,n) + beta(k*M+j1,n)) + (dW(j1,k) - dW(j2,k)) *
                      (pow(beta(k*M+j1,n),2) - pow(beta(k*M+j2,n),2)) / 2;
                  }
                  rss /= temp(n);

                  rss += (C(j1,k) - C(j2,k)) * (tau(k*M+j2,n) - tau(k*M+j1,n));

                  log_res_j = R_NegInf;
                  for(SparseMatrix<double>::InnerIterator itz2(mu,j1);
                      itz2; ++itz2)
                  {
                    k2 = itz2.row();
                    if(k2 != k)
                    {
                      log_res_j = lpexp(log_res_j, tau(k2*M+j1,n));
                    }
                  }
                  rss += a_C(j1) * (lpexp(log_res_j, tau(k*M+j2,n)) -
                    lpexp(log_res_j, tau(k*M+j1,n)));

                  log_res_j = R_NegInf;
                  for(SparseMatrix<double>::InnerIterator itz2(mu,j2);
                      itz2; ++itz2)
                  {
                    k2 = itz2.row();
                    if(k2 != k)
                    {
                      log_res_j = lpexp(log_res_j, tau(k2*M+j2,n));
                    }
                  }
                  rss += a_C(j2) * (lpexp(log_res_j, tau(k*M+j1,n)) -
                    lpexp(log_res_j, tau(k*M+j2,n)));

                  rss += psi(k*M+j1,n) + psi(k*M+j2,n) - K(j1) * psi(k*M+j1,n) / K(j2) -
                    K(j2) * psi(k*M+j2,n) / K(j1);

                  acc_prob = exp(rss) * p_21 / p_12;

                  if((acc_prob > 1) || (R::runif(0,1) < acc_prob))
                  {
                    if(W_12 > 0)
                    {
                      ch = beta(k*M+j1,n);
                      beta(k*M+j1,n) = beta(k*M+j2,n);
                      beta(k*M+j2,n) = ch;
                    } else {
                      ch = beta(k*M+j1,n);
                      beta(k*M+j1,n) = -beta(k*M+j2,n);
                      beta(k*M+j2,n) = -ch;
                    }

                    ch = tau(k*M+j1,n);
                    tau(k*M+j1,n) = tau(k*M+j2,n);
                    tau(k*M+j2,n) = ch;

                    ch = psi(k*M+j1,n);
                    psi(k*M+j1,n) = K(j2) * psi(k*M+j2,n) / K(j1);
                    psi(k*M+j2,n) = K(j1) * ch / K(j2);
                  }
                }
              }
            }
          }
        } else {
          for(n = n_chain-2; n >= 0; --n)
          {
            rss = 0;
            for(j = 0; j < M; ++j)
            {
              for(SparseMatrix<double>::InnerIterator itz(mu,j); itz; ++itz)
              {
                k = itz.row();
                rss += itz.value() * beta(k*M+j,n) -
                  dW(j,k) * pow(beta(k*M+j,n),2) / 2;
                for(SparseMatrix<double>::InnerIterator it(W, k*M+j); it; ++it)
                {
                  l = it.row();
                  if(l > j)
                  {
                    rss -= it.value() * beta(k*M+j,n) * beta(k*M+l,n);
                  }
                }
              }
            }
            acc_prob = rss * (1 / temp(n+1) - 1 / temp(n));

            rss = 0;
            for(j = 0; j < M; ++j)
            {
              for(SparseMatrix<double>::InnerIterator itz(mu,j); itz; ++itz)
              {
                k = itz.row();
                rss += itz.value() * beta(k*M+j,n+1) -
                  dW(j,k) * pow(beta(k*M+j,n+1),2) / 2;
                for(SparseMatrix<double>::InnerIterator it(W, k*M+j); it; ++it)
                {
                  l = it.row();
                  if(l > j)
                  {
                    rss -= it.value() * beta(k*M+j,n+1) * beta(k*M+l,n+1);
                  }
                }
              }
            }
            acc_prob -= rss * (1 / temp(n+1) - 1 / temp(n));
            acc_prob = exp(acc_prob);

            if((acc_prob > 1) || (R::runif(0,1) < acc_prob))
            {
              ch_v = beta.col(n);
              beta.col(n) = beta.col(n+1);
              beta.col(n+1) = ch_v;

              ch_v = tau.col(n);
              tau.col(n) = tau.col(n+1);
              tau.col(n+1) = ch_v;

              ch_v = psi.col(n);
              psi.col(n) = psi.col(n+1);
              psi.col(n+1) = ch_v;

              ch = log_res(n);
              log_res(n) = log_res(n+1);
              log_res(n+1) = ch;
            }
          }
        }

        for(k = 0; k < P; ++k)
        {
          h2_beta(k) = 0;
          for(j = 0; j < M; ++j)
          {
            for(SparseMatrix<double>::InnerIterator it(R,k*M+j); it; ++it)
            {
              l = it.row();
              if(l > j)
              {
                h2_beta(k) += 2 * it.value() * beta(k*M+j,0) * beta(k*M+l,0);
              }
            }
            h2_beta(k) += pow(beta(k*M+j,0), 2);
          }
        }

        if(ii >= pre_burn_in)
        {
          keep_h2_beta(ii - pre_burn_in) = mean(h2_beta);
        }

        ++ii;
      }
    }

    estimate = mean(keep_h2_beta);
    if(estimate >= 1e-6 && estimate < 1)
    {
      test = t_test(_["x"] = keep_h2_beta, _["mu"] = h2_est);
      pval = test["p.value"];
      if(i == pre_miniter) dir = estimate > h2_est;

      cout << "Pretrain step " << i << ", b=" << b << ", h2_est=" <<
        h2_est << ", mean=" << estimate << ", p=" << pval <<
          "." << endl;

      if(((pval >= pre_p) || (estimate > h2_est) != dir) && i >= pre_miniter)
      {
        break;
      } else {
        b = asum_pop / estimate - asum;
        h2_est = estimate;
      }
    } else {
      estimate = 1e-6;
      b = asum_pop / estimate - asum;
      h2_est = estimate;
    }

    ++i;
  }

  cout << "End pretrain. b=" << b << "." << endl;

  MACHINE.slot("b") = b;

  List last_sample = List::create(beta, tau, psi);
  return(last_sample);
}

// [[Rcpp::export]]
List MACHINE_sampling(S4 MACHINE,
                      List sample,
                      const LogicalMatrix obs,
                      const Eigen::SparseMatrix<double> R,
                      const Eigen::MatrixXd dW,
                      const Eigen::SparseMatrix<double> W,
                      const Eigen::SparseMatrix<double> mu,
                      const Eigen::SparseMatrix<double> LD_pairs,
                      const Eigen::SparseMatrix<double> LD_pairs_pop,
                      const int n_chain,
                      const NumericVector temp,
                      const int mcmc_n = 100, 
                      const int thin = 1,
                      const double stepsize = 2, 
                      const unsigned int seed = 428)
{
  Function set_seed("set.seed");
  set_seed(seed);
  
  // Data
  const int P = MACHINE.slot("P");
  const int M = MACHINE.slot("M");
  
  // Hyper parameters
  const VectorXd a = MACHINE.slot("a");
  const double b = MACHINE.slot("b");
  const MatrixXd C = MACHINE.slot("c");
  const double sqrt2 = sqrt(2);
  NumericVector a_C(M);
  
  // Initial samples
  int i = 0; // index for mcmc iteration
  int ii; // index for thin
  int j, l, j1;//, j1; // index for SNP
  int j2 = 0;
  int k, k2; // index for pop
  int n; // index for chain
  
  for(j = 0; j < M; ++j)
  {
    a_C(j) = a(j);
    for(SparseMatrix<double>::InnerIterator itz(mu,j); itz; ++itz)
    {
      a_C(j) -= C(j,itz.row());
    }
  }
  
  MatrixXd beta(M*P, n_chain), tau(M*P, n_chain), psi(M*P, n_chain);
  VectorXd log_sigma2(M);
  VectorXd h2_beta(P), sqrt_h2_beta(P);
  VectorXd log_res(n_chain);
  VectorXd ch_v(M*P);
  
  NumericVector K(M);
  for(j = 0; j < M; ++j)
  {
    K(j) = mu.col(j).nonZeros();
  }
  
  beta = sample[0];
  tau = sample[1];
  psi = sample[2];
  
  for(n = 0; n < n_chain; ++n)
  {
    log_res(n) = 0;
    for(j = 0; j < M; ++j)
    {
      for(SparseMatrix<double>::InnerIterator itz(mu,j); itz; ++itz)
      {
        k = itz.row();
        log_res(n) = lmexp(log_res(n), tau(k*M+j, n));
      }
    }
  }
  
  for(k = 0; k < P; ++k)
  {
    h2_beta(k) = 0;
    for(j = 0; j < M; ++j)
    {
      for(SparseMatrix<double>::InnerIterator it(R,k*M+j); it; ++it)
      {
        l = it.row();
        if(l > j)
        {
          h2_beta(k) += 2 * it.value() * beta(k*M+j,0) * beta(k*M+l,0);
        }
      }
      h2_beta(k) += pow(beta(k*M+j,0), 2);
    }
    sqrt_h2_beta(k) = 1;
  }
  
  MatrixXd keep_beta = MatrixXd::Zero(P*M,mcmc_n);
  MatrixXd keep_tau = MatrixXd::Zero(P*M,mcmc_n);
  MatrixXd keep_psi = MatrixXd::Zero(P*M,mcmc_n);
  MatrixXd keep_h2_beta = MatrixXd::Zero(P,mcmc_n);
  VectorXd keep_h2 = VectorXd::Zero(mcmc_n);
  
  // mcmc
  double sigma2_jk, tau_jk_new, sigma2_jk_new, tsigma2_jk, tsigma2_jk_new;
  double u_jk, acc_prob;
  double log_res_j, nu_jk, rss, ch, r, u, q, W_12, p_21;
  double p_12 = 0;
  NumericVector order_sigma2;
  NumericVector s(P);
  
  NumericMatrix Z = MACHINE.slot("Z");
  NumericVector zmax(M);
  for(j = 0; j < M; ++j)
  {
    for(SparseMatrix<double>::InnerIterator itz(mu,j); itz; ++itz)
    {
      k = itz.row();
      if(abs(Z(j,k)) > zmax(j))
      {
        zmax(j) = abs(Z(j,k));
      }
    }
  }
  Function order("order");
  NumericVector order_SNP = order(zmax, _["decreasing"] = true);
  order_SNP = order_SNP - 1;
  
  while(i < mcmc_n)
  {
    for(ii = 0; ii < thin; ++ii)
    {
      for(l = 0; l < M; ++l)
      {
        j = order_SNP(l);
        
        for(SparseMatrix<double>::InnerIterator itz(mu,j); itz; ++itz)
        {
          k = itz.row();
          
          for(n = 0; n < n_chain; ++n)
          {
            //1. sigma2
            
            sigma2_jk = exp(tau(k*M+j,n));
            log_res(n) = lpexp(log_res(n), tau(k*M+j,n));
            
            tau_jk_new = propose_lognormal(
              tau(k*M+j,n) - log_res(n), R::runif(0,stepsize)) + log_res(n);
            acc_prob = (tau_jk_new - log_res(n)) / (tau(k*M+j,n) - log_res(n));
            sigma2_jk_new = exp(tau_jk_new);
            
            u_jk = itz.value();
            for(SparseMatrix<double>::InnerIterator it(W, k*M+j); it; ++it)
            {
              u_jk -= it.value() * beta(k*M+it.row(),n);
            }
            u_jk /= temp(n);
            
            tsigma2_jk = 1 / (dW(j,k) / temp(n) + 1 / 
              (K(j) * sigma2_jk * psi(k*M+j,n)));
            tsigma2_jk_new = 1 / (dW(j,k) / temp(n) + 1 / 
              (K(j) * sigma2_jk_new * psi(k*M+j,n)));
            
            acc_prob *= sqrt((K(j) * dW(j,k) * sigma2_jk * psi(k*M+j,n) + temp(n)) /
              (K(j) * dW(j,k) * sigma2_jk_new * psi(k*M+j,n) + temp(n))) *
                exp(pow(u_jk, 2) * (tsigma2_jk_new - tsigma2_jk) / 2 +
                (b-1) * (lmexp(log_res(n), tau_jk_new) - 
                lmexp(log_res(n), tau(k*M+j,n))));
            
            if(K(j) == 1)
            {
              acc_prob *= exp(a(j) * (tau_jk_new - tau(k*M+j,n)));
            } else {
              log_res_j = R_NegInf;
              for(SparseMatrix<double>::InnerIterator itz2(mu,j); 
                  itz2; ++itz2)
              {
                k2 = itz2.row();
                if(k2 != k)
                {
                  log_res_j = lpexp(log_res_j, tau(k2*M+j,n));
                }
              }
              acc_prob *= exp(C(j,k) * (tau_jk_new - tau(k*M+j,n)) +
                a_C(j) * (lpexp(log_res_j, tau_jk_new) -
                lpexp(log_res_j, tau(k*M+j,n))));
            }
            
            if((acc_prob > 1) || (R::runif(0,1) < acc_prob))
            {
              tau(k*M+j,n) = tau_jk_new;
              // sigma2_jk = sigma2_jk_new;
              tsigma2_jk = tsigma2_jk_new;
            }
            log_res(n) = lmexp(log_res(n), tau(k*M+j,n));
            
            //2. beta
            beta(k*M+j,n) = R::rnorm(tsigma2_jk * u_jk, sqrt(tsigma2_jk));
            // if(abs(beta(k*M+j,n)) > 2*sqrt_h2_beta(k) && abs(beta(k*M+j,n)) > 0.01)
            // {
            //   beta(k*M+j,n) = 0;
            // }
            
            //3. psi
            if(beta(k*M+j,n) == 0)
            {
              psi(k*M+j,n) = 1;
            } else {
              nu_jk = sqrt2 * exp(tau(k*M+j,n)/2) / abs(beta(k*M+j,n));
              psi(k*M+j,n) = 1 / rinvGauss(nu_jk, 2);
            }
          }
        }
      }
      
      for(k = 0; k < P; ++k)
      {
        h2_beta(k) = 0;
        for(j = 0; j < M; ++j)
        {
          for(SparseMatrix<double>::InnerIterator it(R,k*M+j); it; ++it)
          {
            l = it.row();
            if(l > j)
            {
              h2_beta(k) += 2 * it.value() * beta(k*M+j,0) * beta(k*M+l,0);
            }
          }
          h2_beta(k) += pow(beta(k*M+j,0), 2);
        }
        sqrt_h2_beta(k) = sqrt(h2_beta(k));
      }
    }
    
    keep_beta.col(i) = beta.col(0);
    keep_tau.col(i) = tau.col(0);
    keep_psi.col(i) = psi.col(0);
    keep_h2_beta.col(i) = h2_beta;
    keep_h2(i) = 1 - exp(log_res(0));
    
    ++i;
    
    if(((i + 1) % 5 == 0) && (i < mcmc_n))
    {
      if(((i + 1) % 15 == 5))
      {
        for(n = 0; n < n_chain; ++n)
        {
          for(j = 0; j < M; ++j)
          {
            log_sigma2(j) = R_NegInf;
            for(SparseMatrix<double>::InnerIterator itz(mu,j); itz; ++itz)
            {
              log_sigma2(j) = lpexp(log_sigma2(j), tau(itz.row()*M+j,n));
            }
          }
          
          //Choose the index for switching.
          order_sigma2 = order(log_sigma2, _["decreasing"] = true);
          order_sigma2 = order_sigma2 - 1;
          
          l = 10 + R::rgeom(0.3);
          if(l > M) l = M - 1;
          
          for(j = 0; j < l; ++j)
          {
            j1 = order_sigma2(j); //index of SNP 1
            if(LD_pairs.col(j1).nonZeros() >= 1)
            {
              r = R::runif(0,1);
              u = 0;
              for(SparseMatrix<double>::InnerIterator it(LD_pairs,j1); it; ++it)
              {
                u += it.value();
                if(u > r)
                {
                  j2 = it.row();
                  p_12 = it.value();
                  break;
                }
              }
              p_21 = LD_pairs.coeff(j1,j2);
              
              rss = 0;
              
              for(SparseMatrix<double>::InnerIterator itz(mu,j1); itz; ++itz)
              {
                k = itz.row();
                q = 0;
                W_12 = W.coeff(j2,k*M+j1);
                if(W_12 > 0)
                {
                  s(k) = 1;
                  for(SparseMatrix<double>::InnerIterator it(W,k*M+j1); it; ++it)
                  {
                    q += it.value() * beta(k*M+it.row(),n);
                  }
                  for(SparseMatrix<double>::InnerIterator it(W,k*M+j2); it; ++it)
                  {
                    q -= it.value() * beta(k*M+it.row(),n);
                  }
                  q += W_12 * (beta(k*M+j1,n) - beta(k*M+j2,n));
                  q *= (beta(k*M+j1,n) - beta(k*M+j2,n));
                  
                  q += (itz.value() - mu.coeff(k,j2)) *
                    (beta(k*M+j2,n) - beta(k*M+j1,n)) + (dW(j1,k) - dW(j2,k)) *
                    (pow(beta(k*M+j1,n),2) - pow(beta(k*M+j2,n),2)) / 2;
                } else {
                  s(k) = -1;
                  for(SparseMatrix<double>::InnerIterator it(W,k*M+j1); it; ++it)
                  {
                    q += it.value() * beta(k*M+it.row(),n);
                  }
                  for(SparseMatrix<double>::InnerIterator it(W,k*M+j2); it; ++it)
                  {
                    q += it.value() * beta(k*M+it.row(),n);
                  }
                  q -= W_12 * (beta(k*M+j1,n) + beta(k*M+j2,n));
                  q *= (beta(k*M+j1,n) + beta(k*M+j2,n));
                  
                  q += -(itz.value() + mu.coeff(k,j2)) *
                    (beta(k*M+j2,n) + beta(k*M+j1,n)) + (dW(j1,k) - dW(j2,k)) *
                    (pow(beta(k*M+j1,n),2) - pow(beta(k*M+j2,n),2)) / 2;
                }
                
                rss += q / temp(n);
                
                if(K(j1) > 1)
                {
                  rss += (C(j1,k) - C(j2,k)) * (tau(k*M+j2,n) - tau(k*M+j1,n));
                }
              }
              
              if(K(j1) > 1)
              {
                rss += (a_C(j1) - a_C(j2)) * (log_sigma2(j2) - log_sigma2(j1));
              } else {
                rss += (a(j1) - a(j2)) * (log_sigma2(j2) - log_sigma2(j1));
              }
              
              acc_prob = exp(rss) * p_21 / p_12;
              
              if((acc_prob > 1) || (R::runif(0,1) < acc_prob))
              {
                for(SparseMatrix<double>::InnerIterator itz(mu,j1); itz; ++itz)
                {
                  k = itz.row();
                  
                  ch = beta(k*M+j1,n);
                  beta(k*M+j1,n) = beta(k*M+j2,n) * s(k);
                  beta(k*M+j2,n) = ch * s(k);
                  
                  ch = tau(k*M+j1,n);
                  tau(k*M+j1,n) = tau(k*M+j2,n);
                  tau(k*M+j2,n) = ch;
                  
                  ch = psi(k*M+j1,n);
                  psi(k*M+j1,n) = psi(k*M+j2,n);
                  psi(k*M+j2,n) = ch;
                }
                
                ch = log_sigma2(j1);
                log_sigma2(j1) = log_sigma2(j2);
                log_sigma2(j2) = ch;
              }
            }
          }
        }
      } else if(((i + 1) % 15 == 10)) {
        for(n = 0; n < n_chain; ++n)
        {
          for(k = 0; k < P; ++k)
          {
            for(j = 0; j < M; ++j)
            {
              log_sigma2(j) = tau(k*M+j,n) + log(K(j));
            }
            
            //Choose the index for switching.
            order_sigma2 = order(log_sigma2, _["decreasing"] = true);
            order_sigma2 = order_sigma2 - 1;
            
            l = 10 + R::rgeom(0.3);
            if(l > M) l = M - 1;
            
            for(j = 0; j < l; ++j)
            {
              j1 = order_sigma2(j); //index of SNP 1
              if(LD_pairs_pop.col(k*M+j1).nonZeros() >= 1)
              {
                r = R::runif(0,1);
                u = 0;
                for(SparseMatrix<double>::InnerIterator it(LD_pairs_pop,k*M+j1);
                    it; ++it)
                {
                  u += it.value();
                  if(u > r)
                  {
                    j2 = it.row();
                    p_12 = it.value();
                    break;
                  }
                }
                
                p_21 = LD_pairs_pop.coeff(j1,k*M+j2);
                
                rss = 0;
                W_12 = W.coeff(j2,k*M+j1);
                
                if(W_12 > 0)
                {
                  for(SparseMatrix<double>::InnerIterator it(W,k*M+j1); it; ++it)
                  {
                    rss += it.value() * beta(k*M+it.row(),n);
                  }
                  for(SparseMatrix<double>::InnerIterator it(W,k*M+j2); it; ++it)
                  {
                    rss -= it.value() * beta(k*M+it.row(),n);
                  }
                  rss += W_12 * (beta(k*M+j1,n) - beta(k*M+j2,n));
                  rss *= (beta(k*M+j1,n) - beta(k*M+j2,n));
                  
                  rss += (mu.coeff(k,j1) - mu.coeff(k,j2)) *
                    (beta(k*M+j2,n) - beta(k*M+j1,n)) + (dW(j1,k) - dW(j2,k)) *
                    (pow(beta(k*M+j1,n),2) - pow(beta(k*M+j2,n),2)) / 2;
                } else {
                  for(SparseMatrix<double>::InnerIterator it(W,k*M+j1); it; ++it)
                  {
                    rss += it.value() * beta(k*M+it.row(),n);
                  }
                  for(SparseMatrix<double>::InnerIterator it(W,k*M+j2); it; ++it)
                  {
                    rss += it.value() * beta(k*M+it.row(),n);
                  }
                  rss -= W_12 * (beta(k*M+j1,n) + beta(k*M+j2,n));
                  rss *= (beta(k*M+j1,n) + beta(k*M+j2,n));
                  
                  rss += -(mu.coeff(k,j1) + mu.coeff(k,j2)) *
                    (beta(k*M+j2,n) + beta(k*M+j1,n)) + (dW(j1,k) - dW(j2,k)) *
                    (pow(beta(k*M+j1,n),2) - pow(beta(k*M+j2,n),2)) / 2;
                }
                rss /= temp(n);
                
                rss += (C(j1,k) - C(j2,k)) * (tau(k*M+j2,n) - tau(k*M+j1,n));
                
                log_res_j = R_NegInf;
                for(SparseMatrix<double>::InnerIterator itz2(mu,j1);
                    itz2; ++itz2)
                {
                  k2 = itz2.row();
                  if(k2 != k)
                  {
                    log_res_j = lpexp(log_res_j, tau(k2*M+j1,n));
                  }
                }
                rss += a_C(j1) * (lpexp(log_res_j, tau(k*M+j2,n)) -
                  lpexp(log_res_j, tau(k*M+j1,n)));
                
                log_res_j = R_NegInf;
                for(SparseMatrix<double>::InnerIterator itz2(mu,j2);
                    itz2; ++itz2)
                {
                  k2 = itz2.row();
                  if(k2 != k)
                  {
                    log_res_j = lpexp(log_res_j, tau(k2*M+j2,n));
                  }
                }
                rss += a_C(j2) * (lpexp(log_res_j, tau(k*M+j1,n)) -
                  lpexp(log_res_j, tau(k*M+j2,n)));
                
                rss += psi(k*M+j1,n) + psi(k*M+j2,n) - K(j1) * psi(k*M+j1,n) / K(j2) -
                  K(j2) * psi(k*M+j2,n) / K(j1);
                
                acc_prob = exp(rss) * p_21 / p_12;
                
                if((acc_prob > 1) || (R::runif(0,1) < acc_prob))
                {
                  if(W_12 > 0)
                  {
                    ch = beta(k*M+j1,n);
                    beta(k*M+j1,n) = beta(k*M+j2,n);
                    beta(k*M+j2,n) = ch;
                  } else {
                    ch = beta(k*M+j1,n);
                    beta(k*M+j1,n) = -beta(k*M+j2,n);
                    beta(k*M+j2,n) = -ch;
                  }
                  
                  ch = tau(k*M+j1,n);
                  tau(k*M+j1,n) = tau(k*M+j2,n);
                  tau(k*M+j2,n) = ch;
                  
                  ch = psi(k*M+j1,n);
                  psi(k*M+j1,n) = K(j2) * psi(k*M+j2,n) / K(j1);
                  psi(k*M+j2,n) = K(j1) * ch / K(j2);
                }
              }
            }
          }
        }
      } else {
        for(n = n_chain-2; n >= 0; --n)
        {
          rss = 0;
          for(j = 0; j < M; ++j)
          {
            for(SparseMatrix<double>::InnerIterator itz(mu,j); itz; ++itz)
            {
              k = itz.row();
              rss += itz.value() * beta(k*M+j,n) - 
                dW(j,k) * pow(beta(k*M+j,n),2) / 2;
              for(SparseMatrix<double>::InnerIterator it(W, k*M+j); it; ++it)
              {
                l = it.row();
                if(l > j)
                {
                  rss -= it.value() * beta(k*M+j,n) * beta(k*M+l,n);
                }
              }
            }
          }
          acc_prob = rss * (1 / temp(n+1) - 1 / temp(n));
          
          rss = 0;
          for(j = 0; j < M; ++j)
          {
            for(SparseMatrix<double>::InnerIterator itz(mu,j); itz; ++itz)
            {
              k = itz.row();
              rss += itz.value() * beta(k*M+j,n+1) - 
                dW(j,k) * pow(beta(k*M+j,n+1),2) / 2;
              for(SparseMatrix<double>::InnerIterator it(W, k*M+j); it; ++it)
              {
                l = it.row();
                if(l > j)
                {
                  rss -= it.value() * beta(k*M+j,n+1) * beta(k*M+l,n+1);
                }
              }
            }
          }
          acc_prob -= rss * (1 / temp(n+1) - 1 / temp(n));
          acc_prob = exp(acc_prob);
          
          if((acc_prob > 1) || (R::runif(0,1) < acc_prob))
          {
            // cout << "Switch chain " << n << " and chain " << n+1 << 
            //   ", acc_prob=" << acc_prob << "." << endl;
            
            ch_v = beta.col(n);
            beta.col(n) = beta.col(n+1);
            beta.col(n+1) = ch_v;
            
            ch_v = tau.col(n);
            tau.col(n) = tau.col(n+1);
            tau.col(n+1) = ch_v;
            
            ch_v = psi.col(n);
            psi.col(n) = psi.col(n+1);
            psi.col(n+1) = ch_v;
            
            ch = log_res(n);
            log_res(n) = log_res(n+1);
            log_res(n+1) = ch;
          }
        }
      }

      for(k = 0; k < P; ++k)
      {
        h2_beta(k) = 0;
        for(j = 0; j < M; ++j)
        {
          for(SparseMatrix<double>::InnerIterator it(R,k*M+j); it; ++it)
          {
            l = it.row();
            if(l > j)
            {
              h2_beta(k) += 2 * it.value() * beta(k*M+j,0) * beta(k*M+l,0);
            }
          }
          h2_beta(k) += pow(beta(k*M+j,0), 2);
        }
        sqrt_h2_beta(k) = sqrt(h2_beta(k));
      }

      keep_beta.col(i) = beta.col(0);
      keep_tau.col(i) = tau.col(0);
      keep_psi.col(i) = psi.col(0);
      keep_h2_beta.col(i) = h2_beta;
      keep_h2(i) = 1 - exp(log_res(0));

      ++i;
    }
  }
  
  List last_sample = List::create(beta, tau, psi);
  
  return(List::create(keep_beta, keep_tau, keep_psi, keep_h2_beta, keep_h2,
                      last_sample));
}