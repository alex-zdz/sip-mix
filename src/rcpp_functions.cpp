#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// Univariate functions:

const double log2pi = std::log(2.0 * M_PI);

// Function to compute univariate normal density
// [[Rcpp::export]]
arma::vec dnorm_arma_uv(arma::vec const &x, 
                     double const mean, 
                     double const sigma, 
                     bool const logd = false) { 
  arma::vec const distval = (x - mean) / sigma;
  arma::vec const logretval = -(log2pi + std::log(sigma) + arma::square(distval)) / 2;
  
  if (logd)
    return logretval;
  return arma::exp(logretval);
}

// [[Rcpp::export]]
List sample_allocs_uv(int N, int M_samp, int M_a_samp, int M_na_samp, 
                   arma::vec allocs_samp, arma::vec weights_samp, 
                   arma::vec y, arma::vec mus_a_samp, 
                   arma::vec sigmas_a_samp, arma::vec mus_na_samp, 
                   arma::vec sigmas_na_samp, arma::vec n_curr) {
  
  arma::vec alloc_rates(M_samp);
  for (int i = 0; i < N; ++i) {
    int old_alloc = allocs_samp[i];
    for (int c = 0; c < M_a_samp; ++c) {
      alloc_rates[c] = weights_samp[c] * dnorm_arma_uv(y.subvec(i, i), mus_a_samp[c], sigmas_a_samp[c])(0);
    }
    if (M_na_samp != 0) {
      for (int c = 0; c < M_na_samp; ++c) {
        alloc_rates[M_a_samp + c] = weights_samp[c] * dnorm_arma_uv(y.subvec(i, i), mus_na_samp[c], sigmas_na_samp[c])(0);
      }
    }
    alloc_rates /= arma::sum(alloc_rates);
    arma::vec cum_probs = arma::cumsum(alloc_rates);
    double u = R::runif(0, 1);
    int new_alloc = arma::sum(u > cum_probs);
    allocs_samp[i] = new_alloc + 1;
    n_curr[old_alloc - 1]--;
    n_curr[allocs_samp[i] - 1]++;
  }
  
  return List::create(Named("allocs_samp") = allocs_samp,
                      Named("n_curr") = n_curr);
}

// Multivariate functions:

// Function to compute Mahalanobis distance
arma::vec Mahalanobis(arma::mat const &x, 
                      arma::vec const &center, 
                      arma::mat const &cov) {
  arma::mat x_cen = x.t();
  x_cen.each_col() -= center;
  arma::solve(x_cen, arma::trimatl(arma::chol(cov).t()), x_cen);
  x_cen.for_each([](arma::mat::elem_type& val) { val = val * val; });
  return arma::sum(x_cen, 0).t();
}

// Function to compute multivariate normal density
// [[Rcpp::export]]
arma::vec dmvnorm_arma_mv(arma::mat const &x, 
                       arma::vec const &mean, 
                       arma::mat const &sigma, 
                       bool const logd = false) { 
  arma::vec const distval = Mahalanobis(x,  mean, sigma);
  double const logdet = arma::sum(arma::log(arma::eig_sym(sigma)));
  arma::vec const logretval = -( (x.n_cols * log2pi + logdet + distval)/2 );
  
  if (logd)
    return logretval;
  return arma::exp(logretval);
}





// [[Rcpp::export]]
List sample_allocs_mv(int N, int M_samp, int M_a_samp, int M_na_samp, 
                   arma::vec allocs_samp, arma::vec weights_samp, 
                   arma::mat y, arma::mat mus_a_samp, 
                   arma::cube Sigmas_a_samp, arma::mat mus_na_samp, 
                   arma::cube Sigmas_na_samp, arma::vec n_curr) {
  
  arma::vec alloc_rates(M_samp);
  for (int i = 0; i < N; ++i) {
    int old_alloc = allocs_samp[i];
    for (int c = 0; c < M_a_samp; ++c) {
      alloc_rates[c] = weights_samp[c] * dmvnorm_arma_mv(y.row(i), mus_a_samp.col(c), Sigmas_a_samp.slice(c))(0);
    }
    if (M_na_samp != 0) {
      for (int c = 0; c < M_na_samp; ++c) {
        alloc_rates[M_a_samp + c] = weights_samp[c] * dmvnorm_arma_mv(y.row(i), mus_na_samp.col(c), Sigmas_na_samp.slice(c))(0);
      }
    }
    alloc_rates /= arma::sum(alloc_rates);
    arma::vec cum_probs = arma::cumsum(alloc_rates);
    double u = R::runif(0, 1);
    int new_alloc = arma::sum(u > cum_probs);
    allocs_samp[i] = new_alloc + 1;
    n_curr[old_alloc - 1]--;
    n_curr[allocs_samp[i] - 1]++;
  }
  
  return List::create(Named("allocs_samp") = allocs_samp,
                      Named("n_curr") = n_curr);
}



// [[Rcpp::export]]
arma::vec post_density(arma::mat const &y_grid, 
                       arma::vec const &weights_samp, 
                       arma::mat const &mus_samp, 
                       arma::cube const &Sigmas_samp) {
  int n = y_grid.n_rows;
  int M_samp = weights_samp.n_elem;
  arma::vec post_dens(n, arma::fill::zeros);
  
  for (int i = 0; i < n; ++i) {
    for (int c = 0; c < M_samp; ++c) {
      // post_dens[i] += weights_samp[c] * dmvnorm_arma(y_grid.row(i), mus_samp.col(c), Sigmas_samp.slice(c));
      
      arma::vec density = weights_samp(c) * dmvnorm_arma_mv(y_grid.row(i), mus_samp.col(c), Sigmas_samp.slice(c));
      post_dens(i) += density(0);  // Sum over the first element (scalar result)
      
    }
  }
  
  return post_dens;
  
}


