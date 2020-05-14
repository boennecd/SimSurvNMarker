#define ARMA_NO_DEBUG
#include <RcppArmadillo.h>

// [[Rcpp::export(name = "eval_marker_cpp", rng = false)]]
arma::vec eval_marker(arma::mat const &B, arma::vec const &m){
  arma::vec out(B.n_cols, arma::fill::zeros);
  double const *b = B.begin();
  for(auto o = out.begin(); o != out.end(); ++o)
    for(auto mi = m.begin(); mi != m.end(); ++mi)
      *o += *mi * *b++;

  return out;
}
