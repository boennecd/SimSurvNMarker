// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector glq(
    Rcpp::NumericVector const lb, Rcpp::NumericVector const ub,
    Rcpp::NumericVector const nodes, Rcpp::NumericVector const weights,
    Rcpp::Function f){
  using Rcpp::NumericVector;

  size_t const n_out = lb.length(),
             n_nodes = weights.size();
  if(ub.size() != lb.size())
    throw std::invalid_argument("lb length != ub length");
  if(nodes.size() != weights.size())
    throw std::invalid_argument("nodes length != weights length");

  NumericVector out(n_out),
               node(1L);
  double &nv = node[0L];
  for(size_t i = 0; i < n_out; ++i){
    double &o_i = out[i];
    double const d1 = (ub[i] - lb[i]) / 2.,
                 d2 = (ub[i] + lb[i]) / 2.;

    o_i = 0.;
    double const *n = &nodes[0L],
                 *w = &weights[0L];
    for(unsigned j = 0; j < n_nodes; ++j, ++n, ++w){
      nv = d1 * *n + d2;
      o_i += *w * NumericVector(f(node))[0L];
    }

    o_i *= d1;
  }

  return out;
}
