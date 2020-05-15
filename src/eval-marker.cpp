#include <Rcpp.h>

// [[Rcpp::export(name = "eval_marker_cpp", rng = false)]]
Rcpp::NumericVector eval_marker(SEXP B, SEXP m){
  if (!Rf_isMatrix(B))
    throw std::invalid_argument("eval_marker: B must be a matrix");
  if (!Rf_isVector(m))
    throw std::invalid_argument("eval_marker: m must be a vector");
  size_t const nr = Rf_nrows(B),
               nc = Rf_ncols(B),
               nm = XLENGTH(m);
  if(nr != nm)
    throw std::invalid_argument("eval_marker: dims do not match");

  Rcpp::NumericVector out(nc);
  double const *b = REAL(B),
         *m_start = REAL(m),
         *m_end   = m_start + nm;
  for(double *o = out.begin(); o != out.end(); ++o)
    for(auto mi = m_start; mi != m_end; ++mi)
      *o += *mi * *b++;

  return out;
}
