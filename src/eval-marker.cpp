#include <Rcpp.h>

// [[Rcpp::export(name = "eval_marker_cpp", rng = false)]]
Rcpp::NumericMatrix eval_marker(SEXP B, SEXP m){
  if(__builtin_expect(!Rf_isMatrix(B), 0))
    throw std::invalid_argument("eval_marker: B must be a matrix");

  if(Rf_isMatrix(m)){
    size_t const nr = Rf_nrows(B),
                 nc = Rf_ncols(B),
          n_col_out = Rf_nrows(m),
                 nm = Rf_ncols(m);
    if(__builtin_expect(nr != nm, 0))
      throw std::invalid_argument("eval_marker: dims do not match");

    Rcpp::NumericMatrix out(nc, n_col_out);

    double * o = out.begin();
    double const * const m_start = REAL(m),
                 * const b_start = REAL(B);

    for(size_t i = 0; i < n_col_out; ++i){
      double const * const o_end = o + nc,
                   *           b = b_start;
      for(; o != o_end; ++o){
        double const *           mi = m_start + i,
                     * const mi_end = mi + n_col_out * nm;

        for(; mi != mi_end; mi += n_col_out, ++b)
          *o += *mi * *b;
      }
    }

    return out;

  } else if(Rf_isVector(m)){
    size_t const nr = Rf_nrows(B),
                 nc = Rf_ncols(B),
                 nm = XLENGTH(m);
    if(__builtin_expect(nr != nm, 0))
      throw std::invalid_argument("eval_marker: dims do not match");

    Rcpp::NumericMatrix out(nc, 1L);
    double const *b = REAL(B),
           *m_start = REAL(m),
           *m_end   = m_start + nm;
      for(double *o = out.begin(); o != out.end(); ++o)
        for(auto mi = m_start; mi != m_end; ++mi)
          *o += *mi * *b++;

      return out;
  }

  throw std::invalid_argument("eval_marker: m is not a vector or a matrix");
  return Rcpp::NumericMatrix();
}
