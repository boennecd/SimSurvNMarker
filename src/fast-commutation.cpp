#include <Rcpp.h>

// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix get_commutation(unsigned const m){
  unsigned const mm = m * m, 
                mmm = mm * m, 
             mmm_p1 = mmm + 1L, 
              mm_pm = mm + m;
  Rcpp::NumericMatrix out(mm, mm);
  double *o = &out[0];
  unsigned inc_i(0L);
  
  for(unsigned i = 0; i < m; ++i, inc_i += m){
    double *o1 = o + inc_i + i * mm, 
           *o2 = o + i     + inc_i * mm;
    for(unsigned j = 0; j < i; ++j, o1 += mmm_p1, o2 += mm_pm){
      *o1 = 1.;
      *o2 = 1.;
    }
    *o1 += 1.;
  }
  return out;
}

/*** R
options(digits = 3)

library(matrixcalc)
get_commutation_R <- function(m){
  out <- matrix(0., nrow = m * m, ncol = m * m)
  for(i in 1:m)
    for(j in 1:m)
      out[(i - 1L) * m + j, (j - 1L) * m + i] <- 1.

  return(out)
}

for(i in 2:10){
  stopifnot(all.equal(commutation.matrix(i), get_commutation_R(i)))
  stopifnot(all.equal(commutation.matrix(i), get_commutation  (i)))
}

library(microbenchmark)
microbenchmark(
  matrixcalc = commutation.matrix(4L),
  R          = get_commutation_R (4L),
  cpp        = get_commutation   (4L),
  times = 10)
#R> Unit: nanoseconds
#R>        expr    min     lq   mean median     uq    max neval
#R>  matrixcalc 430595 432219 446759 440879 446573 518064    10
#R>           R   4590   4662   8042   5743   6627  29880    10
#R>         cpp    670    892   1902   1486   1791   6625    10

microbenchmark(
  matrixcalc = commutation.matrix(20L),
  R          = get_commutation_R (20L),
  cpp        = get_commutation   (20L),
  times = 10)
#R> Unit: microseconds
#R>        expr      min       lq   mean median       uq    max neval
#R>  matrixcalc 565638.8 572990.5 578593 578227 583098.4 592623    10
#R>           R    164.5    170.3    180    175    185.4    220    10
#R>         cpp     43.9     45.5    189     48     54.7   1450    10
*/
