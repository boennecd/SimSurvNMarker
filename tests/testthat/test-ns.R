test_that("cpp ns functions gives the same as the R version", {
  skip_if_not_installed("splines")
  require(splines)

  b_ks <- seq(1, 10, length.out = 4)
  b_func_R <- function(x)
    ns(x, knots = b_ks[-c(1L, length(b_ks))],
       Boundary.knots = b_ks[ c(1L, length(b_ks))],
       intercept = TRUE)

  b_func <- get_ns_spline(b_ks, do_log = FALSE)

  xs <- seq(0, 20, length.out = 1000)
  expect_equal(c(b_func_R(xs)), drop(b_func(xs)),
               check.attributes = FALSE)

  b_ks <- log(b_ks)
  b_func_R <- function(x)
    ns(log(x), knots = b_ks[-c(1L, length(b_ks))],
       Boundary.knots = b_ks[ c(1L, length(b_ks))],
       intercept = TRUE)

  b_func <- get_ns_spline(b_ks, do_log = TRUE)

  xs <- seq(1e-4, 20, length.out = 1000)
  expect_equal(c(b_func_R(xs)), drop(b_func(xs)),
               check.attributes = FALSE)
})


