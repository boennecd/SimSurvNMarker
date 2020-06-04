
Show configurations

``` r
dput(alpha)
#> c(-0.4, 0.5)
dput(omega)
#> c(-1.68, -1.93, -0.62, -5.11, -0.26)
dput(delta)
#> c(-0.09, -0.22)
dput(gamma)
#> structure(c(-0.15, 0.19, 0.06, -0.23), .Dim = c(2L, 2L))
dput(B)
#> structure(c(-0.43, 0.33, 0.04, 0.47, -0.37, -0.88, 0.48, -0.16, 
#> -0.06, -0.14), .Dim = c(5L, 2L))
dput(sig) # Sigma
#> structure(c(0.02, 0, 0, 0.14), .Dim = c(2L, 2L))
dput(Psi)
#> structure(c(5.52, 1.09, 0.24, 0.28, 1.09, 3.52, 0.56, -0.37, 
#> 0.24, 0.56, 1.85, 0.5, 0.28, -0.37, 0.5, 1.96), .Dim = c(4L, 
#> 4L))
dput(n_obs)
#> 2000L
```

Define sampling functions

``` r
r_n_marker <- function(id)
  rpois(1, 10) + 1L
r_obs_time <- function(id, n_markes)
  sort(runif(n_markes, 0, 10))
r_z <- function(id)
  as.numeric(runif(d_z) > .5)
r_x <- function(id)
  as.numeric(runif(d_x) > .5)
r_left_trunc <- function(id)
   rbeta(1, 1, 2) * 3
r_right_cens <- function(id)
  rbeta(1, 2, 1) * 6 + 4
```

Get splines

``` r
b_func <- get_ns_spline(b_ks, do_log = TRUE)
m_func <- get_ns_spline(m_ks, do_log = FALSE)
g_func <- get_ns_spline(g_ks, do_log = FALSE)
```

Get the Gauss-Legendre quadrature nodes we need

``` r
gl_dat <- get_gl_rule(30L)
```

Plot baseline hazard and survival function without the marker

``` r
library(SimSurvNMarker)
```

``` r
# hazard function without marker
par(mar = c(5, 5, 1, 1))
plot(function(x) exp(drop(b_func(x) %*% omega)),
     xlim = c(1e-8, 10), ylim = c(0, .61), xlab = "Time",
     ylab = "Hazard (no marker)", xaxs = "i", bty = "l")
```

<img src="fig/w-all-plot_wo_marker-1.png" width="100%" />

``` r

# survival function without marker
plot(function(x) eval_surv_base_fun(x, omega = omega, b_func = b_func, 
                                    gl_dat = gl_dat, delta = NULL), 
     xlim = c(1e-4, 10),
     xlab = "Time", ylab = "Survival probability (no marker)", xaxs = "i",
     yaxs = "i", bty = "l", ylim = c(0, 1.01))
abline(h = .75, lty = 3)
abline(h = .25, lty = 3)
```

<img src="fig/w-all-plot_wo_marker-2.png" width="100%" />

Simulate a few markers as an example

``` r
set.seed(1)
show_mark_mean <- function(B, Psi, sigma, m_func, g_func){
  tis <- seq(0, 10, length.out = 100)
  Psi_chol <- chol(Psi)
  
  par_old <- par(no.readonly = TRUE)
  on.exit(par(par_old))
  par(mar = c(4, 3, 1, 1), mfcol = c(4, 4))
  
  sigma_chol <- chol(sigma)
  n_y <- NCOL(sigma_chol)
  for(i in 1:16){
    U <- draw_U(Psi_chol, n_y = n_y)
    y_non_rng <- eval_marker(tis, B = B, g_func = g_func, U = NULL, 
                             offset = NULL, m_func = m_func)
    y_rng     <- eval_marker(tis, B = B, g_func = g_func, U = U, 
                             offset = NULL, m_func = m_func)
    if(length(B) == 0L){
      y_non_rng <- y_rng
      y_non_rng[] <- 0.
    }
    
    if(!is.vector(y_non_rng)){
      y_non_rng <- t(y_non_rng)
      y_rng     <- t(y_rng)
    } else {
      y_non_rng <- as.matrix(y_non_rng)
      y_rng     <- as.matrix(y_rng)
    }

    sds <- sapply(tis, function(ti){
      M <- (diag(n_y) %x% m_func(ti))
      G <- (diag(n_y) %x% g_func(ti))
      sds <- sqrt(diag(tcrossprod(M %*% Psi, M)))
      if(length(B) > 0)
        cbind(drop(G %*% c(B)) - 1.96 * sds,
              drop(G %*% c(B)) + 1.96 * sds)
      else 
        cbind(- 1.96 * sds, 1.96 * sds)
    }, simplify = "array")
    lbs <- sds[, 1, ]
    ubs <- sds[, 2, ]
    if(!is.vector(lbs)){
      lbs <- t(lbs)
      ubs <- t(ubs)
    } else {
      lbs <- as.matrix(lbs)
      ubs <- as.matrix(ubs)
    }

    y_obs <- sim_marker(B = B, U = U, sigma_chol = sigma_chol, 
                        m_func = m_func, r_n_marker = r_n_marker, 
                        r_obs_time = r_obs_time, g_func = g_func, 
                        offset = NULL)

    matplot(tis, y_non_rng, type = "l", lty = 2, ylab = "", xlab = "Time",
            ylim = range(y_non_rng, y_rng, lbs, ubs, y_obs$y_obs))
    matplot(tis, y_rng    , type = "l", lty = 1, add = TRUE)
    matplot(y_obs$obs_time, y_obs$y_obs, type = "p", add = TRUE, pch = 3:4)

    for(i in 1:NCOL(y_non_rng)){
      rg <- col2rgb(i) / 255
      polygon(c(tis, rev(tis)), c(lbs[, i], rev(ubs[, i])), border = NA,
              col = rgb(rg[1], rg[2], rg[3], .1))
    }

  }
  invisible()
}
show_mark_mean(B = B, Psi = Psi, sigma = sig, m_func = m_func, 
               g_func = g_func)
```

<img src="fig/w-all-show_sim_marker-1.png" width="100%" />

Illustrate a few conditional hazard functions and survival functions

``` r
set.seed(1)
local({
  par_old <- par(no.readonly = TRUE)
  on.exit(par(par_old))
  par(mfcol = c(1, 2))

  # hazard functions
  tis <- seq(1e-4, 10, length.out = 50)
  n_y <- NCOL(sig)
  Us <- replicate(100, draw_U(chol(Psi), n_y = n_y), 
                  simplify = "array")

  hz <- apply(Us, 3L, function(U)
    vapply(tis, function(x)
      exp(drop(b_func(x) %*% omega +
                 alpha %*% eval_marker(ti = x, B = B, m_func = m_func, 
                                       g_func = g_func, U = U, 
                                       offset = NULL))),
      FUN.VALUE = numeric(1L)))

  matplot(tis, hz, lty = 1, type = "l", col = rgb(0, 0, 0, .1),
          xaxs = "i", bty = "l", yaxs = "i", ylim = c(0, 1),
          xlab = "time", ylab = "Hazard")

  # survival functions
  ys <- apply(Us, 3L, surv_func_joint,
              ti = tis, B = B, omega = omega, delta = NULL,
              alpha = alpha, b_func = b_func, m_func = m_func, 
              gl_dat = gl_dat, g_func = g_func, offset = NULL)

  matplot(tis, ys, lty = 1, type = "l", col = rgb(0, 0, 0, .1),
          xaxs = "i", bty = "l", yaxs = "i", ylim = c(0, 1),
          xlab = "time", ylab = "Survival probability")
  abline(h = .75, lty = 3)
  abline(h = .25, lty = 3)
})
```

<img src="fig/w-all-show_draw_surv_curves-1.png" width="100%" />

Simulate a data set

``` r
set.seed(1)
system.time(dat <- sim_joint_data_set(
  n_obs = n_obs, B = B, Psi = Psi, omega = omega, delta = delta, 
  alpha = alpha, sigma = sig, gamma = gamma, b_func = b_func, 
  m_func = m_func, g_func = g_func, gl_dat = gl_dat, r_z = r_z, 
  r_left_trunc = r_left_trunc, r_right_cens = r_right_cens, 
  r_n_marker = r_n_marker, r_x = r_x, r_obs_time = r_obs_time, y_max = 10))
#>    user  system elapsed 
#>   0.941   0.055   0.995
```

Show stats

``` r
# survival data
head(dat$survival_data)
#>   Z1 Z2 left_trunc    y event id
#> 1  0  0      1.735 4.62  TRUE  1
#> 2  0  0      0.338 3.74  TRUE  2
#> 3  1  0      1.320 3.58  TRUE  3
#> 4  0  1      2.150 2.31  TRUE  4
#> 5  1  1      0.727 6.32 FALSE  5
#> 6  1  0      0.627 2.27  TRUE  6

# marker data
head(dat$marker_data, 10)
#>    obs_time     Y1     Y2 X1 X2 id
#> 1      1.86 -1.611  0.989  1  1  1
#> 2      4.11 -0.522  0.717  1  1  1
#> 3      1.29  1.093 -0.390  1  0  2
#> 4      2.13  0.698 -1.006  1  0  2
#> 5      2.70  0.804 -0.209  1  0  2
#> 6      3.53  0.644 -0.609  1  0  2
#> 7      1.81  0.744 -0.619  0  0  3
#> 8      2.19  0.658 -0.548  0  0  3
#> 9      2.68  0.814 -0.114  0  0  3
#> 10     2.69  0.744 -0.359  0  0  3

# rate of observed events
mean(dat$survival_data$event) 
#> [1] 0.531

# mean event time
mean(subset(dat$survival_data, event)$y)
#> [1] 3.73

# quantiles of the event time
quantile(subset(dat$survival_data, event)$y)
#>     0%    25%    50%    75%   100% 
#> 0.0711 1.9541 3.3230 5.3027 9.8106

# fraction of observed markers per individual
NROW(dat$marker_data) / NROW(dat$survival_data)
#> [1] 5.41
```

Fit linear mixed model and see that we get estimates which are close to
the true values

``` r
library(lme4)
library(reshape2)
library(splines)
.GlobalEnv$ns_func <- function(x, knots){
  is_bk <- c(1L, length(knots))
  ns(x, knots = knots[-is_bk], Boundary.knots = knots[is_bk], 
     intercept = TRUE)
}

local({
  m_dat <- dat$marker_data
  
  Y_names <- paste0("Y", 1:n_y)
  id_vars <- c("id", "obs_time")
  if(d_x > 0)
    id_vars <- c(id_vars, paste0("X", seq_len(d_x)))
  
  lme_dat <- melt(m_dat, id.vars = id_vars, measure.vars = Y_names, 
                  variable.name = "XXTHEVARIABLEXX", 
                  value.name = "XXTHEVALUEXX")
  
  if(length(alpha) > 1){
    if(length(B) > 0L)
      frm <- substitute(
        XXTHEVALUEXX ~
          XXTHEVARIABLEXX : ns_func(ti, g_ks) - 1L +
          (XXTHEVARIABLEXX : ns_func(ti, m_ks) - 1L | i),
        list(ti = as.name("obs_time"), i = as.name("id"), 
             g_ks = as.name("g_ks"), m_ks = as.name("m_ks")))
    else 
      frm <- substitute(
        XXTHEVALUEXX ~
          (XXTHEVARIABLEXX : ns_func(ti, m_ks) - 1L | i),
        list(ti = as.name("obs_time"), i = as.name("id"), 
             m_ks = as.name("m_ks")))
    frm <- eval(frm)
    
    if(d_x > 0)
      for(i in rev(seq_len(d_x))){
        frm_call <- substitute(
          update(frm, . ~ XXTHEVARIABLEXX : x_var + .),
          list(x_var = as.name(paste0("X", i))))
        frm <- eval(frm_call)
      }
    
  } else {
    if(length(B) > 0L)
      frm <- substitute(
        XXTHEVALUEXX ~
          ns_func(ti, g_ks) - 1L +
          (ns_func(ti, m_ks) - 1L | i),
        list(ti = as.name("obs_time"), i = as.name("id"), 
             g_ks = as.name("g_ks"), m_ks = as.name("m_ks")))
    else 
      frm <- substitute(
        XXTHEVALUEXX ~
          (ns_func(ti, m_ks) - 1L | i),
        list(ti = as.name("obs_time"), i = as.name("id"), 
             m_ks = as.name("m_ks")))
    frm <- eval(frm)
    
    if(d_x > 0)
      for(i in rev(seq_len(d_x))){
        frm_call <- substitute(
          update(frm, . ~ x_var + .),
          list(x_var = as.name(paste0("X", i))))
        frm <- eval(frm_call)
      }
  }
        
  fit <- lmer(frm, lme_dat, control = lmerControl(
    check.conv.grad = .makeCC("ignore", tol = 1e-3, relTol = NULL)))
  
  gamma <- t(matrix(fixef(fit)[seq_len(d_x * n_y)], nr = n_y))
  
  B <- t(matrix(fixef(fit)[seq_len(d_g * n_y) + (d_x * n_y)], nr = n_y))
  vc <- VarCorr(fit)
  Psi <- vc$id
  attr(Psi, "correlation") <- attr(Psi, "stddev") <- NULL
  dimnames(Psi) <- NULL
  K <- SimSurvNMarker:::get_commutation(n_y, d_m)
  Psi <- tcrossprod(K %*% Psi, K)

  Sigma <- diag(attr(vc, "sc")^2, n_y)

  list(gamma = gamma, B = B, Psi = Psi, Sigma = Sigma)
})
#> $gamma
#>        [,1]    [,2]
#> [1,] -0.114  0.0524
#> [2,]  0.178 -0.2113
#> 
#> $B
#>         [,1]   [,2]
#> [1,] -0.3220 -0.940
#> [2,]  0.4245  0.436
#> [3,]  0.0848 -0.241
#> [4,]  0.6978 -0.218
#> [5,] -0.3817 -0.194
#> 
#> $Psi
#>       [,1]   [,2]  [,3]   [,4]
#> [1,] 5.464  1.045 0.418  0.376
#> [2,] 1.045  3.037 0.558 -0.310
#> [3,] 0.418  0.558 1.789  0.463
#> [4,] 0.376 -0.310 0.463  2.146
#> 
#> $Sigma
#>        [,1]   [,2]
#> [1,] 0.0792 0.0000
#> [2,] 0.0000 0.0792
```

Compare with the true values

``` r
gamma
#>       [,1]  [,2]
#> [1,] -0.15  0.06
#> [2,]  0.19 -0.23
B
#>       [,1]  [,2]
#> [1,] -0.43 -0.88
#> [2,]  0.33  0.48
#> [3,]  0.04 -0.16
#> [4,]  0.47 -0.06
#> [5,] -0.37 -0.14
Psi
#>      [,1]  [,2] [,3]  [,4]
#> [1,] 5.52  1.09 0.24  0.28
#> [2,] 1.09  3.52 0.56 -0.37
#> [3,] 0.24  0.56 1.85  0.50
#> [4,] 0.28 -0.37 0.50  1.96
sig
#>      [,1] [,2]
#> [1,] 0.02 0.00
#> [2,] 0.00 0.14
```

Fit Cox model with only the observed markers (likely biased)

``` r
local({
  library(survival)
  tdat <- tmerge(dat$survival_data, dat$survival_data, id = id, 
                 tstart = left_trunc, tstop = y, ev = event(y, event))
  
  for(i in seq_along(alpha)){
    new_call <- substitute(tmerge(
      tdat, dat$marker_data, id = id, tdc(obs_time, YVAR)),
      list(YVAR = as.name(paste0("Y", i))))
    names(new_call)[length(new_call)] <- paste0("Y", i)
    tdat <- eval(new_call)
  }
  tdat <- na.omit(tdat)
  
  sformula <- Surv(left_trunc, y, event) ~ 1
  for(i in seq_along(delta)){
    new_call <- substitute(update(sformula, . ~ . + XVAR), 
                           list(XVAR = as.name(paste0("Z", i))))
    sformula <- eval(new_call)
  }
  for(i in seq_along(alpha)){
    new_call <- substitute(update(sformula, . ~ . + XVAR), 
                           list(XVAR = as.name(paste0("Y", i))))
    sformula <- eval(new_call)
  }
  
  fit <- coxph(sformula, tdat)
  print(summary(fit))  
  invisible(fit)
})
#> Call:
#> coxph(formula = sformula, data = tdat)
#> 
#>   n= 10816, number of events= 3676 
#> 
#>       coef exp(coef) se(coef)      z Pr(>|z|)    
#> Z1  0.0145    1.0146   0.0331   0.44     0.66    
#> Z2 -0.2413    0.7856   0.0331  -7.29  3.1e-13 ***
#> Y1 -0.3755    0.6870   0.0132 -28.47  < 2e-16 ***
#> Y2  0.2272    1.2551   0.0188  12.08  < 2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#>    exp(coef) exp(-coef) lower .95 upper .95
#> Z1     1.015      0.986     0.951     1.083
#> Z2     0.786      1.273     0.736     0.838
#> Y1     0.687      1.456     0.669     0.705
#> Y2     1.255      0.797     1.210     1.302
#> 
#> Concordance= 0.644  (se = 0.005 )
#> Likelihood ratio test= 953  on 4 df,   p=<2e-16
#> Wald test            = 950  on 4 df,   p=<2e-16
#> Score (logrank) test = 951  on 4 df,   p=<2e-16
```

Compare with the true value

``` r
delta
#> [1] -0.09 -0.22
alpha
#> [1] -0.4  0.5
```
