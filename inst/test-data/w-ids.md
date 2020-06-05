
Assign parameters

``` r
alpha <- c(-0.42, 0.23)
omega <- c(-2.28, -1.69)
delta <- c(0.03, 0.23, 0.11, 0.23)
gamma <- structure(c(0.06, 0.2, -0.18, -0.27, 0.23, 0.16), .Dim = 3:2)
B <- structure(c(-0.68, 0.99, 0.06, 0.03, -0.97, 0.39), .Dim = 3:2)
sig <- structure(c(0.01, 0, 0, 0.04), .Dim = c(2L, 2L))
Psi <- structure(c(1.97, -0.24, -0.63, 0.29, -0.24, 0.35, -0.04, 0.16, 
                   -0.63, -0.04, 3.22, -1.39, 0.29, 0.16, -1.39, 3.1),
                 .Dim = c(4L, 4L))
n_obs <- 2000L

n_y <- length(alpha)
d_m <- NROW(Psi) / n_y
d_g <- NROW(B)
d_b <- length(omega)
d_z <- length(delta)
d_x <- NROW(gamma)
```

Define sampling functions

``` r
r_n_marker <- function(id)
  as.integer((id %% 4) * 2L + 1L)
r_obs_time <- function(id, n_markes)
  sort(runif(n_markes, 0, 10))
r_z <- function(id)
  as.numeric((id %% 5) == 1:4)
r_x <- function(id)
  as.numeric((id %% 4) == 1:3)
r_left_trunc <- function(id)
   id / 1000
r_right_cens <- function(id)
  10 - id / 1000
```

Get splines

``` r
b_ks <- seq(log(1e-1), log(10), length.out = d_b)
m_ks <- seq(       0 ,     10 , length.out = d_m)
g_ks <- seq(       0 ,     10 , length.out = d_g)

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

<img src="fig/w-ids-plot_wo_marker-1.png" width="100%" />

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

<img src="fig/w-ids-plot_wo_marker-2.png" width="100%" />

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

<img src="fig/w-ids-show_sim_marker-1.png" width="100%" />

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

<img src="fig/w-ids-show_draw_surv_curves-1.png" width="100%" />

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
#>   1.504   0.036   1.540
```

Show stats

``` r
# survival data
head(dat$survival_data, 10)
#>    Z1 Z2 Z3 Z4 left_trunc      y event id
#> 1   1  0  0  0      0.001 2.8080  TRUE  1
#> 2   0  1  0  0      0.002 1.1583  TRUE  2
#> 3   0  0  1  0      0.003 3.5103  TRUE  3
#> 4   0  0  0  1      0.004 0.0139  TRUE  4
#> 5   0  0  0  0      0.005 8.9899  TRUE  5
#> 6   1  0  0  0      0.006 0.7387  TRUE  6
#> 7   0  1  0  0      0.007 9.9930 FALSE  7
#> 8   0  0  1  0      0.008 9.9920 FALSE  8
#> 9   0  0  0  1      0.009 0.4170  TRUE  9
#> 10  0  0  0  0      0.010 8.4682  TRUE 10
tail(dat$survival_data, 10)
#>      Z1 Z2 Z3 Z4 left_trunc    y event   id
#> 1991  1  0  0  0       1.99 8.01 FALSE 1991
#> 1992  0  1  0  0       1.99 7.18  TRUE 1992
#> 1993  0  0  1  0       1.99 8.01 FALSE 1993
#> 1994  0  0  0  1       1.99 8.01 FALSE 1994
#> 1995  0  0  0  0       2.00 8.00 FALSE 1995
#> 1996  1  0  0  0       2.00 8.00 FALSE 1996
#> 1997  0  1  0  0       2.00 2.12  TRUE 1997
#> 1998  0  0  1  0       2.00 2.41  TRUE 1998
#> 1999  0  0  0  1       2.00 4.94  TRUE 1999
#> 2000  0  0  0  0       2.00 2.44  TRUE 2000

# marker data
head(dat$marker_data, 10)
#>    obs_time       Y1     Y2 X1 X2 X3 id
#> 1     0.618  0.24181 -3.324  1  0  0  1
#> 2     1.766  0.05588 -2.393  1  0  0  1
#> 3     2.060 -0.00941 -2.473  1  0  0  1
#> 4     1.079  0.96718  1.997  0  1  0  2
#> 5     0.842  0.44840 -1.639  0  0  1  3
#> 6     3.338  0.29634 -1.025  0  0  1  3
#> 7     3.391  0.25971 -0.797  0  0  1  3
#> 8     3.467  0.16297 -0.698  0  0  1  3
#> 9     0.012  0.18357 -3.618  0  0  0  4
#> 10    3.742  1.00674 -1.947  1  0  0  5
tail(dat$marker_data, 10)
#>      obs_time     Y1     Y2 X1 X2 X3   id
#> 4187     2.03  1.079 -1.841  0  0  1 1995
#> 4188     2.48  0.862 -1.387  0  0  1 1995
#> 4189     3.25  0.625 -1.575  0  0  1 1995
#> 4190     3.62  0.371 -1.501  0  0  1 1995
#> 4191     4.68  0.264 -1.074  0  0  1 1995
#> 4192     3.29  0.420 -4.769  0  0  0 1996
#> 4193     2.01  1.425 -1.724  1  0  0 1997
#> 4194     2.08  0.245 -0.115  0  1  0 1998
#> 4195     2.78 -1.549 -0.831  0  0  1 1999
#> 4196     2.13  1.717  2.154  0  0  0 2000

# rate of observed events
mean(dat$survival_data$event) 
#> [1] 0.667

# mean event time
mean(subset(dat$survival_data, event)$y)
#> [1] 3.72

# quantiles of the event time
quantile(subset(dat$survival_data, event)$y)
#>     0%    25%    50%    75%   100% 
#> 0.0139 1.9506 3.2577 5.3060 9.5956

# fraction of observed markers per individual
NROW(dat$marker_data) / NROW(dat$survival_data)
#> [1] 2.1
```

We get a special structure for the covariates, left-truncation,  
right-censoring time, etc. because we use the `id` argument in the
simulation function as we see above.

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
#> boundary (singular) fit: see ?isSingular
#> $gamma
#>          [,1]   [,2]
#> [1,]  0.00112 -0.320
#> [2,]  0.18671  0.255
#> [3,] -0.15542  0.149
#> 
#> $B
#>         [,1]    [,2]
#> [1,] -0.6371 -0.0557
#> [2,]  1.3095 -1.3566
#> [3,]  0.0514  0.3628
#> 
#> $Psi
#>        [,1]    [,2]    [,3]    [,4]
#> [1,]  1.887 -0.2623 -0.4520  0.2434
#> [2,] -0.262  0.0415 -0.0581  0.0131
#> [3,] -0.452 -0.0581  3.0266 -1.1894
#> [4,]  0.243  0.0131 -1.1894  2.8835
#> 
#> $Sigma
#>        [,1]   [,2]
#> [1,] 0.0353 0.0000
#> [2,] 0.0000 0.0353
```

Compare with the true values

``` r
gamma
#>       [,1]  [,2]
#> [1,]  0.06 -0.27
#> [2,]  0.20  0.23
#> [3,] -0.18  0.16
B
#>       [,1]  [,2]
#> [1,] -0.68  0.03
#> [2,]  0.99 -0.97
#> [3,]  0.06  0.39
Psi
#>       [,1]  [,2]  [,3]  [,4]
#> [1,]  1.97 -0.24 -0.63  0.29
#> [2,] -0.24  0.35 -0.04  0.16
#> [3,] -0.63 -0.04  3.22 -1.39
#> [4,]  0.29  0.16 -1.39  3.10
sig
#>      [,1] [,2]
#> [1,] 0.01 0.00
#> [2,] 0.00 0.04
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
#>   n= 4196, number of events= 2148 
#> 
#>       coef exp(coef) se(coef)     z Pr(>|z|)    
#> Z1  0.0362    1.0368   0.0703  0.51  0.60665    
#> Z2  0.1428    1.1535   0.0691  2.07  0.03880 *  
#> Z3  0.1764    1.1929   0.0688  2.57  0.01029 *  
#> Z4  0.2299    1.2585   0.0687  3.35  0.00082 ***
#> Y1 -0.1225    0.8847   0.0298 -4.11  4.0e-05 ***
#> Y2  0.0898    1.0940   0.0202  4.45  8.5e-06 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#>    exp(coef) exp(-coef) lower .95 upper .95
#> Z1     1.037      0.964     0.903     1.190
#> Z2     1.153      0.867     1.007     1.321
#> Z3     1.193      0.838     1.043     1.365
#> Z4     1.258      0.795     1.100     1.440
#> Y1     0.885      1.130     0.835     0.938
#> Y2     1.094      0.914     1.052     1.138
#> 
#> Concordance= 0.545  (se = 0.007 )
#> Likelihood ratio test= 55.9  on 6 df,   p=3e-10
#> Wald test            = 55.5  on 6 df,   p=4e-10
#> Score (logrank) test = 55.5  on 6 df,   p=4e-10
```

Compare with the true value

``` r
delta
#> [1] 0.03 0.23 0.11 0.23
alpha
#> [1] -0.42  0.23
```
