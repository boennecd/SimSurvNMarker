
Show configurations

``` r
dput(alpha)
#> c(0.23, -0.07, -0.15)
dput(omega)
#> c(-1.5, -1.35, -1.9, -4.8, -0.17)
dput(delta)
#> c(-0.1, -0.08, 0.03, -0.21)
dput(gamma)
#> NULL
dput(B)
#> structure(c(-0.86, -0.28, -0.58, 0.04, 0.18, -0.12, 0.99, -0.23, 
#> -0.48, -0.24, 0.33, 0.06, 0.3, -0.57, -0.25), .Dim = c(5L, 3L
#> ))
dput(sig) # Sigma
#> structure(c(0.02, 0, 0, 0, 0.15, 0, 0, 0, 0.03), .Dim = c(3L, 
#> 3L))
dput(Psi)
#> structure(c(2.54, -0.64, 0.95, -0.35, -0.73, 0.79, -0.64, 1.91, 
#> -0.66, 0.26, -0.78, 1.08, 0.95, -0.66, 4.13, 1.99, 0.09, -0.94, 
#> -0.35, 0.26, 1.99, 3.86, 0.76, -1.34, -0.73, -0.78, 0.09, 0.76, 
#> 1.99, -1, 0.79, 1.08, -0.94, -1.34, -1, 3.06), .Dim = c(6L, 6L
#> ))
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

<img src="fig/no-gamma-plot_wo_marker-1.png" width="100%" />

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

<img src="fig/no-gamma-plot_wo_marker-2.png" width="100%" />

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

<img src="fig/no-gamma-show_sim_marker-1.png" width="100%" />

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

<img src="fig/no-gamma-show_draw_surv_curves-1.png" width="100%" />

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
#>   0.770   0.036   0.804
```

Show stats

``` r
# survival data
head(dat$survival_data)
#>   Z1 Z2 Z3 Z4 left_trunc    y event id
#> 1  0  1  1  0      1.017 3.84  TRUE  1
#> 2  1  1  0  0      2.725 6.92 FALSE  2
#> 3  1  0  1  0      1.076 7.00 FALSE  3
#> 4  1  0  0  0      2.163 9.01 FALSE  4
#> 5  0  1  1  0      0.280 5.41 FALSE  5
#> 6  1  1  1  0      0.135 7.24 FALSE  6

# marker data
head(dat$marker_data, 10)
#>    obs_time       Y1      Y2      Y3 id
#> 1      3.08  0.22624 -0.0776  0.0971  1
#> 2      3.27 -1.47036  0.3358  1.0065  2
#> 3      5.07 -1.55685 -0.4044  0.6314  2
#> 4      6.81 -1.78375 -1.6058  0.7136  2
#> 5      2.22 -0.01439 -0.3541 -0.7333  3
#> 6      2.30 -0.00835  0.4943 -0.5462  3
#> 7      3.34  0.15898  0.4300  0.2036  3
#> 8      4.20 -0.19284  0.6713  0.4609  3
#> 9      4.37  0.05193  0.9941  0.0635  3
#> 10     4.64  0.18234  0.2669  0.4269  3

# rate of observed events
mean(dat$survival_data$event) 
#> [1] 0.451

# mean event time
mean(subset(dat$survival_data, event)$y)
#> [1] 4.3

# quantiles of the event time
quantile(subset(dat$survival_data, event)$y)
#>     0%    25%    50%    75%   100% 
#> 0.0398 2.3445 4.1200 6.1927 9.9280

# fraction of observed markers per individual
NROW(dat$marker_data) / NROW(dat$survival_data)
#> [1] 5.83
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
#>      [,1] [,2] [,3]
#> 
#> $B
#>         [,1]   [,2]    [,3]
#> [1,] -0.8450 -0.131  0.3363
#> [2,] -0.2777  1.000  0.0347
#> [3,] -0.5814 -0.276  0.3053
#> [4,]  0.0368 -0.394 -0.6100
#> [5,]  0.1677 -0.148 -0.2311
#> 
#> $Psi
#>        [,1]   [,2]    [,3]   [,4]    [,5]   [,6]
#> [1,]  2.449 -0.675  0.8594 -0.390 -0.6532  0.795
#> [2,] -0.675  1.844 -0.6235  0.359 -0.8314  1.082
#> [3,]  0.859 -0.624  4.1381  1.975  0.0257 -0.971
#> [4,] -0.390  0.359  1.9753  4.033  0.6385 -1.263
#> [5,] -0.653 -0.831  0.0257  0.638  1.8987 -0.942
#> [6,]  0.795  1.082 -0.9706 -1.263 -0.9416  3.000
#> 
#> $Sigma
#>        [,1]   [,2]   [,3]
#> [1,] 0.0665 0.0000 0.0000
#> [2,] 0.0000 0.0665 0.0000
#> [3,] 0.0000 0.0000 0.0665
```

Compare with the true values

``` r
gamma
#> NULL
B
#>       [,1]  [,2]  [,3]
#> [1,] -0.86 -0.12  0.33
#> [2,] -0.28  0.99  0.06
#> [3,] -0.58 -0.23  0.30
#> [4,]  0.04 -0.48 -0.57
#> [5,]  0.18 -0.24 -0.25
Psi
#>       [,1]  [,2]  [,3]  [,4]  [,5]  [,6]
#> [1,]  2.54 -0.64  0.95 -0.35 -0.73  0.79
#> [2,] -0.64  1.91 -0.66  0.26 -0.78  1.08
#> [3,]  0.95 -0.66  4.13  1.99  0.09 -0.94
#> [4,] -0.35  0.26  1.99  3.86  0.76 -1.34
#> [5,] -0.73 -0.78  0.09  0.76  1.99 -1.00
#> [6,]  0.79  1.08 -0.94 -1.34 -1.00  3.06
sig
#>      [,1] [,2] [,3]
#> [1,] 0.02 0.00 0.00
#> [2,] 0.00 0.15 0.00
#> [3,] 0.00 0.00 0.03
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
  
  sformula <- Surv(left_trunc, y, ev) ~ 1
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
#>   n= 11651, number of events= 901 
#> 
#>       coef exp(coef) se(coef)     z Pr(>|z|)    
#> Z1 -0.1951    0.8228   0.0670 -2.91   0.0036 ** 
#> Z2 -0.0560    0.9456   0.0668 -0.84   0.4018    
#> Z3  0.1489    1.1606   0.0667  2.23   0.0256 *  
#> Z4 -0.3282    0.7203   0.0671 -4.89  1.0e-06 ***
#> Y1  0.2678    1.3071   0.0357  7.50  6.3e-14 ***
#> Y2 -0.0643    0.9378   0.0271 -2.37   0.0179 *  
#> Y3 -0.1499    0.8608   0.0345 -4.35  1.4e-05 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#>    exp(coef) exp(-coef) lower .95 upper .95
#> Z1     0.823      1.215     0.722     0.938
#> Z2     0.946      1.058     0.830     1.078
#> Z3     1.161      0.862     1.018     1.323
#> Z4     0.720      1.388     0.631     0.822
#> Y1     1.307      0.765     1.219     1.402
#> Y2     0.938      1.066     0.889     0.989
#> Y3     0.861      1.162     0.804     0.921
#> 
#> Concordance= 0.587  (se = 0.01 )
#> Likelihood ratio test= 114  on 7 df,   p=<2e-16
#> Wald test            = 115  on 7 df,   p=<2e-16
#> Score (logrank) test = 115  on 7 df,   p=<2e-16
```

Compare with the true value

``` r
delta
#> [1] -0.10 -0.08  0.03 -0.21
alpha
#> [1]  0.23 -0.07 -0.15
```
