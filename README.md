
# SimSurvNMarker: Simulate Survival time and Markers

[![R-CMD-check](https://github.com/boennecd/SimSurvNMarker/workflows/R-CMD-check/badge.svg)](https://github.com/boennecd/SimSurvNMarker/actions)
[![](https://www.r-pkg.org/badges/version/SimSurvNMarker)](https://CRAN.R-project.org/package=SimSurvNMarker)
[![CRAN RStudio mirror
downloads](http://cranlogs.r-pkg.org/badges/SimSurvNMarker)](https://CRAN.R-project.org/package=SimSurvNMarker)

The `SimSurvNMarker` package reasonably fast simulates from a joint
survival and marker model. The package uses a combination of
Gauss–Legendre quadrature and one dimensional root finding to simulate
the event times as described by Crowther and Lambert (2013).
Specifically, the package can simulate from the model

<!-- $$\begin{align*} -->

<!-- \vec Y_{ij} \mid \vec U_i = \vec u_i -->

<!--   &\sim N^{(r)}(\vec \mu_i(s_{ij}, \vec u_i), \Sigma) -->

<!--   \\ -->

<!-- \vec\mu_i(s, \vec u) &= -->

<!--   \Gamma^\top \vec x_i + B^\top\vec g(s) + U^\top\vec m(s) -->

<!--   \\ -->

<!-- &= \left(I \otimes \vec x_i^\top\right)\text{vec}\Gamma -->

<!--      + \left(I \otimes \vec g(s)^\top\right)\text{vec} B -->

<!--      + \left(I \otimes \vec m(s)^\top\right) \vec u -->

<!--   \\ -->

<!-- \vec U_i &\sim N^{(K)}(\vec 0, \Psi) -->

<!--   \\ -->

<!-- h(t\mid \vec u) &= \exp\left( -->

<!--   \vec\omega^\top\vec b(t) + -->

<!--   \vec z_i^\top\vec\delta + -->

<!--   \vec\alpha^\top\vec\mu_i(t, \vec u) -->

<!--   \right) -->

<!--   \\ -->

<!-- &= \exp\Bigg( -->

<!--   \vec\omega^\top\vec b(t) + -->

<!--   \vec z_i^\top\vec\delta -->

<!--   + \vec 1^\top\left( -->

<!--   \text{diag}(\vec \alpha) \otimes \vec x_i^\top\right)\text{vec}\Gamma -->

<!--   + \vec 1^\top\left( -->

<!--   \text{diag}(\vec \alpha) \otimes \vec g(t)^\top\right)\text{vec} B \\ -->

<!-- &\hspace{50pt}+ \vec 1^\top\left( -->

<!--   \text{diag}(\vec \alpha) \otimes \vec m(t)^\top\right)\vec u -->

<!--   \Bigg) -->

<!-- \end{align*}$$ -->

  
![\\begin{align\*}\\vec Y\_{ij}\\mid\\vec U\_i =\\vec u\_i &\\sim
N^{(r)}(\\vec \\mu\_i(s\_{ij}, \\vec u\_i), \\Sigma)\\\\\\vec\\mu\_i(s,
\\vec u)&=\\Gamma^\\top\\vec x\_i+B^\\top\\vec g(s)+U^\\top\\vec
m(s)\\\\&=\\left(I\\otimes\\vec
x\_i^\\top\\right)\\text{vec}\\Gamma+\\left(I\\otimes\\vec
g(s)^\\top\\right)\\text{vec}B+\\left(I\\otimes\\vec
m(s)^\\top\\right)\\vec u\\\\\\vec U\_i &\\sim
N^{(K)}(\\vec 0,\\Psi)\\end{align\*}](https://render.githubusercontent.com/render/math?math=%5Cbegin%7Balign%2A%7D%5Cvec%20Y_%7Bij%7D%5Cmid%5Cvec%20U_i%20%3D%5Cvec%20u_i%20%26%5Csim%20N%5E%7B%28r%29%7D%28%5Cvec%20%5Cmu_i%28s_%7Bij%7D%2C%20%5Cvec%20u_i%29%2C%20%5CSigma%29%5C%5C%5Cvec%5Cmu_i%28s%2C%20%5Cvec%20u%29%26%3D%5CGamma%5E%5Ctop%5Cvec%20x_i%2BB%5E%5Ctop%5Cvec%20g%28s%29%2BU%5E%5Ctop%5Cvec%20m%28s%29%5C%5C%26%3D%5Cleft%28I%5Cotimes%5Cvec%20x_i%5E%5Ctop%5Cright%29%5Ctext%7Bvec%7D%5CGamma%2B%5Cleft%28I%5Cotimes%5Cvec%20g%28s%29%5E%5Ctop%5Cright%29%5Ctext%7Bvec%7DB%2B%5Cleft%28I%5Cotimes%5Cvec%20m%28s%29%5E%5Ctop%5Cright%29%5Cvec%20u%5C%5C%5Cvec%20U_i%20%26%5Csim%20N%5E%7B%28K%29%7D%28%5Cvec%200%2C%5CPsi%29%5Cend%7Balign%2A%7D
"\\begin{align*}\\vec Y_{ij}\\mid\\vec U_i =\\vec u_i &\\sim N^{(r)}(\\vec \\mu_i(s_{ij}, \\vec u_i), \\Sigma)\\\\\\vec\\mu_i(s, \\vec u)&=\\Gamma^\\top\\vec x_i+B^\\top\\vec g(s)+U^\\top\\vec m(s)\\\\&=\\left(I\\otimes\\vec x_i^\\top\\right)\\text{vec}\\Gamma+\\left(I\\otimes\\vec g(s)^\\top\\right)\\text{vec}B+\\left(I\\otimes\\vec m(s)^\\top\\right)\\vec u\\\\\\vec U_i &\\sim N^{(K)}(\\vec 0,\\Psi)\\end{align*}")  

  
![\\begin{align\*}h(t\\mid\\vec u)&=\\exp\\left(\\vec\\omega^\\top\\vec
b(t)+\\vec z\_i^\\top\\vec\\delta+\\vec\\alpha^\\top\\vec\\mu\_i(t,
\\vec u)\\right)\\\\&=\\exp\\Bigg(\\vec\\omega^\\top\\vec b(t)+\\vec
z\_i^\\top\\vec\\delta+\\vec 1^\\top\\left(\\text{diag}(\\vec\\alpha)\\otimes\\vec
x\_i^\\top\\right)\\text{vec}\\Gamma+\\vec 1^\\top\\left(\\text{diag}(\\vec\\alpha)\\otimes\\vec
g(t)^\\top\\right)\\text{vec}
B\\\\&\\hspace{50pt}+\\vec 1^\\top\\left(\\text{diag}(\\vec\\alpha)\\otimes\\vec
m(t)^\\top\\right)\\vec
u\\Bigg)\\end{align\*}](https://render.githubusercontent.com/render/math?math=%5Cbegin%7Balign%2A%7Dh%28t%5Cmid%5Cvec%20u%29%26%3D%5Cexp%5Cleft%28%5Cvec%5Comega%5E%5Ctop%5Cvec%20b%28t%29%2B%5Cvec%20z_i%5E%5Ctop%5Cvec%5Cdelta%2B%5Cvec%5Calpha%5E%5Ctop%5Cvec%5Cmu_i%28t%2C%20%5Cvec%20u%29%5Cright%29%5C%5C%26%3D%5Cexp%5CBigg%28%5Cvec%5Comega%5E%5Ctop%5Cvec%20b%28t%29%2B%5Cvec%20z_i%5E%5Ctop%5Cvec%5Cdelta%2B%5Cvec%201%5E%5Ctop%5Cleft%28%5Ctext%7Bdiag%7D%28%5Cvec%5Calpha%29%5Cotimes%5Cvec%20x_i%5E%5Ctop%5Cright%29%5Ctext%7Bvec%7D%5CGamma%2B%5Cvec%201%5E%5Ctop%5Cleft%28%5Ctext%7Bdiag%7D%28%5Cvec%5Calpha%29%5Cotimes%5Cvec%20g%28t%29%5E%5Ctop%5Cright%29%5Ctext%7Bvec%7D%20B%5C%5C%26%5Chspace%7B50pt%7D%2B%5Cvec%201%5E%5Ctop%5Cleft%28%5Ctext%7Bdiag%7D%28%5Cvec%5Calpha%29%5Cotimes%5Cvec%20m%28t%29%5E%5Ctop%5Cright%29%5Cvec%20u%5CBigg%29%5Cend%7Balign%2A%7D
"\\begin{align*}h(t\\mid\\vec u)&=\\exp\\left(\\vec\\omega^\\top\\vec b(t)+\\vec z_i^\\top\\vec\\delta+\\vec\\alpha^\\top\\vec\\mu_i(t, \\vec u)\\right)\\\\&=\\exp\\Bigg(\\vec\\omega^\\top\\vec b(t)+\\vec z_i^\\top\\vec\\delta+\\vec 1^\\top\\left(\\text{diag}(\\vec\\alpha)\\otimes\\vec x_i^\\top\\right)\\text{vec}\\Gamma+\\vec 1^\\top\\left(\\text{diag}(\\vec\\alpha)\\otimes\\vec g(t)^\\top\\right)\\text{vec} B\\\\&\\hspace{50pt}+\\vec 1^\\top\\left(\\text{diag}(\\vec\\alpha)\\otimes\\vec m(t)^\\top\\right)\\vec u\\Bigg)\\end{align*}")  

where ![\\vec Y\_{ij}\\in\\mathbb
R^{n\_y}](https://render.githubusercontent.com/render/math?math=%5Cvec%20Y_%7Bij%7D%5Cin%5Cmathbb%20R%5E%7Bn_y%7D
"\\vec Y_{ij}\\in\\mathbb R^{n_y}") is individual
![i](https://render.githubusercontent.com/render/math?math=i "i")’s
![j](https://render.githubusercontent.com/render/math?math=j "j")th
observed marker at time
![s\_{ij}](https://render.githubusercontent.com/render/math?math=s_%7Bij%7D
"s_{ij}"), ![U\_i\\in\\mathbb
R^K](https://render.githubusercontent.com/render/math?math=U_i%5Cin%5Cmathbb%20R%5EK
"U_i\\in\\mathbb R^K") is individual
![i](https://render.githubusercontent.com/render/math?math=i "i")’s
random effect, and
![h](https://render.githubusercontent.com/render/math?math=h "h") is the
instantaneous hazard rate for the time-to-event outcome.
![\\vec\\alpha](https://render.githubusercontent.com/render/math?math=%5Cvec%5Calpha
"\\vec\\alpha") is the so-called association parameter. It shows the
strength of the relation between the latent mean function,
![\\vec\\mu\_i(t,\\vec
u)](https://render.githubusercontent.com/render/math?math=%5Cvec%5Cmu_i%28t%2C%5Cvec%20u%29
"\\vec\\mu_i(t,\\vec u)"), and the log of the instantaneous rate,
![h(t\\mid \\vec
u)](https://render.githubusercontent.com/render/math?math=h%28t%5Cmid%20%5Cvec%20u%29
"h(t\\mid \\vec u)"). ![\\vec
m(t)](https://render.githubusercontent.com/render/math?math=%5Cvec%20m%28t%29
"\\vec m(t)"), ![\\vec
g(t)](https://render.githubusercontent.com/render/math?math=%5Cvec%20g%28t%29
"\\vec g(t)") and ![\\vec
b(t)](https://render.githubusercontent.com/render/math?math=%5Cvec%20b%28t%29
"\\vec b(t)") are basis expansions of time. As an example, these can be
a polynomial, a B-spline, or a natural cubic spline. The expansion for
the baseline hazard, ![\\vec
b(t)](https://render.githubusercontent.com/render/math?math=%5Cvec%20b%28t%29
"\\vec b(t)"), is typically made on ![\\log
t](https://render.githubusercontent.com/render/math?math=%5Clog%20t
"\\log t") instead of
![t](https://render.githubusercontent.com/render/math?math=t "t"). One
reason is that the model reduces to a Weibull distribution when a first
polynomial is used and ![\\vec\\alpha =
\\vec 0](https://render.githubusercontent.com/render/math?math=%5Cvec%5Calpha%20%3D%20%5Cvec%200
"\\vec\\alpha = \\vec 0"). ![\\vec
x\_i](https://render.githubusercontent.com/render/math?math=%5Cvec%20x_i
"\\vec x_i") and ![\\vec
z\_i](https://render.githubusercontent.com/render/math?math=%5Cvec%20z_i
"\\vec z_i") are individual specific known time-invariant covariates.

We provide an example of how to use the package here, in
[inst/test-data](https://github.com/boennecd/SimSurvNMarker/tree/master/inst/test-data)
directory on Github, and at
[rpubs.com/boennecd/SimSurvNMarker-ex](https://rpubs.com/boennecd/SimSurvNMarker-ex)
where we show that we simulate from the correct model by using a
simulation study. The former examples includes

  - using [polynomials](#example-polynomial) as the basis functions.
  - using [polynomials with derivatives in the log
    hazard](#using-derivatives) as the basis functions.
  - using [natural cubic splines](#example-natural-cubic-splines) as the
    basis functions.

The purpose/goal of this package is to

  - allow the user to specify all components in R including the basis
    functions.
  - simulate reasonably fast from the model given this requirement.

## Installation

The package can be installed from Github using the`remotes` package:

``` r
stopifnot(require(remotes)) # needs the remotes package
install_github("boennecd/SimSurvNMarker")
```

It can also be installed from CRAN using `install.packages`:

``` r
install.packages("SimSurvNMarker")
```

## Example with Polynomials

We start with an example where we use polynomials as the basis
functions. First, we assign the polynomial functions we will use.

``` r
b_func <- function(x){
  x <- x - 1
  cbind(x^3, x^2, x)
}
g_func <- function(x){
  x <- x - 1
  cbind(x^3, x^2, x)
}
m_func <- function(x){
  x <- x - 1
  cbind(x^2, x, 1)
}
```

We use a third order polynomial for the two fixed terms, ![\\vec
b](https://render.githubusercontent.com/render/math?math=%5Cvec%20b
"\\vec b") and ![\\vec
g](https://render.githubusercontent.com/render/math?math=%5Cvec%20g
"\\vec g"), and a second order random polynomial for the random term,
![U\_i^\\top \\vec
m(s)](https://render.githubusercontent.com/render/math?math=U_i%5E%5Ctop%20%5Cvec%20m%28s%29
"U_i^\\top \\vec m(s)"), in the latent mean function. We choose the
following parameters for the baseline hazard.

``` r
library(SimSurvNMarker)
```

``` r
omega <- c(1.4, -1.2, -2.1)
delta <- -.5 # the intercept

# quadrature nodes
gl_dat <- get_gl_rule(30L)

# hazard function without marker
par(mar = c(5, 5, 1, 1), mfcol = c(1, 2))
plot(function(x) exp(drop(b_func(x) %*% omega) + delta),
     xlim = c(0, 2), ylim = c(0, 1.5), xlab = "Time",
     ylab = "Hazard (no marker)", xaxs = "i",  yaxs = "i", bty = "l")
grid()

# survival function without marker
plot(function(x) eval_surv_base_fun(x, omega = omega, b_func = b_func, 
                                    gl_dat = gl_dat, delta = delta), 
     xlim = c(1e-4, 2.5),
     xlab = "Time", ylab = "Survival probability (no marker)", xaxs = "i",
     yaxs = "i", bty = "l", ylim = c(0, 1.01))
grid()
```

<img src="man/figures/README-poly_base_haz-1.png" width="100%" />

Then we set the following parameters for the random effect, ![\\vec
U\_i](https://render.githubusercontent.com/render/math?math=%5Cvec%20U_i
"\\vec U_i"), and the parameters for the marker process. We also
simulate a number of latent markers’ mean curves and observed marker
values and plot the result. The dashed curve is the mean,
![\\vec\\mu\_i(s,
\\vec 0)](https://render.githubusercontent.com/render/math?math=%5Cvec%5Cmu_i%28s%2C%20%5Cvec%200%29
"\\vec\\mu_i(s, \\vec 0)"), the fully drawn curve is the individual
specific curve, ![\\vec\\mu\_i(s, \\vec
U\_i)](https://render.githubusercontent.com/render/math?math=%5Cvec%5Cmu_i%28s%2C%20%5Cvec%20U_i%29
"\\vec\\mu_i(s, \\vec U_i)"), the shaded areas are a pointwise 95%
interval for each mean curve, and the points are observed markers,
![\\vec
y\_{ij}](https://render.githubusercontent.com/render/math?math=%5Cvec%20y_%7Bij%7D
"\\vec y_{ij}").

``` r
r_n_marker <- function(id)
  # the number of markers is Poisson distributed
  rpois(1, 10) + 1L
r_obs_time <- function(id, n_markes)
  # the observations times are uniform distributed 
  sort(runif(n_markes, 0, 2))

Psi <- structure(c(0.18, 0.05, -0.05, 0.1, -0.02, 0.06, 0.05, 0.34, -0.25, 
                   -0.06, -0.03, 0.29, -0.05, -0.25, 0.24, 0.04, 0.04, 
                   -0.12, 0.1, -0.06, 0.04, 0.34, 0, -0.04, -0.02, -0.03, 
                   0.04, 0, 0.1, -0.08, 0.06, 0.29, -0.12, -0.04, -0.08, 
                   0.51), .Dim = c(6L, 6L))
B <- structure(c(-0.57, 0.17, -0.48, 0.58, 1, 0.86), .Dim = 3:2)
sig <- diag(c(.6, .3)^2)

# function to simulate a given number of individuals' markers' latent means
# and observed values
show_mark_mean <- function(B, Psi, sigma, m_func, g_func, ymax){
  tis <- seq(0, ymax, length.out = 100)
  Psi_chol <- chol(Psi)
  
  par_old <- par(no.readonly = TRUE)
  on.exit(par(par_old))
  par(mar = c(4, 3, 1, 1), mfcol = c(4, 4))
  
  sigma_chol <- chol(sigma)
  n_y <- NCOL(sigma_chol)
  for(i in 1:16){
    U <- draw_U(Psi_chol, n_y = n_y)
    y_non_rng <- t(eval_marker(tis, B = B, g_func = g_func, U = NULL, 
                               offset = NULL, m_func = m_func))
    y_rng     <- t(eval_marker(tis, B = B, g_func = g_func, U = U, 
                               offset = NULL, m_func = m_func))

    sds <- sapply(tis, function(ti){
      M <- (diag(n_y) %x% m_func(ti))
      G <- (diag(n_y) %x% g_func(ti))
      sds <- sqrt(diag(tcrossprod(M %*% Psi, M)))
      cbind(drop(G %*% c(B)) - 1.96 * sds,
            drop(G %*% c(B)) + 1.96 * sds)
    }, simplify = "array")
    lbs <- t(sds[, 1, ])
    ubs <- t(sds[, 2, ])

    y_obs <- sim_marker(B = B, U = U, sigma_chol = sigma_chol, 
                        m_func = m_func, r_n_marker = r_n_marker, 
                        r_obs_time = r_obs_time, g_func = g_func, 
                        offset = NULL)

    matplot(tis, y_non_rng, type = "l", lty = 2, ylab = "", xlab = "Time",
            ylim = range(y_non_rng, y_rng, lbs, ubs, y_obs$y_obs))
    matplot(tis, y_rng    , type = "l", lty = 1, add = TRUE)
    matplot(y_obs$obs_time, y_obs$y_obs, type = "p", add = TRUE, pch = 3:4)

    polygon(c(tis, rev(tis)), c(lbs[, 1], rev(ubs[, 1])), border = NA,
            col = rgb(0, 0, 0, .1))
    polygon(c(tis, rev(tis)), c(lbs[, 2], rev(ubs[, 2])), border = NA,
            col = rgb(1, 0, 0, .1))

  }
  invisible()
}
set.seed(1)
show_mark_mean(B = B, Psi = Psi, sigma = sig, m_func = m_func, 
               g_func = g_func, ymax = 2)
```

<img src="man/figures/README-poly_marker_setup-1.png" width="100%" />

As an example, we simulate the random effects and plot the conditional
hazards and survival functions. We start by assigning the association
parameter,
![\\vec\\alpha](https://render.githubusercontent.com/render/math?math=%5Cvec%5Calpha
"\\vec\\alpha").

``` r
alpha <- c(.5, .9)

# function to plot simulated conditional hazards and survival 
# functions
sim_surv_curves <- function(sig, Psi, delta, omega, alpha, B, m_func, 
                            g_func, b_func, ymax) {
  par_old <- par(no.readonly = TRUE)
  on.exit(par(par_old))
  par(mfcol = c(1, 2), mar = c(5, 5, 1, 1))

  # hazard functions
  tis <- seq(0, ymax, length.out = 50)
  n_y <- NCOL(sig)
  Us <- replicate(100, draw_U(chol(Psi), n_y = n_y), 
                  simplify = "array")

  hz <- apply(Us, 3L, function(U)
    vapply(tis, function(x)
      exp(drop(delta + b_func(x) %*% omega +
                 alpha %*% eval_marker(ti = x, B = B, m_func = m_func, 
                                       g_func = g_func, U = U, 
                                       offset = NULL))),
      FUN.VALUE = numeric(1L)))

  matplot(tis, hz, lty = 1, type = "l", 
          col = rgb(0, 0, 0, .1), xaxs = "i", bty = "l", yaxs = "i", 
          ylim = c(0, max(hz, na.rm = TRUE)), xlab = "time", 
          ylab = "Hazard")
  grid()

  # survival functions
  ys <- apply(Us, 3L, surv_func_joint,
              ti = tis, B = B, omega = omega, delta = delta,
              alpha = alpha, b_func = b_func, m_func = m_func, 
              gl_dat = gl_dat, g_func = g_func, offset = NULL)

  matplot(tis, ys, lty = 1, type = "l", col = rgb(0, 0, 0, .1),
          xaxs = "i", bty = "l", yaxs = "i", ylim = c(0, 1),
          xlab = "time", ylab = "Survival probability")
  grid()
}

set.seed(1)
sim_surv_curves(sig = sig, Psi = Psi, delta = delta, omega = omega, 
                alpha = alpha, B = B, m_func = m_func, g_func = g_func, 
                b_func = b_func, ymax = 2)
```

<img src="man/figures/README-poly_show_cond_surv-1.png" width="100%" />

We end by assigning the functions to get the covariates, coefficients
for the fixed effects, the left-truncation function, and right-censoring
function.

``` r
r_z <- function(id)
  # returns a design matrix for a dummy setup
  cbind(1, (id %% 3) == 1, (id %% 3) == 2)
r_z(1:6) # covariates for the first six individuals
#>      [,1] [,2] [,3]
#> [1,]    1    1    0
#> [2,]    1    0    1
#> [3,]    1    0    0
#> [4,]    1    1    0
#> [5,]    1    0    1
#> [6,]    1    0    0

# same covariates for the fixed time-invariant effects for the marker
r_x <- r_z

r_left_trunc <- function(id)
  # no left-truncation
  0
r_right_cens <- function(id)
  # right-censoring time is exponentially distributed
  rexp(1, rate = .5)

# fixed effect coefficients in the hazard
delta_vec <- c(delta, -.5, .5)
# fixed time-invariant effect coefficients in the marker process 
gamma <- matrix(c(.25, .5, 0, -.4, 0, .3), 3, 2)
```

A full data set can now be simulated as follows. We consider the time it
takes in seconds by using the `system.time` function.

``` r
set.seed(70483614)
system.time(dat <- sim_joint_data_set(
  n_obs = 1000L, B = B, Psi = Psi, omega = omega, delta = delta_vec,
  alpha = alpha, sigma = sig, b_func = b_func, g_func = g_func,
  m_func = m_func, gl_dat = gl_dat, r_z = r_z, r_left_trunc = r_left_trunc, 
  r_right_cens = r_right_cens, r_n_marker = r_n_marker, 
  r_obs_time = r_obs_time, y_max = 2, gamma = gamma, r_x = r_x))
#>    user  system elapsed 
#>   0.464   0.047   0.511
```

The first entries of the survival data and the observed markers looks as
follows.

``` r
# survival data
head(dat$survival_data)
#>   Z1 Z2 Z3 left_trunc     y event id
#> 1  1  1  0          0 0.869  TRUE  1
#> 2  1  0  1          0 0.320  TRUE  2
#> 3  1  0  0          0 1.535  TRUE  3
#> 4  1  1  0          0 0.129  TRUE  4
#> 5  1  0  1          0 2.000 FALSE  5
#> 6  1  0  0          0 0.369 FALSE  6

# marker data
head(dat$marker_data, 10)
#>    obs_time    Y1      Y2 X1 X2 X3 id
#> 1    0.4392 1.984  0.8920  1  1  0  1
#> 2    0.6009 1.651  0.2385  1  1  0  1
#> 3    0.7747 0.819  0.1867  1  1  0  1
#> 4    0.0702 3.796 -0.5575  1  0  1  2
#> 5    0.1186 3.122  0.2074  1  0  1  2
#> 6    0.2737 2.977 -0.2365  1  0  1  2
#> 7    0.2829 2.506 -0.2323  1  0  1  2
#> 8    0.1121 2.211  0.1531  1  0  0  3
#> 9    0.3337 1.595  0.0196  1  0  0  3
#> 10   0.3430 0.294 -0.1353  1  0  0  3
```

To illustrate that we simulate from the correct model, we can estimate a
linear mixed models for the markers as follows.

``` r
library(lme4)

# estimate the linear mixed model (skip this if you want and look at the 
# estimates in the end)
local({
  m_dat <- dat$marker_data
  n_y <- NCOL(gamma)
  d_x <- NROW(gamma)
  d_g <- NROW(B)
  d_m <- NROW(Psi) / n_y
  
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
          XXTHEVARIABLEXX : g_func(ti) - 1L +
          (XXTHEVARIABLEXX : m_func(ti) - 1L | i),
        list(ti = as.name("obs_time"), i = as.name("id")))
    else 
      frm <- substitute(
        XXTHEVALUEXX ~
          (XXTHEVARIABLEXX : m_func(ti, m_ks) - 1L | i),
        list(ti = as.name("obs_time"), i = as.name("id")))
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
#>         [,1]      [,2]
#> [1,]  0.1684 -0.384104
#> [2,]  0.5224  0.000604
#> [3,] -0.0767  0.219848
#> 
#> $B
#>        [,1]  [,2]
#> [1,] -0.569 0.484
#> [2,]  0.282 0.988
#> [3,] -0.466 0.893
#> 
#> $Psi
#>         [,1]    [,2]     [,3]     [,4]    [,5]    [,6]
#> [1,]  0.5247  0.0687 -0.16088  0.12023 -0.0285  0.1125
#> [2,]  0.0687  0.3572 -0.24080 -0.03208 -0.0450  0.2571
#> [3,] -0.1609 -0.2408  0.25861 -0.00238  0.0464 -0.1342
#> [4,]  0.1202 -0.0321 -0.00238  0.16060 -0.0057  0.0470
#> [5,] -0.0285 -0.0450  0.04644 -0.00570  0.0533 -0.0832
#> [6,]  0.1125  0.2571 -0.13418  0.04703 -0.0832  0.4165
#> 
#> $Sigma
#>       [,1]  [,2]
#> [1,] 0.225 0.000
#> [2,] 0.000 0.225
```

Although we assume equal noise variance, the estimates are close to the
true values.

``` r
gamma
#>      [,1] [,2]
#> [1,] 0.25 -0.4
#> [2,] 0.50  0.0
#> [3,] 0.00  0.3
B
#>       [,1] [,2]
#> [1,] -0.57 0.58
#> [2,]  0.17 1.00
#> [3,] -0.48 0.86
Psi
#>       [,1]  [,2]  [,3]  [,4]  [,5]  [,6]
#> [1,]  0.18  0.05 -0.05  0.10 -0.02  0.06
#> [2,]  0.05  0.34 -0.25 -0.06 -0.03  0.29
#> [3,] -0.05 -0.25  0.24  0.04  0.04 -0.12
#> [4,]  0.10 -0.06  0.04  0.34  0.00 -0.04
#> [5,] -0.02 -0.03  0.04  0.00  0.10 -0.08
#> [6,]  0.06  0.29 -0.12 -0.04 -0.08  0.51
sig
#>      [,1] [,2]
#> [1,] 0.36 0.00
#> [2,] 0.00 0.09
```

We then fit a Cox model with only the observed markers (likely biased
but gives us an idea about whether we are using the correct model).

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
  for(i in seq_along(delta_vec)){
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
#>   n= 4042, number of events= 535 
#> 
#>       coef exp(coef) se(coef)     z Pr(>|z|)    
#> Z1      NA        NA   0.0000    NA       NA    
#> Z2 -0.8688    0.4194   0.1168 -7.44  1.0e-13 ***
#> Z3  0.4566    1.5787   0.1007  4.53  5.8e-06 ***
#> Y1  0.5474    1.7287   0.0419 13.06  < 2e-16 ***
#> Y2  0.4519    1.5712   0.0489  9.23  < 2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#>    exp(coef) exp(-coef) lower .95 upper .95
#> Z1        NA         NA        NA        NA
#> Z2     0.419      2.384     0.334     0.527
#> Z3     1.579      0.633     1.296     1.923
#> Y1     1.729      0.578     1.592     1.877
#> Y2     1.571      0.636     1.428     1.729
#> 
#> Concordance= 0.693  (se = 0.012 )
#> Likelihood ratio test= 270  on 4 df,   p=<2e-16
#> Wald test            = 271  on 4 df,   p=<2e-16
#> Score (logrank) test = 278  on 4 df,   p=<2e-16
```

This is close-ish to the true values.

``` r
delta_vec
#> [1] -0.5 -0.5  0.5
alpha
#> [1] 0.5 0.9
```

### Using Derivatives

It is possible to use derivatives of the latent mean, ![\\vec\\mu\_i(s,
\\vec
u)](https://render.githubusercontent.com/render/math?math=%5Cvec%5Cmu_i%28s%2C%20%5Cvec%20u%29
"\\vec\\mu_i(s, \\vec u)"), with respect to time in the hazard. As an
example, we consider the first-order derivative below and plot
conditional hazards and survival functions simulated from the new model.

``` r
g_func_surv <- function(x){
  x <- x - 1
  cbind(3 * x^2, 2 * x, 1)
}
m_func_surv <- function(x){
  x <- x - 1
  cbind(2 * x, 1, 0)
}

set.seed(1)
sim_surv_curves(sig = sig, Psi = Psi, delta = delta, omega = omega, 
                alpha = alpha, B = B, m_func = m_func_surv, 
                g_func = g_func_surv, b_func = b_func, ymax = 2)
```

<img src="man/figures/README-poly_derivs-1.png" width="100%" />

A new data set can now be simulated as follows.

``` r
set.seed(70483614)
system.time(dat <- sim_joint_data_set(
  n_obs = 1000L, B = B, Psi = Psi, omega = omega, delta = delta_vec,
  alpha = alpha, sigma = sig, b_func = b_func, g_func = g_func,
  m_func = m_func, gl_dat = gl_dat, r_z = r_z, r_left_trunc = r_left_trunc, 
  r_right_cens = r_right_cens, r_n_marker = r_n_marker, 
  r_obs_time = r_obs_time, y_max = 2, gamma = gamma, r_x = r_x, 
  # the additions
  m_func_surv = m_func_surv, g_func_surv = g_func_surv, 
  use_fixed_latent = FALSE))
#>    user  system elapsed 
#>   0.567   0.035   0.601
```

The first entries of the new data looks as follows.

``` r
# survival data
head(dat$survival_data)
#>   Z1 Z2 Z3 left_trunc     y event id
#> 1  1  1  0          0 1.038 FALSE  1
#> 2  1  0  1          0 1.704  TRUE  2
#> 3  1  0  0          0 1.752  TRUE  3
#> 4  1  1  0          0 1.170  TRUE  4
#> 5  1  0  1          0 1.343  TRUE  5
#> 6  1  0  0          0 0.369 FALSE  6

# marker data
head(dat$marker_data, 10)
#>    obs_time    Y1     Y2 X1 X2 X3 id
#> 1    0.4392 1.984  0.892  1  1  0  1
#> 2    0.6009 1.651  0.239  1  1  0  1
#> 3    0.7747 0.819  0.187  1  1  0  1
#> 4    0.8931 0.853  0.521  1  1  0  1
#> 5    0.9442 1.522 -0.216  1  1  0  1
#> 6    0.0702 3.796 -0.558  1  0  1  2
#> 7    0.1186 3.122  0.207  1  0  1  2
#> 8    0.2737 2.977 -0.237  1  0  1  2
#> 9    0.2829 2.506 -0.232  1  0  1  2
#> 10   0.4253 2.098 -0.124  1  0  1  2
```

## Example with Natural Cubic Splines

In this section, we will use natural cubic splines for the time-varying
basis functions. We start by assigning all the variables that we will
pass to the functions in the package.

``` r
# quadrature nodes
gl_dat <- get_gl_rule(30L)

# knots for the spline functions
b_ks <- seq(log(1), log(10), length.out = 4)
m_ks <- seq(0, 10, length.out = 3)
g_ks <- m_ks

# simulation functions
r_n_marker <- function(id)
  rpois(1, 10) + 1L
r_obs_time <- function(id, n_markes)
  sort(runif(n_markes, 0, 10))
r_z <- function(id)
  as.numeric(runif(1L) > .5)
r_x <- function(id)
  numeric()
r_left_trunc <- function(id)
   rbeta(1, 1, 2) * 3
r_right_cens <- function(id)
  rbeta(1, 2, 1) * 6 + 4

# model parameters
omega <- c(-0.96, -2.26, -3.04, .45)
Psi <- structure(c(1.08, 0.12, -0.36, -0.48, 0.36, -0.12, 0.12, 0.36,
                   0, 0.12, 0, -0.12, -0.36, 0, 0.84, 0.12, 0.12, 0.12, -0.48, 0.12,
                   0.12, 0.84, -0.12, 0.24, 0.36, 0, 0.12, -0.12, 0.84, -0.12, -0.12,
                   -0.12, 0.12, 0.24, -0.12, 0.6), .Dim = c(6L, 6L))
B <- structure(c(0.97, 0.01, -0.07, -0.78, -0.86, 0.98), .Dim = 3:2)
sig <- diag(c(.2, .1)^2)
alpha <- c(0.7, 0.6)
delta <- 0.
gamma <- numeric()
```

``` r
# spline functions
b_func <- get_ns_spline(b_ks, do_log = TRUE)
m_func <- get_ns_spline(m_ks, do_log = FALSE)
g_func <- get_ns_spline(g_ks, do_log = FALSE)
```

Notice that we use the `get_ns_spline` functions which reduces the
computation time a lot. We show the baseline hazard function and the
survival function without the marker below.

``` r
# hazard function without marker
par(mar = c(5, 5, 1, 1), mfcol = c(1, 2))
plot(function(x) exp(drop(b_func(x) %*% omega)),
     xlim = c(1e-8, 10), ylim = c(0, .61), xlab = "Time",
     ylab = "Hazard (no marker)", xaxs = "i", bty = "l")
grid()

# survival function without marker
plot(function(x) eval_surv_base_fun(x, omega = omega, b_func = b_func, 
                                    gl_dat = gl_dat, delta = delta), 
     xlim = c(1e-4, 10),
     xlab = "Time", ylab = "Survival probability (no marker)", xaxs = "i",
     yaxs = "i", bty = "l", ylim = c(0, 1.01))
grid()
```

<img src="man/figures/README-plot_wo_marker-1.png" width="100%" />

Next, we simulate individual specific markers. Each plot is for a given
individual.

``` r
set.seed(1)
show_mark_mean(B = B, Psi = Psi, sigma = sig, m_func = m_func, 
               g_func = g_func, ymax = 10)
```

<img src="man/figures/README-show_sim_marker-1.png" width="100%" />

We sample a number of random effects and plot the hazard curves and
survival functions given these random effects below.

``` r
set.seed(1)
sim_surv_curves(sig = sig, Psi = Psi, delta = delta, omega = omega, 
                alpha = alpha, B = B, m_func = m_func, g_func = g_func, 
                b_func = b_func, ymax = 10)
```

<img src="man/figures/README-show_draw_surv_curves-1.png" width="100%" />

We end by drawing a data set.

``` r
set.seed(70483614)
delta_vec <- 1
system.time(dat <- sim_joint_data_set(
  n_obs = 1000L, B = B, Psi = Psi, omega = omega, delta = delta_vec,
  alpha = alpha, sigma = sig, b_func = b_func, g_func = g_func,
  m_func = m_func, gl_dat = gl_dat, r_z = r_z, r_left_trunc = r_left_trunc, 
  r_right_cens = r_right_cens, r_n_marker = r_n_marker, 
  r_obs_time = r_obs_time, y_max = 10, gamma = gamma, r_x = r_x))
#>    user  system elapsed 
#>   0.593   0.039   0.632
```

Finally, we show a few of the first rows along with some summary
statistics.

``` r
# survival data
head(dat$survival_data)
#>   Z1 left_trunc    y event id
#> 1  0     0.3830 7.85 FALSE  1
#> 2  1     1.0779 2.10  TRUE  2
#> 3  0     0.8923 1.37  TRUE  3
#> 4  0     0.9055 8.04 FALSE  4
#> 5  1     0.8378 6.22  TRUE  5
#> 6  1     0.0649 1.84  TRUE  6

# marker data
head(dat$marker_data, 10)
#>    obs_time    Y1     Y2 id
#> 1      1.32 0.824 -2.165  1
#> 2      1.85 0.163 -2.132  1
#> 3      2.20 1.052 -1.914  1
#> 4      3.00 0.831 -1.941  1
#> 5      3.59 1.303 -1.840  1
#> 6      4.91 1.061 -1.660  1
#> 7      4.92 0.803 -1.770  1
#> 8      5.74 0.595 -1.520  1
#> 9      5.87 0.432 -1.364  1
#> 10     1.77 1.285 -0.497  2

# rate of observed events
mean(dat$survival_data$event) 
#> [1] 0.742

# mean event time
mean(subset(dat$survival_data, event           )$y)
#> [1] 3.49

# mean event time for the two group
mean(subset(dat$survival_data, event & Z1 == 1L)$y)
#> [1] 3.08
mean(subset(dat$survival_data, event & Z1 == 0L)$y)
#> [1] 4.16

# quantiles of the event time
quantile(subset(dat$survival_data, event)$y)
#>     0%    25%    50%    75%   100% 
#> 0.0707 1.6318 2.7531 5.2749 9.5271

# fraction of observed markers per individual
NROW(dat$marker_data) / NROW(dat$survival_data)
#> [1] 4.08
```

## References

<div id="refs" class="references">

<div id="ref-Crowther13">

Crowther, Michael J., and Paul C. Lambert. 2013. “Simulating
Biologically Plausible Complex Survival Data.” *Statistics in Medicine*
32 (23): 4118–34. <https://doi.org/10.1002/sim.5823>.

</div>

</div>
