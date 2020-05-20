.time_eps <- sqrt(.Machine$double.eps)

#' Faster Pointwise function than ns
#'
#' @param knots numeric vector with boundary and interior knots.
#' @param intercept logical for whether to include an intercept.
#' @param do_log logical for whether to evaluate the spline at \code{log(x)}
#'               or \code{x}.
#'
#' @examples
#' # compare with splines
#' library(splines)
#' library(SimSurvNMarker)
#' xs <- seq(1, 5, length.out = 10L)
#' bks <- c(1, 5)
#' iks <- 2:4
#'
#' # we get the same
#' stopifnot(isTRUE(all.equal(
#'   unclass(ns(xs, knots = iks, Boundary.knots = bks, intercept = TRUE)),
#'   get_ns_spline(knots = sort(c(iks, bks)), intercept = TRUE,
#'                 do_log = FALSE)(xs),
#'   check.attributes = FALSE)))
#'
#' stopifnot(isTRUE(all.equal(
#'   unclass(ns(log(xs), knots = log(iks), Boundary.knots = log(bks),
#'              intercept = TRUE)),
#'   get_ns_spline(knots = log(sort(c(iks, bks))), intercept = TRUE,
#'                 do_log = TRUE)(xs),
#'   check.attributes = FALSE)))
#'
#' # the latter is faster
#' system.time(
#'   replicate(100,
#'             ns(xs, knots = iks, Boundary.knots = bks, intercept = TRUE)))
#' system.time(
#'   replicate(100,
#'             get_ns_spline(knots = sort(c(iks, bks)), intercept = TRUE,
#'                           do_log = FALSE)(xs)))
#' func <- get_ns_spline(knots = sort(c(iks, bks)), intercept = TRUE,
#'                       do_log = FALSE)
#' system.time(replicate(100, func(xs)))
#'
#' @export
get_ns_spline <- function(knots, intercept = TRUE, do_log = TRUE)
  with(new.env(), {
    ptr_obj <- get_ns_ptr(knots          = knots[-c(1L, length(knots))],
                          boundary_knots = knots[ c(1L, length(knots))],
                          intercept = intercept)
    out <- if(do_log)
      function(x)
        ns_cpp(x = log(x), ns_ptr = ptr_obj)
    else
      function(x)
        ns_cpp(x =     x , ns_ptr = ptr_obj)

    attributes(out) <- list(
      knots = knots, intercept = intercept, do_log = do_log)
    out
  })

.surv_func_inner <- function(ti, omega, b_func, gl_dat){
  lb <- .time_eps
  ub <- ti
  if(ti < lb)
    stop("too small ti")

  nodes <- (ub - lb) / 2 * gl_dat$node + (ub + lb) / 2
  f <- exp(drop(b_func(nodes) %*% omega))
  -(ub - lb) / 2 * drop(gl_dat$weight %*% f)
}

#' Evaluates the Survival Function without a Marker
#'
#' @param ti time points.
#' @param omega coefficient for baseline hazard.
#' @param b_func basis function for baseline hazard like \code{\link{poly}}.
#' @param gl_dat Gauss–Legendre quadrature data.
#' @param delta offset on the log hazard scale. Use \code{NULL} if there is no effect.
#'
#' @examples
#' # Example of a hazard function
#' b_func <- function(x)
#'   cbind(1, sin(2 * pi * x), x)
#' omega <- c(-3, 3, .25)
#' haz_fun <- function(x)
#'   exp(drop(b_func(x) %*% omega))
#'
#' plot(haz_fun, xlim = c(0, 10))
#'
#' # plot the hazard
#' library(SimSurvNMarker)
#' gl_dat <- get_gl_rule(60L)
#' plot(function(x) eval_surv_base_fun(ti = x, omega = omega,
#'                                     b_func = b_func, gl_dat = gl_dat),
#'      xlim = c(1e-4, 10), ylim = c(0, 1), bty = "l", xlab = "time",
#'      ylab = "Survival", yaxs = "i")
#' @export
eval_surv_base_fun <- function(
  ti, omega, b_func, gl_dat = get_gl_rule(30L), delta = 0){
  cum_haz_fac <- vapply(
    ti, .surv_func_inner, FUN.VALUE = numeric(1L),
    omega = omega, b_func = b_func, gl_dat = gl_dat)

  if(length(delta) > 0)
    exp(exp(delta) * cum_haz_fac)
  else
    exp(cum_haz_fac)
}

.eval_marker_inner <- function(ti, B, m_func)
  eval_marker_cpp(B = B, m = m_func(ti))

#' Fast Evaluation of Time-varying Marker Mean Term
#'
#' @param ti time points.
#' @param B coefficient matrix for time-varying fixed effects. Use \code{NULL} if there is no effect.
#' @param g_func basis function for \code{B} like \code{\link{poly}}.
#' @param U random effect matrix for time-varying random effects. Use \code{NULL} if there is no effect.
#' @param m_func basis function for \code{U} like \code{\link{poly}}.
#' @param offset vector with non-time-varying fixed effects.
#'
#' @examples
#' # compare R version with this function
#' library(SimSurvNMarker)
#' set.seed(1)
#' n <- 100L
#' n_y <- 3L
#'
#' ti <- seq(0, 1, length.out = n)
#' offset <- runif(n_y)
#' B <- matrix(runif(5L * n_y), nr = 5L)
#' g_func <- function(x)
#'   cbind(1, x, x^2, x^3, x^4)
#' U <- matrix(runif(3L * n_y), nr = 3L)
#' m_func <- function(x)
#'   cbind(1, x, x^2)
#'
#' r_version <- function(ti, B, g_func, U, m_func, offset){
#'   n_y <- NCOL(B)
#'   B <- c(B)
#'   U <- c(U)
#'   func <- function(ti)
#'     drop((diag(n_y) %x% g_func(ti)) %*% B) +
#'       drop((diag(n_y) %x% m_func(ti)) %*% U)
#'
#'   vapply(ti, func, numeric(n_y)) + offset
#' }
#'
#' # check that we get the same
#' stopifnot(isTRUE(all.equal(
#'   c(r_version  (ti[1], B, g_func, U, m_func, offset)),
#'     eval_marker(ti[1], B, g_func, U, m_func, offset))))
#' stopifnot(isTRUE(all.equal(
#'   r_version  (ti, B, g_func, U, m_func, offset),
#'   eval_marker(ti, B, g_func, U, m_func, offset))))
#'
#' # check the computation time
#' system.time(replicate(100, r_version  (ti, B, g_func, U, m_func, offset)))
#' system.time(replicate(100, eval_marker(ti, B, g_func, U, m_func, offset)))
#'
#' @export
eval_marker <- function(ti, B, g_func, U, m_func, offset){
  out <- 0.
  if(length(B) > 0)
    out <- out + eval_marker_cpp(B = B, m = g_func(ti))
  if(length(U) > 0)
    out <- out + eval_marker_cpp(B = U, m = m_func(ti))
  if(length(offset) > 0)
    out <- out + offset

  drop(out)
}

#' Samples from a Multivariate Normal Distribution
#'
#' @param Psi_chol Cholesky decomposition of the covariance matrix.
#' @param n_y number of markers.
#'
#' @examples
#' library(SimSurvNMarker)
#' set.seed(1)
#' n_y <- 2L
#' K <- 3L * n_y
#' Psi <- drop(rWishart(1, K, diag(K)))
#' Psi_chol <- chol(Psi)
#'
#' # example
#' dim(draw_U(Psi_chol, n_y))
#' samples <- replicate(100, draw_U(Psi_chol, n_y))
#' samples <- t(apply(samples, 3, c))
#'
#' colMeans(samples) # ~ zeroes
#' cov(samples) # ~ Psi
#'
#' @importFrom stats rnorm
#' @export
draw_U <- function(Psi_chol, n_y){
  K <- NCOL(Psi_chol)
  structure(rnorm(K) %*%  Psi_chol, dim = c(K %/% n_y, n_y))
}

#' Simulate a Number of Observed Marker for an Individual
#'
#' @param sigma_chol Cholesky decomposition of the noise's covariance matrix.
#' @param r_n_marker function to generate the number of observed markers.
#' @param r_obs_time function to generate the observations times given the
#'                   number of observed markers. It should take a single
#'                   integer.
#' @inheritParams eval_marker
#'
#' @importFrom stats rnorm
#' @export
sim_marker <- function(B, U, sigma_chol, r_n_marker, r_obs_time, m_func,
                       g_func, offset){
  n_markes <- r_n_marker()
  obs_time <- r_obs_time(n_markes)
  n_y <- NCOL(sigma_chol)

  noise <- matrix(rnorm(n_markes * n_y), ncol = n_y)  %*% sigma_chol
  mu <- eval_marker(
    ti = obs_time, B = B, g_func = g_func, m_func = m_func, offset = offset,
    U = U)
  y_obs <- if(is.vector(mu))
      mu  + noise
  else
    t(mu) + noise

  list(obs_time = obs_time, y_obs = y_obs)
}

.surv_func_joint_inner <- function(ti, B, U, omega, alpha, b_func, m_func,
                                   g_func, gl_dat){
  lb <- .time_eps
  ub <- ti
  if(ti < lb)
    stop("too small ti")

  nodes <- (ub - lb) / 2 * gl_dat$node + (ub + lb) / 2
  f <- exp(
    drop(b_func(nodes) %*% omega +
           drop(alpha %*% eval_marker(
             ti = nodes, B = B, U = U, m_func = m_func, g_func = g_func,
             offset = NULL))))
  -(ub - lb) / 2 * drop(gl_dat$weight %*% f)
}

#' Evaluate the Conditional Survival Function Given the Random Effect
#'
#' @inheritParams eval_surv_base_fun
#' @inheritParams eval_marker
#' @param alpha numeric vector with association parameters.
#'
#' @export
surv_func_joint <- function(ti, B, U, omega, delta, alpha, b_func, m_func,
                            gl_dat = get_gl_rule(30L), g_func, offset){
  cum_haz_fac <- vapply(
    ti, .surv_func_joint_inner, FUN.VALUE = numeric(1L),
    B = B, U = U, omega = omega, alpha = alpha, b_func = b_func,
    m_func = m_func, g_func = g_func, gl_dat = gl_dat)

  if(length(delta) > 0)
    cum_haz_fac <- cum_haz_fac * exp(delta)
  if(length(offset) > 0)
    cum_haz_fac <- cum_haz_fac * exp(drop(offset %*% alpha))
  exp(cum_haz_fac)
}

list_of_lists_to_data_frame <- function(dat)
  as.data.frame(do.call(mapply, c(list(FUN = function(...){
    if(is.matrix(...elt(1L)))
      do.call(rbind, list(...))
    else
      do.call(c, list(...))
  }, SIMPLIFY = FALSE), dat)))

#' Simulate Individuals from a Joint Survival and Marker Model
#'
#' @param n_obs integer with the number of individuals to draw.
#' @param Psi the random effect covariance matrix.
#' @param sigma the noise's covariance matrix.
#' @param gamma coefficient matrix for the non-time-varying fixed effects. Use \code{NULL} if there is no effect.
#' @param r_z generator for the covariates in the log hazard.
#' @param r_left_trunc generator for the left-truncation time.
#' @param r_right_cens generator for the right-censoring time.
#' @param r_x generator for the covariates in for the markers.
#' @param y_max maximum survival time before administrative censoring.
#' @inheritParams surv_func_joint
#' @inheritParams sim_marker
#'
#' @importFrom stats uniroot runif
#' @export
sim_joint_data_set <- function(
  n_obs, B, Psi, omega, delta, alpha, sigma, gamma, b_func, m_func, g_func,
  gl_dat = get_gl_rule(30L), r_z, r_left_trunc, r_right_cens,
  r_n_marker, r_x, r_obs_time, y_max){
  Psi_chol <- chol(Psi)
  y_min <- .time_eps
  sigma_chol <- chol(sigma)

  # checks
  n_y <- NCOL(sigma)
  d_m <- NROW(B)
  d_b <- length(omega)
  d_z <- length(delta)
  d_x <- length(gamma) / n_y
  K <- length(m_func(1L)) * n_y

  local({
    old_seed <- .GlobalEnv$.Random.seed
    on.exit(.GlobalEnv$.Random.seed <- old_seed)

    stopifnot(
      length(B) == 0L || is.matrix(B),
      length(B) == 0L || is.numeric(B),
      length(B) == 0L || NCOL(B) == n_y,
      is.matrix(Psi), is.numeric(Psi), NCOL(Psi) == K,
      is.matrix(sigma), is.numeric(sigma),
      is.numeric(alpha), length(alpha) == n_y,
      is.numeric(b_func(1)), length(b_func(1)) == d_b,
      is.numeric(m_func(1)),
      is.numeric(g_func(1)),
      is.numeric(gl_dat$node), is.numeric(gl_dat$weight),
      length(gl_dat$node) == length(gl_dat$weight),
      is.numeric(r_z()), length(r_z()) == d_z,
      is.numeric(r_x()), length(r_x()) == d_x,
      length(gamma) == 0L || is.matrix(gamma),
      length(gamma) == 0L || is.numeric(gamma),
      length(gamma) == 0L || NCOL(gamma) == n_y,
      is.numeric(r_left_trunc()),
      is.numeric(r_right_cens()),
      is.integer(r_n_marker()),
      is.numeric(r_obs_time(1L)))
  })

  out <- rep(list(list()), n_obs)
  for(i in seq_along(out))
    out[[i]] <- {
      # keep going till we get a none left-truncated event
      has_sample <- FALSE

      z_delta <- if(d_z == 0L){
        z <- numeric()
        NULL
      }
      else {
        z <- r_z()
        drop(z %*% delta)
      }

      mu_offest <- if(d_x == 0L){
        x <- numeric()
        NULL
      } else {
        x <- r_x()
        drop(eval_marker_cpp(B = gamma, m = x))
      }

      while(!has_sample){
        U <- draw_U(Psi_chol, n_y = n_y)
        unif <- runif(1)
        fun <- function(x){
          surv <- surv_func_joint(
            ti = x, B = B, U = U, omega = omega,
            delta = z_delta, alpha = alpha, b_func = b_func, m_func = m_func,
            gl_dat = gl_dat, g_func = g_func, offset = mu_offest)

          surv - unif
        }


        f_lower <- fun(y_min)
        f_upper <- fun(y_max)
        y <- if(f_lower < 0) y_min else
          if(f_upper > 0)
            y_max else {
              root <- uniroot(fun, interval = c(y_min, y_max),
                              f.lower = f_lower, f.upper = f_upper)
              root$root
            }

        left_trunc <- r_left_trunc()
        if(y < left_trunc)
          next

        # simulate right truncation time
        right_cens <- r_right_cens()

        # simulate the marker
        keep <- logical()
        while(sum(keep) < 1L){
          markers <- sim_marker(
            B = B, U = U, sigma_chol = sigma_chol, r_n_marker = r_n_marker,
            r_obs_time = r_obs_time, m_func = m_func, g_func = g_func,
            offset = mu_offest)
          keep <- markers$obs_time < min(y, right_cens) &
            markers$obs_time > left_trunc

          if(sum(keep) < 1L)
            next

          has_sample <- TRUE
          markers <- cbind(obs_time = markers$obs_time, markers$y_obs)
          markers <- markers[keep, , drop = FALSE]
          colnames(markers)[-1] <- paste0("Y", 1:(NCOL(markers) - 1L))
        }
      }

      list(z = z, left_trunc = left_trunc, y = min(right_cens, y),
           event = y < right_cens, U = U, markers = markers, x = x)
    }

  # form data frames for estimation
  ids <- seq_along(out)
  survival_dat <- lapply(ids, function(i){
      dat_i <- out[[i]][c("z", "left_trunc", "y", "event")]
      if(d_z > 0){
        dat_i$z <- matrix(
          dat_i$z, nrow = 1L,
          dimnames = list(NULL, paste0("Z", seq_along(dat_i$z))))
        names(dat_i)[1L] <- ""
      } else
        dat_i$z <- NULL
      dat_i$id <- i
      dat_i
    })
  survival_dat <- list_of_lists_to_data_frame(survival_dat)

  marker_dat <- lapply(ids, function(i){
      markers <- out[[i]][["markers"]]
      n_y <- NROW(markers)
      if(d_x > 0){
        x <- structure(rep(out[[i]]$x, each = n_y),
                       dimnames = list(NULL, paste0("X", seq_len(d_x))),
                       dim = c(n_y, d_x))
        markers <- cbind(markers, x)
      }

      out <- list(markers)
      out$id <- rep(i, NROW(markers))
      out
    })
  marker_dat <- list_of_lists_to_data_frame(marker_dat)

  list(survival_data = survival_dat,
       marker_data   = marker_dat,
       complete_data = out,
       params        = list(
         B = B, Psi = Psi, omega = omega, delta = delta, alpha = alpha,
         sigma = sigma, b_attr = attributes(b_func),
         m_attr = attributes(m_func)))
}
