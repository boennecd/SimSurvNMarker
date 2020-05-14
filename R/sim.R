.time_eps <- sqrt(.Machine$double.eps)

#'  Faster Pointwise function than ns
#' @export
get_ns_spline <- function(knots, intercept = TRUE, do_log = TRUE)
  with(new.env(), {
    ptr_obj <- get_ns_ptr(knots          = knots[-c(1L, length(knots))],
                          boundary_knots = knots[ c(1L, length(knots))],
                          intercept = intercept)
    if(do_log)
      function(x)
        ns_cpp(x = log(x), ns_ptr = ptr_obj)
    else
      function(x)
        ns_cpp(x =     x , ns_ptr = ptr_obj)
  })

.surv_func_inner <- function(x, omega, b_func)
  exp(drop(b_func(x) %*% omega))

.eval_surv_base_fun_outer <- function(ti, omega, b_func, gl_dat){
  lb <- .time_eps
  ub <- ti
  if(ti < lb)
    stop("too small ti")

  nodes <- (ub - lb) / 2 * gl_dat$node + (ub + lb) / 2
  f <- vapply(nodes, .surv_func_inner, FUN.VALUE = numeric(1L),
              omega = omega, b_func)
  -(ub - lb) / 2 * drop(gl_dat$weight %*% f)
}

#' Evalutes the Survival Function without a Marker
#' @export
eval_surv_base_fun <- function(
  ti, omega, b_func, gl_dat = SimSurvNMarker:::get_gl_rule(30L), delta = 0){
  cum_haz_fac <- vapply(
    ti, .eval_surv_base_fun_outer, FUN.VALUE = numeric(1L),
    omega = omega, b_func = b_func, gl_dat = gl_dat)

  exp(exp(delta) * cum_haz_fac)
}

.eval_marker_inner <- function(ti, B, m_func)
  eval_marker_cpp(B = B, m = m_func(ti))

#' Fast Evaluation of Time-varying Marker Mean Term
#' @export
eval_marker <- function(ti, B, m_func)
  vapply(ti, .eval_marker_inner, FUN.VALUE = numeric(NCOL(B)),
         B = B, m_func = m_func)

#' @export
draw_U <- function(Psi_chol)
  drop(rnorm(length(B)) %*%  Psi_chol)

#' @export
sim_marker <- function(B, U, sigma_chol, r_n_marker, r_obs_time, m_func){
  n_markes <- r_n_marker()
  obs_time <- r_obs_time(n_markes)

  B_use <- B + U
  y_obs <- t(eval_marker(obs_time, B_use, m_func) +
               replicate(n_markes, drop(rnorm(NCOL(B)) %*% sigma_chol)))
  list(obs_time = obs_time, y_obs = y_obs)
}

.surv_func_joint_inner <- function(x, B, U, omega, alpha, b_func, m_func)
  exp(drop(b_func(x) %*% omega +
             alpha %*% eval_marker(ti = x, B = B + U, m_func = m_func)))

.surv_func_joint_outer <- function(ti, B, U, omega, alpha, b_func, m_func,
                                   gl_dat){
  lb <- .time_eps
  ub <- ti
  if(ti < lb)
    stop("too small ti")

  nodes <- (ub - lb) / 2 * gl_dat$node + (ub + lb) / 2
  f <- vapply(nodes, .surv_func_joint_inner, FUN.VALUE = numeric(1L),
              B = B, U = U, omega = omega, alpha = alpha, m_func = m_func,
              b_func = b_func)
  -(ub - lb) / 2 * drop(gl_dat$weight %*% f)
}

#' @export
surv_func_joint <- function(ti, B, U, omega, delta, alpha, b_func, m_func,
                            gl_dat = get_gl_rule(30L)){
  cum_haz_fac <- vapply(
    ti, .surv_func_joint_outer, FUN.VALUE = numeric(1L),
    B = B, U = U, omega = omega, alpha = alpha, b_func = b_func,
    m_func = m_func, gl_dat = gl_dat)

  exp(exp(delta) * cum_haz_fac)
}

#' @importFrom stats uniroot
#' @export
sim_joint_data_set <- function(
  n_obs, B, Psi, omega, delta, alpha, sigma, b_func, m_func,
  gl_dat = get_gl_rule(30L), r_z, r_left_trunc, r_right_cens,
  r_n_marker, r_obs_time, y_max){
  Psi_chol <- chol(Psi)
  y_min <- .time_eps
  sigma_chol <- chol(sigma)

  # checks
  local({
    old_seed <- .GlobalEnv$.Random.seed
    on.exit(.GlobalEnv$.Random.seed <- old_seed)
    n_y <- NCOL(sigma)
    d_m <- NROW(B)
    d_b <- length(omega)
    d_z <- length(delta)
    K <- d_m * n_y

    stopifnot(
      is.matrix(B), NCOL(B) == n_y,
      NCOL(Psi) == K,
      is.numeric(alpha), length(alpha) == n_y,
      is.numeric(b_func(1)), length(b_func(1)) == d_b,
      is.numeric(m_func(1)), length(m_func(1)) == d_m,
      is.numeric(gl_dat$node), is.numeric(gl_dat$weight),
      length(gl_dat$node) == length(gl_dat$weight),
      is.numeric(r_z()), length(r_z()) == d_z,
      is.numeric(r_left_trunc()),
      is.numeric(r_right_cens()),
      is.integer(r_n_marker()),
      is.numeric(r_obs_time(1L)))
  })

  out <- replicate(n_obs, {
    # keep going till we get a none left-truncated event
    has_sample <- FALSE
    z <- r_z()
    z_delta <- drop(z %*% delta)

    while(!has_sample){
      U <- draw_U(Psi_chol)
      unif <- runif(1)
      fun <- function(x){
        surv <- surv_func_joint(
          ti = x, B = B, U = U, omega = omega,
          delta = z_delta, alpha = alpha, b_func = b_func, m_func = m_func,
          gl_dat = gl_dat)

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
          r_obs_time = r_obs_time, m_func = m_func)
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
         event = y < right_cens, U = U, markers = markers)
  }, simplify = FALSE)

  # form data frames for estimation
  ids <- seq_along(out)
  survival_dat <- do.call(
    rbind, lapply(ids, function(i){
      dat_i <- out[[i]][c("z", "left_trunc", "y", "event")]
      dat_i$z <- matrix(
        dat_i$z, nr = 1L,
        dimnames = list(NULL, paste0("Z", seq_along(dat_i$z))))
      dat_i$id <- i
      names(dat_i)[1L] <- ""
      do.call(data.frame, dat_i)
    }))

  marker_dat <- do.call(
    rbind, lapply(ids, function(i){
      out <- as.data.frame(out[[i]][["markers"]])
      out$id <- i
      out
    }))

  list(survival_data = survival_dat,
       marker_data   = marker_dat,
       complete_data = out,
       params        = list(
         B = B, Psi = Psi, omega = omega, delta = delta, alpha = alpha,
         sigma = sigma, b_knots = b_ks, m_knots = m_ks))
}
