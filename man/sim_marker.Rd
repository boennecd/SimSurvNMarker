% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim.R
\name{sim_marker}
\alias{sim_marker}
\title{Simulate a Number of Observed Marker for an Individual}
\usage{
sim_marker(
  B,
  U,
  sigma_chol,
  r_n_marker,
  r_obs_time,
  m_func,
  g_func,
  offset,
  id = 1L
)
}
\arguments{
\item{B}{coefficient matrix for time-varying fixed effects.
Use \code{NULL} if there is no effect.}

\item{U}{random effects matrix for time-varying random effects.
Use \code{NULL} if there is no effects.}

\item{sigma_chol}{Cholesky decomposition of the noise's covariance matrix.}

\item{r_n_marker}{function to generate the number of observed markers.
Takes an integer for the individual's id.}

\item{r_obs_time}{function to generate the observations times given the
number of observed markers. Takes an
integer for the number of markers and an integer
for the individual's id.}

\item{m_func}{basis function for \code{U} like \code{\link{poly}}.}

\item{g_func}{basis function for \code{B} like \code{\link{poly}}.}

\item{offset}{numeric vector with non-time-varying fixed effects.}

\item{id}{integer with id passed to \code{r_n_marker} and
\code{r_obs_time}.}
}
\description{
Simulates from

\deqn{\vec U_i \sim N^{(K)}(\vec 0, \Psi)}
\deqn{\vec Y_{ij} \mid \vec U_i = \vec u_i \sim N^{(r)}(\vec \mu(s_{ij}, \vec u_i), \Sigma)}

with

\deqn{\vec\mu(s, \vec u) = \vec o + \left(I \otimes \vec g(s)^\top\right)vec(B) + \left(I \otimes \vec m(s)^\top\right) \vec u.}

The number of observations and the observations times, \eqn{s_{ij}}s, are
determined from the passed generating functions.
}
\examples{
#####
# example with polynomial basis functions
g_func <- function(x){
  x <- x - 1
  cbind(x^3, x^2, x)
}
m_func <- function(x){
  x <- x - 1
  cbind(x^2, x, 1)
}

# parameters
gamma <- matrix(c(.25, .5, 0, -.4, 0, .3), 3, 2)
Psi <- structure(c(0.18, 0.05, -0.05, 0.1, -0.02, 0.06, 0.05, 0.34, -0.25,
                   -0.06, -0.03, 0.29, -0.05, -0.25, 0.24, 0.04, 0.04,
                   -0.12, 0.1, -0.06, 0.04, 0.34, 0, -0.04, -0.02, -0.03,
                   0.04, 0, 0.1, -0.08, 0.06, 0.29, -0.12, -0.04, -0.08,
                   0.51), .Dim = c(6L, 6L))
B <- structure(c(-0.57, 0.17, -0.48, 0.58, 1, 0.86), .Dim = 3:2)
sig <- diag(c(.6, .3)^2)

# generator functions
r_n_marker <- function(id){
  cat(sprintf("r_n_marker: passed id is \%d\n", id))
  # the number of markers is Poisson distributed
  rpois(1, 10) + 1L
}
r_obs_time <- function(id, n_markes){
  cat(sprintf("r_obs_time: passed id is \%d\n", id))
  # the observations times are uniform distributed
  sort(runif(n_markes, 0, 2))
}

# simulate marker
set.seed(1)
U <- draw_U(chol(Psi), NCOL(B))
sim_marker(B = B, U = U, sigma_chol = chol(sig), r_n_marker = r_n_marker,
           r_obs_time = r_obs_time, m_func = m_func, g_func = g_func,
           offset = NULL, id = 1L)

}
\seealso{
\code{\link{draw_U}}, \code{\link{eval_marker}}
}
