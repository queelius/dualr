# MLE helper functions
#
# High-level functions for maximum likelihood estimation workflows.
# These wrap forward-mode AD to compute score vectors, Hessians, and
# observed information matrices.
#
# The loglik function signature convention:
#   loglik(theta) where theta supports theta[1], theta[2], ... indexing
#   This works with both numeric vectors and dual_vector objects.

# -- Internal: build a dual_vector seeding parameter i with deriv=1 -----------
# (kept for score_and_hessian, which differentiates a vector-valued function)

.make_dual_vector <- function(theta, seed_index) {
  p <- length(theta)
  duals <- vector("list", p)
  for (j in seq_len(p)) {
    duals[[j]] <- if (j == seed_index) .dual(theta[j], 1) else .dual(theta[j], 0)
  }
  new("dual_vector", duals)
}

# -- Internal: vector-gradient seeding (1-pass score) -------------------------
# Seeds all p parameters with unit-vector gradients: dualr(theta_k, e_k)
# where e_k is the k-th standard basis vector of length p.
# One forward pass then yields the full gradient in result@deriv.

.make_grad_vector <- function(theta) {
  p <- length(theta)
  duals <- vector("list", p)
  for (k in seq_len(p)) {
    e_k <- numeric(p)
    e_k[k] <- 1
    duals[[k]] <- .dual(theta[k], e_k)
  }
  new("dual_vector", duals)
}

# -- Internal: nested vector-gradient seeding (p-pass Hessian) ----------------
# For Hessian row i: seeds parameter k with a nested dual where the inner
# level carries the full vector gradient and the outer level tracks d/d(theta_i).
# After evaluation, result@deriv@deriv gives the entire i-th row of H.

.make_grad2_vector <- function(theta, i) {
  p <- length(theta)
  duals <- vector("list", p)
  for (k in seq_len(p)) {
    e_k <- numeric(p)
    e_k[k] <- 1
    inner <- .dual(theta[k], e_k)
    outer_d <- if (k == i) .dual(1, numeric(p)) else .dual(0, numeric(p))
    duals[[k]] <- .dual(inner, outer_d)
  }
  new("dual_vector", duals)
}

#' Compute the score (gradient) of a log-likelihood
#'
#' Evaluates the gradient of \code{loglik} at \code{theta} using forward-mode
#' AD with vector-valued derivatives. Computes the full gradient in a single
#' forward pass by seeding each parameter with a unit-vector gradient.
#'
#' The log-likelihood function should use \code{theta[1]}, \code{theta[2]},
#' etc. to access parameters. This is compatible with the standard R
#' convention used by \code{optim}.
#'
#' @param loglik A function taking a parameter vector (numeric or dual) and
#'   returning a scalar. Parameters are accessed via \code{theta[i]}.
#' @param theta A numeric vector of parameter values.
#' @return A numeric vector of length \code{p} containing the gradient.
#' @export
#' @examples
#' # Normal log-likelihood for mu
#' ll <- function(theta) {
#'   mu <- theta[1]
#'   -0.5 * (2 - mu)^2
#' }
#' score(ll, 1)  # should be 1 (= 2 - 1)
score <- function(loglik, theta) {
  p <- length(theta)
  dual_theta <- .make_grad_vector(theta)
  result <- loglik(dual_theta)
  if (is(result, "dualr")) {
    d <- result@deriv
    if (is.numeric(d) && length(d) == p) return(d)
  }
  numeric(p)
}

#' Compute the Hessian of a log-likelihood
#'
#' Uses nested duals with vector-gradient inner derivatives to compute the
#' matrix of second partial derivatives. Each of \code{p} forward passes
#' yields an entire row of the Hessian.
#'
#' @param loglik A function taking a parameter vector and returning a scalar.
#' @param theta A numeric vector of parameter values.
#' @return A \code{p x p} numeric matrix (the Hessian).
#' @export
#' @examples
#' # Normal log-likelihood for mu
#' ll <- function(theta) {
#'   mu <- theta[1]
#'   -0.5 * (2 - mu)^2
#' }
#' hessian(ll, 1)  # should be matrix(-1)
hessian <- function(loglik, theta) {
  p <- length(theta)
  H <- matrix(0, nrow = p, ncol = p)

  for (i in seq_len(p)) {
    dual_theta <- .make_grad2_vector(theta, i)
    result <- loglik(dual_theta)
    if (is(result, "dualr") && is(result@deriv, "dualr")) {
      H[i, ] <- result@deriv@deriv
    }
    # else: row stays zero (constant function)
  }

  H
}

#' Compute the observed information matrix
#'
#' Returns the negative Hessian of the log-likelihood at \code{theta}.
#'
#' @param loglik A function taking a parameter vector and returning a scalar.
#' @param theta A numeric vector of parameter values.
#' @return A \code{p x p} numeric matrix.
#' @examples
#' ll <- function(theta) -(theta[1] - 3)^2
#' observed_information(ll, 2)  # matrix(2)
#' @export
observed_information <- function(loglik, theta) {
  -hessian(loglik, theta)
}

#' Compute score and Hessian from an analytical score function
#'
#' When the score function is provided analytically, differentiates it
#' to obtain the Hessian. This is useful when the score is simpler than
#' the log-likelihood (e.g., after algebraic simplification).
#'
#' The score function should take a parameter vector (supporting
#' \code{theta[i]} indexing) and return a list of scalar values (one
#' per parameter).
#'
#' @param score_fn A function taking a parameter vector and returning a
#'   list of scalar values (one per component of the score).
#' @param theta A numeric vector of parameter values.
#' @return A list with components:
#'   \describe{
#'     \item{score}{Numeric vector of score values.}
#'     \item{hessian}{\code{p x p} numeric Hessian matrix (Jacobian of score).}
#'   }
#' @examples
#' # Analytical score of -(theta - 3)^2
#' s_fn <- function(theta) list(-2 * (theta[1] - 3))
#' result <- score_and_hessian(s_fn, 2)
#' result$score    # 2
#' result$hessian  # matrix(-2)
#' @export
score_and_hessian <- function(score_fn, theta) {
  p <- length(theta)
  s <- numeric(p)
  H <- matrix(0, nrow = p, ncol = p)

  for (i in seq_len(p)) {
    dual_theta <- .make_dual_vector(theta, i)
    result <- score_fn(dual_theta)

    # Extract score values (only need to do this on the first pass)
    if (i == 1L) {
      for (k in seq_len(p)) {
        rk <- result[[k]]
        s[k] <- if (is(rk, "dualr")) rk@value else rk
      }
    }

    # The Jacobian column i: d(score_k)/d(theta_i) for each k
    for (k in seq_len(p)) {
      rk <- result[[k]]
      H[k, i] <- if (is(rk, "dualr")) rk@deriv else 0
    }
  }

  list(score = s, hessian = H)
}
