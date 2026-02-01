# Higher-order derivatives via nested dual numbers
#
# Arbitrary-order derivatives are computed by recursively nesting duals.
# At each nesting level, the dual tracks d/dx of the level below.
#
# For order n, the seeding structure is:
#   n=0: just x (plain numeric)
#   n=1: dual(x, 1)
#   n=2: dual(dual(x, 1), dual(1, 0))
#   n=3: dual(dual(dual(x, 1), dual(1, 0)), dual(dual(1, 0), dual(0, 0)))
#
# Extraction: deriv_n(result, k) applies @deriv k times, then walks
# @value the rest of the way down to get the k-th derivative.

# =============================================================================
# Arbitrary-order constructors
# =============================================================================

#' Create a dual seeded for n-th order differentiation
#'
#' Recursively nests dual numbers to enable exact computation of
#' derivatives up to order \code{n}. The variable is seeded so that
#' after evaluating a function \code{f}, the k-th derivative can be
#' extracted with \code{deriv_n(result, k)}.
#'
#' @param x A numeric value at which to differentiate.
#' @param order A positive integer specifying the derivative order.
#' @return A (possibly nested) dual number.
#' @export
#' @examples
#' # 3rd derivative of x^4 at x=2: 4*3*2*x = 24*2 = 48
#' x <- dual_variable_n(2, order = 3)
#' r <- x^4
#' deriv_n(r, 3)  # 48
dual_variable_n <- function(x, order) {
  order <- as.integer(order)
  if (order < 0L) stop("order must be a non-negative integer")
  if (order == 0L) return(x)
  if (order == 1L) return(.dual(x, 1))
  .dual(dual_variable_n(x, order - 1L),
        dual_constant_n(1, order - 1L))
}

#' Create a constant dual for n-th order differentiation
#'
#' Wraps a numeric value as a nested dual with all derivative components
#' zero, representing a constant with respect to the differentiation
#' variable at nesting depth \code{n}.
#'
#' @param x A numeric value.
#' @param order A non-negative integer specifying the nesting depth.
#' @return A (possibly nested) dual number with zero derivatives.
#' @export
#' @examples
#' k <- dual_constant_n(5, order = 3)
#' deriv_n(k, 1)  # 0
#' deriv_n(k, 2)  # 0
#' deriv_n(k, 3)  # 0
dual_constant_n <- function(x, order) {
  order <- as.integer(order)
  if (order < 0L) stop("order must be a non-negative integer")
  if (order == 0L) return(x)
  if (order == 1L) return(.dual(x, 0))
  .dual(dual_constant_n(x, order - 1L),
        dual_constant_n(0, order - 1L))
}

# =============================================================================
# Extraction
# =============================================================================

#' Extract the k-th derivative from a nested dual result
#'
#' After evaluating a function on a dual created by
#' \code{\link{dual_variable_n}}, use \code{deriv_n} to extract any
#' derivative from 0 (the function value) up to the seeded order.
#'
#' @param d A (possibly nested) dual number, or a numeric.
#' @param k A non-negative integer: 0 for the function value, 1 for the
#'   first derivative, etc.
#' @return A numeric value.
#' @export
#' @examples
#' x <- dual_variable_n(1, order = 3)
#' r <- exp(x)
#' deriv_n(r, 0)  # exp(1) = 2.718...
#' deriv_n(r, 1)  # exp(1)
#' deriv_n(r, 2)  # exp(1)
#' deriv_n(r, 3)  # exp(1)
deriv_n <- function(d, k) {
  k <- as.integer(k)
  if (k < 0L) stop("k must be a non-negative integer")
  if (k == 0L) {
    # Walk all the way down through value slots
    while (is(d, "dualr")) d <- d@value
    return(d)
  }
  # Take one deriv, then recurse
  if (!is(d, "dualr")) return(0)
  deriv_n(d@deriv, k - 1L)
}

# =============================================================================
# Convenience evaluator
# =============================================================================

#' Compute a function value and all derivatives up to order n
#'
#' Evaluates \code{f} at a dual variable seeded for order \code{n},
#' returning the function value and all derivatives from 1 to \code{n}.
#'
#' @param f A function of one numeric argument.
#' @param x A numeric value at which to differentiate.
#' @param order A positive integer: the maximum derivative order.
#' @return A named list with components \code{value}, \code{d1},
#'   \code{d2}, ..., \code{d<order>}.
#' @export
#' @examples
#' # All derivatives of sin(x) at x = pi/4
#' differentiate_n(sin, pi/4, order = 4)
#' # $value = sin(pi/4)
#' # $d1 = cos(pi/4)
#' # $d2 = -sin(pi/4)
#' # $d3 = -cos(pi/4)
#' # $d4 = sin(pi/4)
differentiate_n <- function(f, x, order) {
  order <- as.integer(order)
  if (order < 1L) stop("order must be a positive integer")
  result <- f(dual_variable_n(x, order))
  out <- list(value = deriv_n(result, 0L))
  for (k in seq_len(order)) {
    out[[paste0("d", k)]] <- deriv_n(result, k)
  }
  out
}

# =============================================================================
# Deprecated 2nd-order wrappers (backward compatibility)
# =============================================================================

#' Create a second-order dual variable
#'
#' @description
#' \strong{Deprecated.}
#'
#' Use \code{\link{dual_variable_n}(x, 2)} instead.
#'
#' @param x A numeric value.
#' @return A nested dual: \code{dual(dual(x, 1), dual(1, 0))}.
#' @export
#' @examples
#' x <- dual2_variable(2)
#' result <- x^3
#' second_deriv(result)  # 12 = 6*x at x=2
dual2_variable <- function(x) {
  dual_variable_n(x, 2L)
}

#' Create a second-order dual constant
#'
#' @description
#' \strong{Deprecated.}
#'
#' Use \code{\link{dual_constant_n}(x, 2)} instead.
#'
#' @param x A numeric value.
#' @return A nested dual: \code{dual(dual(x, 0), dual(0, 0))}.
#' @export
#' @examples
#' k <- dual2_constant(5)
#' second_deriv(k^2)  # 0 (constant has no derivatives)
dual2_constant <- function(x) {
  dual_constant_n(x, 2L)
}

#' Extract function value from a second-order dual
#'
#' @description
#' \strong{Deprecated.}
#'
#' Use \code{\link{deriv_n}(d, 0)} instead.
#'
#' @param d A nested dual.
#' @return The numeric function value.
#' @export
#' @examples
#' x <- dual2_variable(3)
#' result <- x^2
#' value2(result)  # 9
value2 <- function(d) {
  deriv_n(d, 0L)
}

#' Extract first derivative from a second-order dual
#'
#' @description
#' \strong{Deprecated.}
#'
#' Use \code{\link{deriv_n}(d, 1)} instead.
#'
#' @param d A nested dual.
#' @return The numeric first derivative.
#' @export
#' @examples
#' x <- dual2_variable(3)
#' result <- x^2
#' first_deriv(result)  # 6 (= 2*x at x=3)
first_deriv <- function(d) {
  deriv_n(d, 1L)
}

#' Extract second derivative from a second-order dual
#'
#' @description
#' \strong{Deprecated.}
#'
#' Use \code{\link{deriv_n}(d, 2)} instead.
#'
#' @param d A nested dual.
#' @return The numeric second derivative.
#' @export
#' @examples
#' x <- dual2_variable(3)
#' result <- x^2
#' second_deriv(result)  # 2
second_deriv <- function(d) {
  deriv_n(d, 2L)
}

#' Compute value, first and second derivatives of a function
#'
#' @description
#' \strong{Deprecated.}
#'
#' Use \code{\link{differentiate_n}(f, x, 2)} instead.
#'
#' @param f A function of one argument.
#' @param x A numeric value at which to differentiate.
#' @return A named list with components \code{value}, \code{first}, and
#'   \code{second}.
#' @export
#' @examples
#' differentiate2(sin, pi/4)
#' # $value = sin(pi/4)
#' # $first = cos(pi/4)
#' # $second = -sin(pi/4)
differentiate2 <- function(f, x) {
  result <- f(dual_variable_n(x, 2L))
  list(
    value  = deriv_n(result, 0L),
    first  = deriv_n(result, 1L),
    second = deriv_n(result, 2L)
  )
}
