# Tests for multi-parameter derivative functions: gradient, hessian, jacobian.
#
# Each test defines a function with known analytical derivatives and
# verifies the AD results.
#
# IMPORTANT: Functions must keep dual parameters as duals throughout.
# Since base R's sum() doesn't dispatch on dual objects, compute sums
# algebraically (e.g., n*mu^2 - 2*mu*sum(data)) rather than element-wise
# (sum((data - mu)^2)).

tol_deriv <- 1e-6

# -- Normal distribution (mu only) --------------------------------------------

test_that("Normal (mu only): gradient and hessian", {
  data <- c(1, 2, 3, 4, 5)
  n <- length(data)
  sum_x <- sum(data)
  sum_x2 <- sum(data^2)

  # Log-likelihood keeping mu as dual:
  # -n/2*log(2pi) - 0.5 * (sum(x^2) - 2*mu*sum(x) + n*mu^2)
  f <- function(x) {
    mu <- x[1]
    -n/2 * log(2 * pi) - 0.5 * (sum_x2 - 2 * mu * sum_x + n * mu^2)
  }

  theta <- c(2)

  g <- gradient(f, theta)
  # Analytical: d/dmu = sum(x) - n*mu = 15 - 10 = 5
  expect_equal(g[1], sum_x - n * theta[1], tolerance = tol_deriv)

  H <- hessian(f, theta)
  # Analytical: d^2/dmu^2 = -n
  expect_equal(H[1,1], -n, tolerance = tol_deriv)

  # observed_information is just -hessian
  expect_equal(-H[1,1], n, tolerance = tol_deriv)
})

# -- Normal distribution (mu and sigma) ----------------------------------------

test_that("Normal (mu, sigma): gradient and hessian", {
  data <- c(1, 2, 3, 4, 5)
  n <- length(data)
  sum_x <- sum(data)
  sum_x2 <- sum(data^2)

  # sum((x-mu)^2) = sum(x^2) - 2*mu*sum(x) + n*mu^2
  f <- function(x) {
    mu <- x[1]
    sigma <- x[2]
    ss <- sum_x2 - 2 * mu * sum_x + n * mu^2  # sum of squared residuals
    -n/2 * log(2 * pi) - n * log(sigma) - 0.5 * ss / sigma^2
  }

  mu0 <- 3
  sigma0 <- sqrt(2)
  theta <- c(mu0, sigma0)

  g <- gradient(f, theta)
  H <- hessian(f, theta)

  # Analytical gradient:
  ss0 <- sum((data - mu0)^2)
  expected_g1 <- sum(data - mu0) / sigma0^2
  expected_g2 <- -n / sigma0 + ss0 / sigma0^3
  expect_equal(g[1], expected_g1, tolerance = tol_deriv)
  expect_equal(g[2], expected_g2, tolerance = tol_deriv)

  # Analytical Hessian:
  expect_equal(H[1,1], -n/sigma0^2, tolerance = tol_deriv)
  expect_equal(H[1,2], -2*sum(data-mu0)/sigma0^3, tolerance = tol_deriv)
  expect_equal(H[2,1], H[1,2], tolerance = 1e-10)
  expect_equal(H[2,2], n/sigma0^2 - 3*ss0/sigma0^4, tolerance = tol_deriv)
})

# -- Poisson distribution (lambda) --------------------------------------------

test_that("Poisson (lambda): gradient and hessian", {
  data <- c(1, 3, 2, 0, 4, 2)
  n <- length(data)
  sum_x <- sum(data)

  f <- function(x) {
    lambda <- x[1]
    sum_x * log(lambda) - n * lambda - sum(lfactorial(data))
  }

  lambda0 <- 2
  theta <- c(lambda0)

  g <- gradient(f, theta)
  H <- hessian(f, theta)

  expect_equal(g[1], sum_x/lambda0 - n, tolerance = tol_deriv)
  expect_equal(H[1,1], -sum_x/lambda0^2, tolerance = tol_deriv)
})

# -- Gamma distribution (shape, rate) -----------------------------------------

test_that("Gamma (shape, rate): gradient matches numerical", {
  data <- c(1.5, 2.1, 0.8, 3.2, 1.9)
  n <- length(data)
  sum_x <- sum(data)
  sum_log_x <- sum(log(data))

  f <- function(x) {
    alpha <- x[1]
    beta <- x[2]
    n * alpha * log(beta) - n * lgamma(alpha) +
      (alpha - 1) * sum_log_x - beta * sum_x
  }

  theta <- c(2, 1.5)

  g <- gradient(f, theta)
  H <- hessian(f, theta)

  # Compare against numerical derivatives of the same function
  num_f <- function(t) {
    n * t[1] * log(t[2]) - n * lgamma(t[1]) +
      (t[1] - 1) * sum_log_x - t[2] * sum_x
  }
  num_g <- numerical_gradient(num_f, theta)
  num_H <- numerical_hessian(num_f, theta)

  expect_equal(g, num_g, tolerance = tol_deriv)
  expect_equal(H, num_H, tolerance = 1e-4)
})

# -- jacobian -----------------------------------------------------------------

test_that("jacobian of vector-valued function", {
  # f: R^2 -> R^2, f(a, b) = (a*b, a^2 + b)
  # J = [[b, a], [2a, 1]]
  f <- function(x) {
    a <- x[1]; b <- x[2]
    list(a * b, a^2 + b)
  }

  theta <- c(3, 4)
  J <- jacobian(f, theta)

  expect_equal(dim(J), c(2, 2))
  expect_equal(J[1, 1], 4, tolerance = 1e-10)   # df1/da = b = 4
  expect_equal(J[1, 2], 3, tolerance = 1e-10)   # df1/db = a = 3
  expect_equal(J[2, 1], 6, tolerance = 1e-10)   # df2/da = 2a = 6
  expect_equal(J[2, 2], 1, tolerance = 1e-10)   # df2/db = 1
})

test_that("jacobian of R^2 -> R^3 function", {
  # f(a, b) = (a*b, a^2, sin(b))
  # J = [[b, a], [2a, 0], [0, cos(b)]]
  f <- function(x) {
    a <- x[1]; b <- x[2]
    list(a * b, a^2, sin(b))
  }

  a0 <- 2; b0 <- pi/4
  J <- jacobian(f, c(a0, b0))

  expect_equal(dim(J), c(3, 2))
  expect_equal(J[1, 1], b0, tolerance = 1e-10)
  expect_equal(J[1, 2], a0, tolerance = 1e-10)
  expect_equal(J[2, 1], 2 * a0, tolerance = 1e-10)
  expect_equal(J[2, 2], 0, tolerance = 1e-10)
  expect_equal(J[3, 1], 0, tolerance = 1e-10)
  expect_equal(J[3, 2], cos(b0), tolerance = 1e-10)
})

test_that("jacobian of scalar function returns gradient vector", {
  f <- function(x) x[1]^2 + x[2]^2

  J <- jacobian(f, c(3, 4))

  # D returns a vector (not 1xp matrix) for scalar-valued f
  expect_equal(length(J), 2)
  expect_equal(J[1], 6, tolerance = 1e-10)
  expect_equal(J[2], 8, tolerance = 1e-10)
})

test_that("jacobian matches numerical for non-linear R^3 -> R^2", {
  # f(a, b, c) = (exp(a*b), sin(b + c))
  f <- function(x) {
    a <- x[1]; b <- x[2]; cc <- x[3]
    list(exp(a * b), sin(b + cc))
  }

  theta <- c(0.5, 1.0, 0.3)
  J <- jacobian(f, theta)

  # Numerical Jacobian via central differences
  eps <- 1e-7
  p <- length(theta)
  J_num <- matrix(0, nrow = 2, ncol = p)
  f_num <- function(t) c(exp(t[1] * t[2]), sin(t[2] + t[3]))
  for (i in seq_len(p)) {
    tp <- tm <- theta
    tp[i] <- tp[i] + eps
    tm[i] <- tm[i] - eps
    J_num[, i] <- (f_num(tp) - f_num(tm)) / (2 * eps)
  }

  expect_equal(J, J_num, tolerance = 1e-6)
})

# -- Works with optim-style function (vector indexing) -------------------------

test_that("gradient and hessian work with optim-style functions", {
  data <- c(2.1, 3.5, 2.8, 4.0, 3.2)
  n <- length(data)
  sum_x <- sum(data)
  sum_x2 <- sum(data^2)

  f <- function(x) {
    mu <- x[1]
    -0.5 * (sum_x2 - 2 * mu * sum_x + n * mu^2)
  }

  g <- gradient(f, c(3))
  expect_true(is.numeric(g))
  expect_equal(length(g), 1)

  H <- hessian(f, c(3))
  expect_true(is.matrix(H))
  expect_equal(dim(H), c(1, 1))
  expect_equal(H[1,1], -n, tolerance = tol_deriv)
})

# -- Multi-parameter indexing --------------------------------------------------

test_that("multi-parameter function with vector indexing", {
  f <- function(x) {
    a <- x[1]
    b <- x[2]
    -(a - 1)^2 - 2 * (b - 3)^2
  }

  theta <- c(0, 0)
  g <- gradient(f, theta)
  # d/da = -2*(a-1) = 2 at a=0
  # d/db = -4*(b-3) = 12 at b=0
  expect_equal(g[1], 2, tolerance = tol_deriv)
  expect_equal(g[2], 12, tolerance = tol_deriv)

  H <- hessian(f, theta)
  expect_equal(H[1,1], -2, tolerance = tol_deriv)
  expect_equal(H[2,2], -4, tolerance = tol_deriv)
  expect_equal(H[1,2], 0, tolerance = tol_deriv)
  expect_equal(H[2,1], 0, tolerance = tol_deriv)
})

# -- Exponential family --------------------------------------------------------

test_that("Exponential distribution: gradient and hessian", {
  data <- c(0.5, 1.2, 0.3, 2.1, 0.8)
  n <- length(data)
  sum_x <- sum(data)

  f <- function(x) {
    rate <- x[1]
    n * log(rate) - rate * sum_x
  }

  rate0 <- 1.5
  theta <- c(rate0)

  g <- gradient(f, theta)
  H <- hessian(f, theta)

  expect_equal(g[1], n/rate0 - sum_x, tolerance = tol_deriv)
  expect_equal(H[1,1], -n/rate0^2, tolerance = tol_deriv)
})

# == Vector-gradient optimization tests (0.4.0) ================================

# -- gradient() matches numerical gradient for all models -----------------------

test_that("1-pass gradient matches numerical_gradient: Normal(mu,sigma)", {
  data <- c(1, 2, 3, 4, 5)
  n <- length(data)
  sum_x <- sum(data)
  sum_x2 <- sum(data^2)

  f <- function(x) {
    mu <- x[1]
    sigma <- x[2]
    ss <- sum_x2 - 2 * mu * sum_x + n * mu^2
    -n/2 * log(2 * pi) - n * log(sigma) - 0.5 * ss / sigma^2
  }

  theta <- c(3, sqrt(2))
  g <- gradient(f, theta)
  num_g <- numerical_gradient(f, theta)
  expect_equal(g, num_g, tolerance = tol_deriv)
})

test_that("1-pass gradient matches numerical_gradient: Poisson(lambda)", {
  data <- c(1, 3, 2, 0, 4, 2)
  n <- length(data)
  sum_x <- sum(data)

  f <- function(x) {
    lambda <- x[1]
    sum_x * log(lambda) - n * lambda - sum(lfactorial(data))
  }

  theta <- c(2)
  g <- gradient(f, theta)
  num_g <- numerical_gradient(f, theta)
  expect_equal(g, num_g, tolerance = tol_deriv)
})

test_that("1-pass gradient matches numerical_gradient: Gamma(shape,rate)", {
  data <- c(1.5, 2.1, 0.8, 3.2, 1.9)
  n <- length(data)
  sum_x <- sum(data)
  sum_log_x <- sum(log(data))

  f <- function(x) {
    alpha <- x[1]
    beta <- x[2]
    n * alpha * log(beta) - n * lgamma(alpha) +
      (alpha - 1) * sum_log_x - beta * sum_x
  }

  theta <- c(2, 1.5)
  g <- gradient(f, theta)
  num_g <- numerical_gradient(f, theta)
  expect_equal(g, num_g, tolerance = tol_deriv)
})

# -- hessian() symmetry and numerical accuracy --------------------------------

test_that("p-pass hessian is symmetric: Normal(mu,sigma)", {
  data <- c(1, 2, 3, 4, 5)
  n <- length(data)
  sum_x <- sum(data)
  sum_x2 <- sum(data^2)

  f <- function(x) {
    mu <- x[1]
    sigma <- x[2]
    ss <- sum_x2 - 2 * mu * sum_x + n * mu^2
    -n/2 * log(2 * pi) - n * log(sigma) - 0.5 * ss / sigma^2
  }

  theta <- c(3, sqrt(2))
  H <- hessian(f, theta)

  expect_equal(H[1,2], H[2,1], tolerance = 1e-10)
})

test_that("p-pass hessian matches numerical_hessian: Normal(mu,sigma)", {
  data <- c(1, 2, 3, 4, 5)
  n <- length(data)
  sum_x <- sum(data)
  sum_x2 <- sum(data^2)

  f <- function(x) {
    mu <- x[1]
    sigma <- x[2]
    ss <- sum_x2 - 2 * mu * sum_x + n * mu^2
    -n/2 * log(2 * pi) - n * log(sigma) - 0.5 * ss / sigma^2
  }

  theta <- c(3, sqrt(2))
  H <- hessian(f, theta)
  num_H <- numerical_hessian(f, theta)
  expect_equal(H, num_H, tolerance = 1e-4)
})

test_that("p-pass hessian matches numerical_hessian: multi-param quadratic", {
  f <- function(x) {
    a <- x[1]
    b <- x[2]
    -(a - 1)^2 - 2 * (b - 3)^2 - 0.5 * a * b
  }

  theta <- c(0, 0)
  H <- hessian(f, theta)
  num_H <- numerical_hessian(f, theta)

  expect_equal(H, num_H, tolerance = 1e-4)
  # Cross-derivative should be symmetric
  expect_equal(H[1,2], H[2,1], tolerance = 1e-10)
  # Diagonal entries: d2f/da2 = -2, d2f/db2 = -4
  expect_equal(H[1,1], -2, tolerance = tol_deriv)
  expect_equal(H[2,2], -4, tolerance = tol_deriv)
  # Off-diagonal: d2f/dadb = -0.5
  expect_equal(H[1,2], -0.5, tolerance = tol_deriv)
})

test_that("p-pass hessian matches numerical: 3-parameter model", {
  f <- function(x) {
    a <- x[1]; b <- x[2]; c <- x[3]
    -(a - 1)^2 - 2 * (b - 2)^2 - 3 * (c - 3)^2 + a * b - 0.5 * b * c
  }

  theta <- c(1, 2, 3)
  H <- hessian(f, theta)
  num_H <- numerical_hessian(f, theta)

  expect_equal(H, num_H, tolerance = 1e-4)
  # Verify symmetry
  expect_equal(H[1,2], H[2,1], tolerance = 1e-10)
  expect_equal(H[1,3], H[3,1], tolerance = 1e-10)
  expect_equal(H[2,3], H[3,2], tolerance = 1e-10)
})

# -- Edge case: constant function returns zeros --------------------------------

test_that("gradient and hessian handle constant function", {
  f <- function(x) 42

  theta <- c(1, 2)

  g <- gradient(f, theta)
  expect_equal(g, c(0, 0))

  H <- hessian(f, theta)
  expect_equal(H, matrix(0, 2, 2))
})

# -- Edge case: single parameter (p=1) ----------------------------------------

test_that("gradient and hessian work with single parameter", {
  f <- function(x) {
    a <- x[1]
    -0.5 * (a - 3)^2
  }

  theta <- c(1)
  g <- gradient(f, theta)
  H <- hessian(f, theta)

  expect_equal(g[1], 2, tolerance = tol_deriv)  # d/da = -(a-3) = 2 at a=1
  expect_equal(H[1,1], -1, tolerance = tol_deriv)
})

# -- -hessian() gives observed information -------------------------------------

test_that("-hessian() gives observed information", {
  data <- c(1, 2, 3, 4, 5)
  n <- length(data)
  sum_x <- sum(data)
  sum_x2 <- sum(data^2)

  f <- function(x) {
    mu <- x[1]
    sigma <- x[2]
    ss <- sum_x2 - 2 * mu * sum_x + n * mu^2
    -n/2 * log(2 * pi) - n * log(sigma) - 0.5 * ss / sigma^2
  }

  theta <- c(3, sqrt(2))
  I <- -hessian(f, theta)
  H <- hessian(f, theta)

  expect_equal(I, -H, tolerance = 1e-10)
})

# == D operator tests ==========================================================

tol_D1 <- 1e-10   # first-order (exact AD)
tol_D2 <- 1e-10   # second-order
tol_D3 <- 1e-10   # third-order (AD still exact; only numerical comparison is loose)
tol_num <- 1e-5   # for comparing against numerical finite differences

# -- D returns a function ------------------------------------------------------

test_that("D(f) returns a function", {
  f <- function(x) x[1]^2 + x[2]^2
  g <- D(f)
  expect_true(is.function(g))
  result <- g(c(3, 4))
  expect_equal(result, c(6, 8), tolerance = tol_D1)
})

test_that("D(f, order=2) returns a function when x is NULL", {
  f <- function(x) x[1]^2 + x[2]^2
  g <- D(f, order = 2)
  expect_true(is.function(g))
  result <- g(c(3, 4))
  expect_equal(result, diag(2, 2), tolerance = tol_D2)
})

test_that("D rejects order < 1", {
  expect_error(D(sin, order = 0), "positive integer")
})

# -- First-order, scalar f: matches gradient -----------------------------------

test_that("D(f, x) matches gradient for quadratic", {
  f <- function(x) x[1]^2 * x[2] + x[2]^3
  x <- c(3, 4)
  expect_equal(D(f, x), gradient(f, x), tolerance = tol_D1)
})

test_that("D(f, x) matches gradient for exponential", {
  f <- function(x) exp(x[1] + x[2])
  x <- c(1, 2)
  expect_equal(D(f, x), gradient(f, x), tolerance = tol_D1)
})

test_that("D(f, x) matches gradient for trigonometric", {
  f <- function(x) sin(x[1]) * cos(x[2])
  x <- c(pi/3, pi/4)
  expect_equal(D(f, x), gradient(f, x), tolerance = tol_D1)
})

test_that("D(f, x) matches numerical gradient", {
  f <- function(x) x[1]^2 * x[2] + x[2]^3
  x <- c(3, 4)
  expect_equal(D(f, x), numerical_gradient(f, x), tolerance = tol_num)
})

# -- First-order, vector f: matches old jacobian results -----------------------

test_that("D(f, x) produces Jacobian for R^2 -> R^2", {
  f <- function(x) {
    a <- x[1]; b <- x[2]
    list(a * b, a^2 + b)
  }
  J <- D(f, c(3, 4))
  expect_equal(dim(J), c(2, 2))
  expect_equal(J[1, 1], 4, tolerance = tol_D1)   # df1/da = b = 4
  expect_equal(J[1, 2], 3, tolerance = tol_D1)   # df1/db = a = 3
  expect_equal(J[2, 1], 6, tolerance = tol_D1)   # df2/da = 2a = 6
  expect_equal(J[2, 2], 1, tolerance = tol_D1)   # df2/db = 1
})

test_that("D(f, x) produces Jacobian for R^2 -> R^3", {
  f <- function(x) {
    a <- x[1]; b <- x[2]
    list(a * b, a^2, sin(b))
  }
  a0 <- 2; b0 <- pi/4
  J <- D(f, c(a0, b0))
  expect_equal(dim(J), c(3, 2))
  expect_equal(J[1, 1], b0, tolerance = tol_D1)
  expect_equal(J[1, 2], a0, tolerance = tol_D1)
  expect_equal(J[2, 1], 2 * a0, tolerance = tol_D1)
  expect_equal(J[2, 2], 0, tolerance = tol_D1)
  expect_equal(J[3, 1], 0, tolerance = tol_D1)
  expect_equal(J[3, 2], cos(b0), tolerance = tol_D1)
})

test_that("D(f, x) for scalar f returns vector (not matrix)", {
  f <- function(x) x[1]^2 + x[2]^2
  result <- D(f, c(3, 4))
  expect_true(is.numeric(result))
  expect_equal(length(result), 2)
  expect_null(dim(result))  # vector, not 1-row matrix
  expect_equal(result, c(6, 8), tolerance = tol_D1)
})

# -- Second-order, scalar f: matches hessian -----------------------------------

test_that("D(f, x, order=2) matches hessian for quadratic", {
  f <- function(x) -(x[1] - 3)^2 - (x[2] - 5)^2
  x <- c(1, 2)
  expect_equal(D(f, x, order = 2), hessian(f, x), tolerance = tol_D2)
})

test_that("D(f, x, order=2) matches hessian for cross-term", {
  f <- function(x) x[1]^2 * x[2] + x[2]^3
  x <- c(3, 4)
  H <- D(f, x, order = 2)
  expect_equal(dim(H), c(2, 2))
  # d2f/da2 = 2*b = 8
  expect_equal(H[1, 1], 2 * x[2], tolerance = tol_D2)
  # d2f/dadb = 2*a = 6
  expect_equal(H[1, 2], 2 * x[1], tolerance = tol_D2)
  expect_equal(H[2, 1], 2 * x[1], tolerance = tol_D2)  # symmetry
  # d2f/db2 = 6*b = 24
  expect_equal(H[2, 2], 6 * x[2], tolerance = tol_D2)
})

test_that("D(f, x, order=2) matches numerical Hessian", {
  f <- function(x) x[1]^2 * x[2] + x[2]^3
  x <- c(3, 4)
  expect_equal(D(f, x, order = 2), numerical_hessian(f, x), tolerance = tol_num)
})

# -- Second-order via composition: D(D(f)) == D(f, order=2) -------------------

test_that("D(D(f))(x) equals D(f, x, order=2)", {
  f <- function(x) x[1]^2 * x[2]
  x <- c(3, 4)
  Df <- D(f)
  DDf <- D(Df)
  expect_equal(DDf(x), D(f, x, order = 2), tolerance = tol_D2)
})

# -- Third-order ---------------------------------------------------------------

test_that("D(f, x, order=3) for x[1]^4 at x=c(2,1)", {
  f <- function(x) x[1]^4
  x <- c(2, 1)
  T3 <- D(f, x, order = 3)
  # f = x1^4, only x1 derivatives are non-zero
  # d3f/dx1^3 = 24*x1 = 48
  # All other entries are 0 (no x2 dependence)
  expect_equal(dim(T3), c(2, 2, 2))
  expect_equal(T3[1, 1, 1], 48, tolerance = tol_D3)  # d3f/dx1^3
  # All entries involving x2 should be zero
  expect_equal(T3[2, 1, 1], 0, tolerance = tol_D3)
  expect_equal(T3[1, 2, 1], 0, tolerance = tol_D3)
  expect_equal(T3[1, 1, 2], 0, tolerance = tol_D3)
  expect_equal(T3[2, 2, 1], 0, tolerance = tol_D3)
  expect_equal(T3[2, 1, 2], 0, tolerance = tol_D3)
  expect_equal(T3[1, 2, 2], 0, tolerance = tol_D3)
  expect_equal(T3[2, 2, 2], 0, tolerance = tol_D3)
})

test_that("D(f, x, order=3) for x[1]^2 * x[2] gives correct tensor", {
  f <- function(x) x[1]^2 * x[2]
  x <- c(3, 4)
  T3 <- D(f, x, order = 3)
  expect_equal(dim(T3), c(2, 2, 2))
  # d3f/dx1^2 dx2 = 2 (in all permutations of indices)
  # d3f/dx1^3 = 0
  # d3f/dx2^3 = 0
  expect_equal(T3[1, 1, 1], 0, tolerance = tol_D3)
  expect_equal(T3[1, 1, 2], 2, tolerance = tol_D3)
  expect_equal(T3[1, 2, 1], 2, tolerance = tol_D3)
  expect_equal(T3[2, 1, 1], 2, tolerance = tol_D3)
  expect_equal(T3[2, 2, 2], 0, tolerance = tol_D3)
  expect_equal(T3[2, 2, 1], 0, tolerance = tol_D3)
  expect_equal(T3[2, 1, 2], 0, tolerance = tol_D3)
  expect_equal(T3[1, 2, 2], 0, tolerance = tol_D3)
})

# -- Fourth-order --------------------------------------------------------------

test_that("D(f, x, order=4) for exp(x[1]) at x=c(1)", {
  f <- function(x) exp(x[1])
  x <- c(1)
  T4 <- D(f, x, order = 4)
  # exp has all derivatives equal to exp(x)
  # n=1: 4th-order tensor is (1,1,1,1) array
  expect_equal(dim(T4), c(1, 1, 1, 1))
  expect_equal(T4[1, 1, 1, 1], exp(1), tolerance = tol_D3)
})

test_that("D(f, x, order=4) for x[1]^5 at x=c(2)", {
  f <- function(x) x[1]^5
  x <- c(2)
  T4 <- D(f, x, order = 4)
  # f^(4) = 120*x = 240
  expect_equal(dim(T4), c(1, 1, 1, 1))
  expect_equal(T4[1, 1, 1, 1], 240, tolerance = tol_D3)
})

# -- Second-order, vector f: produces (m,n,n) tensor ---------------------------

test_that("D(f, x, order=2) for f: R^2 -> R^2 gives (2,2,2) array", {
  f <- function(x) list(x[1]^2 * x[2], x[1] + x[2]^3)
  x <- c(2, 3)
  T2 <- D(f, x, order = 2)
  expect_equal(dim(T2), c(2, 2, 2))
  # Component 1: x1^2 * x2
  #   d2/dx1^2 = 2*x2 = 6
  #   d2/dx1dx2 = 2*x1 = 4
  #   d2/dx2^2 = 0
  expect_equal(T2[1, 1, 1], 6, tolerance = tol_D2)
  expect_equal(T2[1, 1, 2], 4, tolerance = tol_D2)
  expect_equal(T2[1, 2, 1], 4, tolerance = tol_D2)
  expect_equal(T2[1, 2, 2], 0, tolerance = tol_D2)
  # Component 2: x1 + x2^3
  #   d2/dx1^2 = 0
  #   d2/dx1dx2 = 0
  #   d2/dx2^2 = 6*x2 = 18
  expect_equal(T2[2, 1, 1], 0, tolerance = tol_D2)
  expect_equal(T2[2, 1, 2], 0, tolerance = tol_D2)
  expect_equal(T2[2, 2, 1], 0, tolerance = tol_D2)
  expect_equal(T2[2, 2, 2], 18, tolerance = tol_D2)
})

# -- Edge cases ----------------------------------------------------------------

test_that("D of constant function gives zero tensor", {
  f <- function(x) 42
  x <- c(1, 2)
  expect_equal(D(f, x), c(0, 0))
  expect_equal(D(f, x, order = 2), matrix(0, 2, 2))
})

test_that("D of linear function gives constant Jacobian, zero Hessian", {
  f <- function(x) 3 * x[1] + 5 * x[2]
  x <- c(1, 2)
  expect_equal(D(f, x), c(3, 5), tolerance = tol_D1)
  expect_equal(D(f, x, order = 2), matrix(0, 2, 2))
})

test_that("D works with n=1 (scalar input)", {
  f <- function(x) x[1]^3
  expect_equal(D(f, c(2)), 12, tolerance = tol_D1)        # 3*x^2
  # n=1 higher-order: returns (1,1) matrix, (1,1,1) array, etc.
  expect_equal(D(f, c(2), order = 2)[1, 1], 12, tolerance = tol_D2)  # 6*x
  expect_equal(D(f, c(2), order = 3)[1, 1, 1], 6, tolerance = tol_D3)   # 6
})

test_that("D works with n=1 for exp", {
  f <- function(x) exp(x[1])
  expect_equal(D(f, c(1)), exp(1), tolerance = tol_D1)
  expect_equal(D(f, c(1), order = 2)[1, 1], exp(1), tolerance = tol_D2)
})

# -- Numerical verification for gradient and hessian via D ---------------------

test_that("gradient via D matches numerical for Normal(mu,sigma)", {
  data <- c(1, 2, 3, 4, 5)
  n <- length(data)
  sum_x <- sum(data)
  sum_x2 <- sum(data^2)

  f <- function(x) {
    mu <- x[1]; sigma <- x[2]
    ss <- sum_x2 - 2 * mu * sum_x + n * mu^2
    -n/2 * log(2 * pi) - n * log(sigma) - 0.5 * ss / sigma^2
  }

  theta <- c(3, sqrt(2))
  expect_equal(gradient(f, theta), numerical_gradient(f, theta), tolerance = tol_num)
})

test_that("hessian via D matches numerical for Normal(mu,sigma)", {
  data <- c(1, 2, 3, 4, 5)
  n <- length(data)
  sum_x <- sum(data)
  sum_x2 <- sum(data^2)

  f <- function(x) {
    mu <- x[1]; sigma <- x[2]
    ss <- sum_x2 - 2 * mu * sum_x + n * mu^2
    -n/2 * log(2 * pi) - n * log(sigma) - 0.5 * ss / sigma^2
  }

  theta <- c(3, sqrt(2))
  expect_equal(hessian(f, theta), numerical_hessian(f, theta), tolerance = 1e-4)
})

test_that("hessian via D is symmetric", {
  f <- function(x) x[1]^2 * x[2] + exp(x[1] * x[2])
  x <- c(1, 0.5)
  H <- hessian(f, x)
  expect_equal(H[1, 2], H[2, 1], tolerance = 1e-12)
})
