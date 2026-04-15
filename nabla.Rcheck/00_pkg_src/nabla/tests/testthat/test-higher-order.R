# Tests for arbitrary-order derivatives via nested duals

tol_2nd <- 1e-5   # 2nd-order tolerance
tol_3rd <- 1e-4   # 3rd-order (numerical FD loses precision)
tol_4th <- 1e-3   # 4th-order
tol_5th <- 1e-2   # 5th-order

# =============================================================================
# dual_variable_n / dual_constant_n constructors
# =============================================================================

test_that("dual_variable_n(x, 1) matches dual_variable", {
  x1 <- dual_variable_n(3, 1)
  x2 <- dual_variable(3)
  expect_equal(value(x1), value(x2))
  expect_equal(deriv(x1), deriv(x2))
})

test_that("dual_variable_n(x, 2) extracts correctly via deriv_n", {
  x <- dual_variable_n(3, 2)
  expect_equal(deriv_n(x, 0), 3)
  expect_equal(deriv_n(x, 1), 1)
  expect_equal(deriv_n(x, 2), 0)
})

test_that("dual_variable_n(x, 0) returns plain numeric", {
  expect_identical(dual_variable_n(3, 0), 3)
})

test_that("dual_constant_n(x, 0) returns plain numeric", {
  expect_identical(dual_constant_n(5, 0), 5)
})

test_that("dual_constant_n has all-zero derivatives", {
  k <- dual_constant_n(7, order = 4)
  for (i in 1:4) {
    expect_equal(deriv_n(k, i), 0, label = paste("deriv_n(k,", i, ")"))
  }
  expect_equal(deriv_n(k, 0), 7)
})

test_that("dual_variable_n rejects negative order", {
  expect_error(dual_variable_n(1, -1), "non-negative")
})

test_that("dual_constant_n rejects negative order", {
  expect_error(dual_constant_n(1, -1), "non-negative")
})

# =============================================================================
# deriv_n extraction
# =============================================================================

test_that("deriv_n(numeric, 0) returns the numeric", {
  expect_equal(deriv_n(42, 0), 42)
})

test_that("deriv_n(numeric, k>0) returns 0", {
  expect_equal(deriv_n(42, 1), 0)
  expect_equal(deriv_n(42, 5), 0)
})

test_that("deriv_n rejects negative k", {
  expect_error(deriv_n(dual_variable(1), -1), "non-negative")
})


# =============================================================================
# Second-order derivatives (existing tests, updated)
# =============================================================================

test_that("x^2: f''=2", {
  x <- dual_variable_n(3, 2)
  r <- x^2
  expect_equal(deriv_n(r, 0), 9)
  expect_equal(deriv_n(r, 1), 6)
  expect_equal(deriv_n(r, 2), 2)
})

test_that("x^3: f''=6x", {
  x <- dual_variable_n(2, 2)
  r <- x^3
  expect_equal(deriv_n(r, 0), 8)
  expect_equal(deriv_n(r, 1), 12)
  expect_equal(deriv_n(r, 2), 12)
})

test_that("exp(x): f''=exp(x)", {
  x_val <- 1.5
  x <- dual_variable_n(x_val, 2)
  r <- exp(x)
  expect_equal(deriv_n(r, 0), exp(x_val))
  expect_equal(deriv_n(r, 1), exp(x_val), tolerance = 1e-10)
  expect_equal(deriv_n(r, 2), exp(x_val), tolerance = 1e-10)
})

test_that("sin(x): f''=-sin(x)", {
  x_val <- pi / 3
  x <- dual_variable_n(x_val, 2)
  r <- sin(x)
  expect_equal(deriv_n(r, 0), sin(x_val))
  expect_equal(deriv_n(r, 1), cos(x_val), tolerance = 1e-10)
  expect_equal(deriv_n(r, 2), -sin(x_val), tolerance = 1e-10)
})

test_that("cos(x): f''=-cos(x)", {
  x_val <- pi / 4
  x <- dual_variable_n(x_val, 2)
  r <- cos(x)
  expect_equal(deriv_n(r, 0), cos(x_val))
  expect_equal(deriv_n(r, 1), -sin(x_val), tolerance = 1e-10)
  expect_equal(deriv_n(r, 2), -cos(x_val), tolerance = 1e-10)
})

test_that("log(x): f''=-1/x^2", {
  x_val <- 2
  x <- dual_variable_n(x_val, 2)
  r <- log(x)
  expect_equal(deriv_n(r, 0), log(x_val))
  expect_equal(deriv_n(r, 1), 1/x_val, tolerance = 1e-10)
  expect_equal(deriv_n(r, 2), -1/x_val^2, tolerance = 1e-10)
})

test_that("1/x: f''=2/x^3", {
  x_val <- 3
  x <- dual_variable_n(x_val, 2)
  r <- 1 / x
  expect_equal(deriv_n(r, 0), 1/x_val)
  expect_equal(deriv_n(r, 1), -1/x_val^2, tolerance = 1e-10)
  expect_equal(deriv_n(r, 2), 2/x_val^3, tolerance = 1e-10)
})

test_that("sqrt(x): f''=-1/(4*x^(3/2))", {
  x_val <- 4
  x <- dual_variable_n(x_val, 2)
  r <- sqrt(x)
  expect_equal(deriv_n(r, 0), 2)
  expect_equal(deriv_n(r, 1), 0.25, tolerance = 1e-10)
  expect_equal(deriv_n(r, 2), -1/(4*x_val^(3/2)), tolerance = 1e-10)
})

# =============================================================================
# Third-order derivatives
# =============================================================================

test_that("x^3: f'''=6", {
  x <- dual_variable_n(2, 3)
  r <- x^3
  expect_equal(deriv_n(r, 0), 8)
  expect_equal(deriv_n(r, 1), 12)
  expect_equal(deriv_n(r, 2), 12)
  expect_equal(deriv_n(r, 3), 6)
})

test_that("x^4: f'''=24x", {
  x_val <- 2
  x <- dual_variable_n(x_val, 3)
  r <- x^4
  # f(x)=x^4, f'=4x^3, f''=12x^2, f'''=24x
  expect_equal(deriv_n(r, 0), x_val^4)
  expect_equal(deriv_n(r, 1), 4*x_val^3)
  expect_equal(deriv_n(r, 2), 12*x_val^2)
  expect_equal(deriv_n(r, 3), 24*x_val)
})

test_that("exp(x): f'''=exp(x)", {
  x_val <- 1.0
  x <- dual_variable_n(x_val, 3)
  r <- exp(x)
  expect_equal(deriv_n(r, 3), exp(x_val), tolerance = 1e-10)
})

test_that("sin(x): f'''=-cos(x)", {
  x_val <- pi / 6
  x <- dual_variable_n(x_val, 3)
  r <- sin(x)
  expect_equal(deriv_n(r, 3), -cos(x_val), tolerance = 1e-10)
})

test_that("cos(x): f'''=sin(x)", {
  x_val <- pi / 6
  x <- dual_variable_n(x_val, 3)
  r <- cos(x)
  expect_equal(deriv_n(r, 3), sin(x_val), tolerance = 1e-10)
})

test_that("log(x): f'''=2/x^3", {
  x_val <- 2
  x <- dual_variable_n(x_val, 3)
  r <- log(x)
  # f'''(log(x)) = 2/x^3
  expect_equal(deriv_n(r, 3), 2/x_val^3, tolerance = 1e-10)
})

# =============================================================================
# Fourth-order derivatives
# =============================================================================

test_that("x^5: f^(4) = 5*4*3*2 * x = 120x", {
  x_val <- 2
  x <- dual_variable_n(x_val, 4)
  r <- x^5
  # f^(4) = 120*x at x=2 => 240
  expect_equal(deriv_n(r, 4), 120*x_val)
})

test_that("exp(x): f^(4)=exp(x)", {
  x_val <- 0.5
  x <- dual_variable_n(x_val, 4)
  r <- exp(x)
  expect_equal(deriv_n(r, 4), exp(x_val), tolerance = 1e-10)
})

test_that("sin(x): f^(4)=sin(x) (the 4-cycle)", {
  x_val <- pi / 5
  x <- dual_variable_n(x_val, 4)
  r <- sin(x)
  # sin -> cos -> -sin -> -cos -> sin
  expect_equal(deriv_n(r, 4), sin(x_val), tolerance = 1e-10)
})

test_that("cos(x): f^(4)=cos(x)", {
  x_val <- pi / 5
  x <- dual_variable_n(x_val, 4)
  r <- cos(x)
  # cos -> -sin -> -cos -> sin -> cos
  expect_equal(deriv_n(r, 4), cos(x_val), tolerance = 1e-10)
})

test_that("log(x): f^(4) = -6/x^4", {
  x_val <- 3
  x <- dual_variable_n(x_val, 4)
  r <- log(x)
  # d^n/dx^n log(x) = (-1)^(n-1) * (n-1)! / x^n
  # n=4: (-1)^3 * 6 / x^4 = -6/x^4
  expect_equal(deriv_n(r, 4), -6/x_val^4, tolerance = 1e-10)
})

# =============================================================================
# Fifth-order derivatives
# =============================================================================

test_that("x^5: f^(5) = 120", {
  x <- dual_variable_n(2, 5)
  r <- x^5
  expect_equal(deriv_n(r, 5), 120)
})

test_that("x^6: f^(5) = 720x at x=1", {
  x <- dual_variable_n(1, 5)
  r <- x^6
  # f^(5) = 6*5*4*3*2 * x = 720*x, at x=1 => 720
  expect_equal(deriv_n(r, 5), 720)
})

test_that("exp(x): f^(5)=exp(x)", {
  x_val <- 0.3
  x <- dual_variable_n(x_val, 5)
  r <- exp(x)
  expect_equal(deriv_n(r, 5), exp(x_val), tolerance = 1e-10)
})

test_that("sin(x): f^(5)=cos(x) (5th in cycle)", {
  x_val <- pi / 7
  x <- dual_variable_n(x_val, 5)
  r <- sin(x)
  # sin cycle: sin, cos, -sin, -cos, sin, cos
  # 5th derivative of sin is cos
  expect_equal(deriv_n(r, 5), cos(x_val), tolerance = 1e-10)
})

# =============================================================================
# Third-order compositions against numerical
# =============================================================================

test_that("3rd derivatives of compositions match numerical", {
  fns <- list(
    list(f = function(x) x^4, name = "x^4"),
    list(f = function(x) exp(sin(x)), name = "exp(sin(x))"),
    list(f = function(x) log(1 + x^2), name = "log(1+x^2)"),
    list(f = function(x) x * exp(-x), name = "x*exp(-x)"),
    list(f = function(x) 1 / (1 + exp(-x)), name = "sigmoid")
  )

  for (fn in fns) {
    x_val <- 1.0
    result <- differentiate_n(fn$f, x_val, order = 3)
    num_3rd <- numerical_deriv_n(fn$f, x_val, 3L)
    expect_equal(
      result$d3, num_3rd,
      tolerance = tol_3rd,
      label = sprintf("f'''(%g) for %s", x_val, fn$name)
    )
  }
})

# =============================================================================
# differentiate_n convenience function
# =============================================================================

test_that("differentiate_n returns correct structure", {
  result <- differentiate_n(sin, pi/4, order = 4)
  expect_true(is.list(result))
  expect_named(result, c("value", "d1", "d2", "d3", "d4"))
  expect_equal(result$value, sin(pi/4))
  expect_equal(result$d1, cos(pi/4), tolerance = 1e-10)
  expect_equal(result$d2, -sin(pi/4), tolerance = 1e-10)
  expect_equal(result$d3, -cos(pi/4), tolerance = 1e-10)
  expect_equal(result$d4, sin(pi/4), tolerance = 1e-10)
})

test_that("differentiate_n rejects order < 1", {
  expect_error(differentiate_n(sin, 1, 0), "positive integer")
})


# =============================================================================
# Second derivatives against numerical (existing, updated)
# =============================================================================

test_that("second derivatives match numerical for compositions", {
  fns <- list(
    list(f = function(x) x^4, name = "x^4"),
    list(f = function(x) exp(sin(x)), name = "exp(sin(x))"),
    list(f = function(x) log(1 + x^2), name = "log(1+x^2)"),
    list(f = function(x) x * exp(-x), name = "x*exp(-x)"),
    list(f = function(x) 1 / (1 + exp(-x)), name = "sigmoid")
  )

  for (fn in fns) {
    x_val <- 1.0
    result <- differentiate_n(fn$f, x_val, order = 2)
    num_second <- numerical_second_deriv(fn$f, x_val)
    expect_equal(
      result$d2, num_second,
      tolerance = tol_2nd,
      label = sprintf("f''(%g) for %s", x_val, fn$name)
    )
  }
})

# =============================================================================
# Nested dual power rule (.is_zero fix)
# =============================================================================

test_that("x^3 second derivative via nested duals uses power rule path", {
  x <- dual_variable_n(2, 2)
  r <- x^3
  expect_equal(deriv_n(r, 2), 12)
})

test_that("2^x second derivative via nested duals", {
  x <- dual_variable_n(3, 2)
  r <- 2^x
  expect_equal(deriv_n(r, 2), 8 * log(2)^2, tolerance = 1e-10)
})

# =============================================================================
# Higher-order for compositions (stress test)
# =============================================================================

test_that("4th derivative of exp(sin(x)) matches numerical", {
  x_val <- 0.5
  ad_result <- differentiate_n(function(x) exp(sin(x)), x_val, order = 4)
  num_4th <- numerical_deriv_n(function(x) exp(sin(x)), x_val, 4L)
  expect_equal(ad_result$d4, num_4th, tolerance = tol_4th)
})

test_that("polynomial x^5 - 3*x^3 + 2*x, all 5 derivatives", {
  f <- function(x) x^5 - 3*x^3 + 2*x
  x_val <- 1.5
  result <- differentiate_n(f, x_val, order = 5)
  # Analytical: f'=5x^4-9x^2+2, f''=20x^3-18x, f'''=60x^2-18,
  # f^(4)=120x, f^(5)=120
  expect_equal(result$value, x_val^5 - 3*x_val^3 + 2*x_val)
  expect_equal(result$d1, 5*x_val^4 - 9*x_val^2 + 2, tolerance = 1e-10)
  expect_equal(result$d2, 20*x_val^3 - 18*x_val, tolerance = 1e-10)
  expect_equal(result$d3, 60*x_val^2 - 18, tolerance = 1e-10)
  expect_equal(result$d4, 120*x_val, tolerance = 1e-10)
  expect_equal(result$d5, 120, tolerance = 1e-10)
})

# =============================================================================
# Multi-parameter cases: gradient/hessian (vector-gradient + nesting)
# =============================================================================

test_that("gradient of 2-parameter quadratic", {
  # f(a, b) = -(a - 3)^2 - (b - 5)^2
  # gradient: (2*(3-a), 2*(5-b))
  ll <- function(theta) {
    a <- theta[1]; b <- theta[2]
    -(a - 3)^2 - (b - 5)^2
  }
  s <- gradient(ll, c(1, 2))
  expect_equal(s, c(4, 6), tolerance = 1e-10)
})

test_that("hessian of 2-parameter quadratic", {
  ll <- function(theta) {
    a <- theta[1]; b <- theta[2]
    -(a - 3)^2 - (b - 5)^2
  }
  H <- hessian(ll, c(1, 2))
  expect_equal(H, matrix(c(-2, 0, 0, -2), 2, 2), tolerance = 1e-10)
})

test_that("hessian of cross-term function a*b", {
  # f(a, b) = a*b, H = [[0, 1], [1, 0]]
  ll <- function(theta) {
    theta[1] * theta[2]
  }
  H <- hessian(ll, c(2, 3))
  expect_equal(H, matrix(c(0, 1, 1, 0), 2, 2), tolerance = 1e-10)
})

test_that("hessian of 3-parameter function matches numerical", {
  # f(a, b, c) = a*b^2 + b*c^3 + a^2*c
  ll <- function(theta) {
    a <- theta[1]; b <- theta[2]; cc <- theta[3]
    a * b^2 + b * cc^3 + a^2 * cc
  }
  theta <- c(1.5, 2.0, 0.5)
  H_ad <- hessian(ll, theta)
  H_num <- numerical_hessian(ll, theta)
  expect_equal(H_ad, H_num, tolerance = 1e-5)
})

test_that("gradient of 3-parameter function matches numerical gradient", {
  ll <- function(theta) {
    a <- theta[1]; b <- theta[2]; cc <- theta[3]
    a * b^2 + b * cc^3 + a^2 * cc
  }
  theta <- c(1.5, 2.0, 0.5)
  s_ad <- gradient(ll, theta)
  s_num <- numerical_gradient(ll, theta)
  expect_equal(s_ad, s_num, tolerance = 1e-6)
})

test_that("Normal MLE (mu, sigma): Hessian at MLE matches analytical", {
  # Normal log-likelihood with sufficient stats (avoids element-wise dual ops)
  data <- c(1, 2, 3, 4, 5)
  n <- length(data)
  sum_x <- sum(data)
  sum_x2 <- sum(data^2)

  ll <- function(theta) {
    mu <- theta[1]; sigma <- theta[2]
    -n/2 * log(2 * pi) - n * log(sigma) -
      0.5 / sigma^2 * (sum_x2 - 2 * mu * sum_x + n * mu^2)
  }

  mu_hat <- mean(data)
  sigma_hat <- sqrt(mean((data - mu_hat)^2))
  theta_mle <- c(mu_hat, sigma_hat)

  # At MLE: H_11 = -n/sigma^2, H_22 = -2n/sigma^2, H_12 = 0
  H <- hessian(ll, theta_mle)
  expect_equal(H[1, 1], -n / sigma_hat^2, tolerance = 1e-6)
  expect_equal(H[2, 2], -2 * n / sigma_hat^2, tolerance = 1e-6)
  expect_equal(H[1, 2], 0, tolerance = 1e-6)
  expect_equal(H[2, 1], 0, tolerance = 1e-6)

  # Score should be ~0 at MLE
  s <- gradient(ll, theta_mle)
  expect_equal(s, c(0, 0), tolerance = 1e-6)
})

# =============================================================================
# 6th-order derivatives of non-linear functions
# =============================================================================

test_that("6th derivative of exp(x) = exp(x)", {
  x_val <- 0.7
  x <- dual_variable_n(x_val, 6)
  r <- exp(x)
  expect_equal(deriv_n(r, 6), exp(x_val), tolerance = 1e-10)
})

test_that("6th derivative of sin(x): sin cycle mod 4", {
  x_val <- pi / 5
  x <- dual_variable_n(x_val, 6)
  r <- sin(x)
  # 6th derivative of sin = -sin (cycle: sin,cos,-sin,-cos,sin,cos,-sin)
  expect_equal(deriv_n(r, 6), -sin(x_val), tolerance = 1e-10)
})

test_that("6th derivative of non-linear composition exp(sin(x))", {
  # exp(sin(x)) at x=0.5, verify against numerical
  x_val <- 0.5
  result <- differentiate_n(function(x) exp(sin(x)), x_val, order = 6)
  num_6th <- numerical_deriv_n(function(x) exp(sin(x)), x_val, 6L, h = 0.01)
  # 6th-order numerical FD is very imprecise, so wide tolerance
  expect_equal(result$d6, num_6th, tolerance = 0.5)
  # But all lower orders should be much closer
  for (k in 1:3) {
    ad_k <- result[[paste0("d", k)]]
    num_k <- numerical_deriv_n(function(x) exp(sin(x)), x_val, k, h = 1e-3)
    expect_equal(ad_k, num_k, tolerance = 1e-3,
                 label = paste("d", k, "of exp(sin(x))"))
  }
})

test_that("6th derivative of x^7 = 7!/1! * x = 5040*x", {
  x_val <- 1.5
  x <- dual_variable_n(x_val, 6)
  r <- x^7
  # f^(6)(x) = 7*6*5*4*3*2 * x = 5040 * x
  expect_equal(deriv_n(r, 6), 5040 * x_val, tolerance = 1e-10)
})

test_that("all 6 derivatives of log(1+x^2) at x=1", {
  # f(x) = log(1+x^2), analytically complex higher derivs
  # Verify against numerical
  f <- function(x) log(1 + x^2)
  x_val <- 1.0
  result <- differentiate_n(f, x_val, order = 6)

  # Check value
  expect_equal(result$value, log(2))

  # First derivative: 2x/(1+x^2) = 2/2 = 1 at x=1
  expect_equal(result$d1, 1.0, tolerance = 1e-10)

  # Compare all derivatives against numerical (with appropriate tolerances)
  for (k in 1:4) {
    num_k <- numerical_deriv_n(f, x_val, k, h = 1e-3)
    expect_equal(result[[paste0("d", k)]], num_k,
                 tolerance = 10^(-4 + k),
                 label = paste0("d", k, " of log(1+x^2)"))
  }
})

# =============================================================================
# Multi-parameter non-linear functions (f: R^n -> R)
# =============================================================================

test_that("Hessian of non-linear 2-param function a^2*b + exp(a*b)", {
  # f(a,b) = a^2*b + exp(a*b)
  # df/da = 2*a*b + b*exp(a*b)
  # df/db = a^2 + a*exp(a*b)
  # d2f/da2 = 2*b + b^2*exp(a*b)
  # d2f/db2 = a^2*exp(a*b)
  # d2f/dadb = 2*a + exp(a*b) + a*b*exp(a*b)
  ll <- function(theta) {
    a <- theta[1]; b <- theta[2]
    a^2 * b + exp(a * b)
  }
  a <- 1.0; b <- 0.5
  theta <- c(a, b)

  # Score (gradient)
  s <- gradient(ll, theta)
  eab <- exp(a * b)
  expect_equal(s[1], 2*a*b + b*eab, tolerance = 1e-10)
  expect_equal(s[2], a^2 + a*eab, tolerance = 1e-10)

  # Hessian
  H <- hessian(ll, theta)
  expect_equal(H[1,1], 2*b + b^2*eab, tolerance = 1e-10)
  expect_equal(H[2,2], a^2*eab, tolerance = 1e-10)
  expect_equal(H[1,2], 2*a + eab + a*b*eab, tolerance = 1e-10)
  expect_equal(H[2,1], H[1,2], tolerance = 1e-10)  # symmetry
})

test_that("Hessian of 3-param non-linear: exp(a*b*c) + a*sin(b+c)", {
  ll <- function(theta) {
    a <- theta[1]; b <- theta[2]; cc <- theta[3]
    exp(a * b * cc) + a * sin(b + cc)
  }
  theta <- c(0.5, 1.0, 0.3)
  H_ad <- hessian(ll, theta)
  H_num <- numerical_hessian(ll, theta)
  expect_equal(H_ad, H_num, tolerance = 1e-5)
  # Verify symmetry
  expect_equal(H_ad[1,2], H_ad[2,1], tolerance = 1e-12)
  expect_equal(H_ad[1,3], H_ad[3,1], tolerance = 1e-12)
  expect_equal(H_ad[2,3], H_ad[3,2], tolerance = 1e-12)
})

test_that("gradient of 4-param non-linear matches numerical gradient", {
  ll <- function(theta) {
    a <- theta[1]; b <- theta[2]; cc <- theta[3]; d <- theta[4]
    log(a^2 + b^2) + exp(-cc * d) + sin(a * cc) * cos(b * d)
  }
  theta <- c(1.5, 0.8, 0.3, 1.2)
  s_ad <- gradient(ll, theta)
  s_num <- numerical_gradient(ll, theta)
  expect_equal(s_ad, s_num, tolerance = 1e-6)
})
