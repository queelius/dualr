## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ----setup--------------------------------------------------------------------
library(nabla)

## -----------------------------------------------------------------------------
f <- function(x) x[1]^3

D(f, 2)                # f'(2) = 3*4 = 12
D(f, 2, order = 2)     # f''(2) = 6*2 = 12
D(f, 2, order = 3)     # f'''(2) = 6

## -----------------------------------------------------------------------------
# Gradient of a 2-parameter function
f <- function(x) x[1]^2 * x[2]
D(f, c(3, 4))  # gradient: c(24, 9)

# Hessian via D^2
D(f, c(3, 4), order = 2)

# Composition works identically
Df <- D(f)
DDf <- D(Df)
DDf(c(3, 4))

## -----------------------------------------------------------------------------
g <- function(x) list(x[1] * x[2], x[1]^2 + x[2])
D(g, c(3, 4))  # 2x2 Jacobian: [[4, 3], [6, 1]]

## -----------------------------------------------------------------------------
gradient(f, c(3, 4))           # == D(f, c(3, 4))
hessian(f, c(3, 4))            # == D(f, c(3, 4), order = 2)
jacobian(g, c(3, 4))           # == D(g, c(3, 4))

## -----------------------------------------------------------------------------
curvature <- function(f, x) {
  d1 <- D(f, x)
  d2 <- D(f, x, order = 2)
  abs(d2) / (1 + d1^2)^(3/2)
}

# Curvature of sin(x) at various points
# Wrap sin for D()'s x[1] convention
sin_f <- function(x) sin(x[1])
xs <- seq(0, 2 * pi, length.out = 7)
kappas <- sapply(xs, function(x) curvature(sin_f, x))
data.frame(x = round(xs, 3), curvature = round(kappas, 6))

## -----------------------------------------------------------------------------
curvature(sin_f, pi / 2)  # should be 1.0

## ----fig-sin-curvature, fig.width=6, fig.height=4-----------------------------
xs_curve <- seq(0, 2 * pi, length.out = 200)
sin_vals <- sin(xs_curve)
kappa_vals <- sapply(xs_curve, function(x) curvature(sin_f, x))

oldpar <- par(mar = c(4, 4.5, 2, 4.5))
plot(xs_curve, sin_vals, type = "l", col = "steelblue", lwd = 2,
     xlab = "x", ylab = "sin(x)",
     main = "sin(x) and its curvature")
par(new = TRUE)
plot(xs_curve, kappa_vals, type = "l", col = "firebrick", lwd = 2, lty = 2,
     axes = FALSE, xlab = "", ylab = "")
axis(4, col.axis = "firebrick")
mtext(expression(kappa(x)), side = 4, line = 2.5, col = "firebrick")
abline(v = c(pi / 2, 3 * pi / 2), col = "grey60", lty = 3)
legend("bottom",
       legend = c("sin(x)", expression(kappa(x))),
       col = c("steelblue", "firebrick"), lty = c(1, 2), lwd = 2,
       bty = "n", horiz = TRUE)
par(oldpar)

## -----------------------------------------------------------------------------
taylor2 <- function(f, x0, x) {
  f0 <- f(x0)
  f1 <- D(f, x0)
  f2 <- D(f, x0, order = 2)
  f0 + f1 * (x - x0) + 0.5 * f2 * (x - x0)^2
}

# Approximate exp(x) near x = 0
f_exp <- function(x) exp(x[1])
xs <- c(-0.1, -0.01, 0, 0.01, 0.1)
data.frame(
  x = xs,
  exact = exp(xs),
  taylor2 = sapply(xs, function(x) taylor2(f_exp, 0, x)),
  error = exp(xs) - sapply(xs, function(x) taylor2(f_exp, 0, x))
)

## ----fig-taylor-vs-exact, fig.width=6, fig.height=4---------------------------
xs_plot <- seq(-2, 3, length.out = 300)
exact_vals <- exp(xs_plot)
taylor_vals <- sapply(xs_plot, function(x) taylor2(f_exp, 0, x))

oldpar <- par(mar = c(4, 4, 2, 1))
plot(xs_plot, exact_vals, type = "l", col = "steelblue", lwd = 2,
     xlab = "x", ylab = "y",
     main = expression("exp(x) vs 2nd-order Taylor at " * x[0] == 0),
     ylim = c(-1, 12))
lines(xs_plot, taylor_vals, col = "firebrick", lwd = 2, lty = 2)
abline(v = 0, col = "grey60", lty = 3)
points(0, 1, pch = 19, col = "grey40", cex = 1.2)
legend("topleft",
       legend = c("exp(x)", expression("Taylor " * T[2](x))),
       col = c("steelblue", "firebrick"), lty = c(1, 2), lwd = 2,
       bty = "n")
par(oldpar)

## -----------------------------------------------------------------------------
set.seed(123)
data_pois <- rpois(50, lambda = 3)
n <- length(data_pois)
sum_x <- sum(data_pois)
sum_lfact <- sum(lfactorial(data_pois))

ll_poisson <- function(theta) {
  lambda <- theta[1]
  sum_x * log(lambda) - n * lambda - sum_lfact
}

lambda0 <- 2.5

## -----------------------------------------------------------------------------
hess_helper <- hessian(ll_poisson, lambda0)
hess_helper

## -----------------------------------------------------------------------------
se <- 1 / sqrt(-hess_helper[1, 1])
cat("SE(lambda) at lambda =", lambda0, ":", se, "\n")

## -----------------------------------------------------------------------------
# Build a dual_variable_n wrapped in a dual_vector
manual_theta <- dual_vector(list(dual_variable_n(lambda0, 2)))
result_manual <- ll_poisson(manual_theta)
manual_hess <- deriv(deriv(result_manual))
manual_hess

## -----------------------------------------------------------------------------
hess_helper[1, 1] - manual_hess  # ~0

## -----------------------------------------------------------------------------
x <- dual_variable_n(2, order = 2)

# Evaluate x^3
result <- x^3

# Extract all three quantities
deriv_n(result, 0)  # f(2) = 8
deriv_n(result, 1)  # f'(2) = 3*4 = 12
deriv_n(result, 2)  # f''(2) = 6*2 = 12

## -----------------------------------------------------------------------------
k <- dual_constant_n(5, order = 2)
deriv_n(k, 0)  # 5
deriv_n(k, 1)  # 0
deriv_n(k, 2)  # 0

## -----------------------------------------------------------------------------
# Differentiate sin(x) at x = pi/4
result <- differentiate_n(sin, pi / 4, order = 2)
result$value   # sin(pi/4)
result$d1      # cos(pi/4) = f'
result$d2      # -sin(pi/4) = f''

## -----------------------------------------------------------------------------
result$value - sin(pi / 4)     # ~0
result$d1 - cos(pi / 4)        # ~0
result$d2 - (-sin(pi / 4))     # ~0

## -----------------------------------------------------------------------------
f <- function(x) x * exp(-x^2)
d2 <- differentiate_n(f, 1, order = 2)

# Analytical: f'(x) = exp(-x^2)(1 - 2x^2)
# f''(x) = exp(-x^2)(-6x + 4x^3)
analytical_f1 <- exp(-1) * (1 - 2)
analytical_f2 <- exp(-1) * (-6 + 4)
d2$d1 - analytical_f1   # ~0
d2$d2 - analytical_f2   # ~0

