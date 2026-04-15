## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ----setup--------------------------------------------------------------------
library(nabla)

## -----------------------------------------------------------------------------
f <- function(x) x[1]^2 + sin(x[1])
D(f, pi / 4)   # exact first derivative at pi/4

## -----------------------------------------------------------------------------
g <- function(x) x[1]^2 * exp(x[2])
gradient(g, c(2, 0))   # c(4, 4)
hessian(g, c(2, 0))    # 2x2 Hessian matrix

## -----------------------------------------------------------------------------
f <- function(x) x[1]^2 + sin(x[1])

# First derivative: f'(x) = 2x + cos(x)
D(f, pi / 4)

# Verify against the analytical formula
2 * (pi / 4) + cos(pi / 4)

# Second derivative: f''(x) = 2 - sin(x)
D(f, pi / 4, order = 2)

## -----------------------------------------------------------------------------
gradient(f, pi / 4)

## ----fig-function-derivative, fig.width=6, fig.height=4-----------------------
f <- function(x) x[1]^2 + sin(x[1])
xs <- seq(0, 2 * pi, length.out = 200)

# Compute f(x) and f'(x) at each grid point using D()
fx <- sapply(xs, function(xi) f(xi))
fpx <- sapply(xs, function(xi) D(f, xi))

# Mark the evaluation point x = pi/4
x_mark <- pi / 4

oldpar <- par(mar = c(4, 4, 2, 1))
plot(xs, fx, type = "l", col = "steelblue", lwd = 2,
     xlab = "x", ylab = "y",
     main = expression(f(x) == x^2 + sin(x) ~ "and its derivative"),
     ylim = range(c(fx, fpx)))
lines(xs, fpx, col = "firebrick", lwd = 2, lty = 2)
points(x_mark, f(x_mark), pch = 19, col = "steelblue", cex = 1.3)
points(x_mark, D(f, x_mark), pch = 19, col = "firebrick", cex = 1.3)
abline(v = x_mark, col = "grey60", lty = 3)
legend("topleft", legend = c("f(x)", "f'(x) via AD"),
       col = c("steelblue", "firebrick"), lty = c(1, 2), lwd = 2,
       bty = "n")
par(oldpar)

## -----------------------------------------------------------------------------
# if/else branching
safe_log <- function(x) {
  if (x[1] > 0) log(x[1]) else -Inf
}
D(safe_log, 2)    # 1/2 = 0.5
D(safe_log, -1)   # 0 (constant branch)

# for loop
poly <- function(x) {
  result <- 0
  for (i in 1:5) result <- result + x[1]^i
  result
}
D(poly, 2)   # 1 + 2*2 + 3*4 + 4*8 + 5*16 = 129

# Reduce
f_reduce <- function(x) Reduce("+", lapply(1:4, function(i) x[1]^i))
D(f_reduce, 2)

## -----------------------------------------------------------------------------
g <- function(x) exp(-x[1]^2 / 2) / sqrt(2 * pi)  # standard normal PDF

# D(g, 1) should equal -x * dnorm(x) at x = 1
D(g, 1)
-1 * dnorm(1)
D(g, 1) - (-1 * dnorm(1))  # ~0

## -----------------------------------------------------------------------------
# Gradient of a 2-parameter function
f2 <- function(x) x[1]^2 * x[2]
gradient(f2, c(3, 4))  # c(2*3*4, 3^2) = c(24, 9)

# Hessian
hessian(f2, c(3, 4))

# Jacobian of a vector-valued function
f_vec <- function(x) list(x[1] * x[2], x[1]^2 + x[2])
jacobian(f_vec, c(3, 4))

## -----------------------------------------------------------------------------
f2 <- function(x) x[1]^2 * x[2]
D(f2, c(3, 4), order = 2)           # Hessian
D(f2, c(3, 4), order = 3)           # 2x2x2 third-order tensor

## -----------------------------------------------------------------------------
f <- function(x) x[1]^3 * sin(x[1])
f_num <- function(x) x^3 * sin(x)  # plain numeric version

x0 <- 2

# 1. Analytical: f'(x) = 3x^2 sin(x) + x^3 cos(x)
analytical <- 3 * x0^2 * sin(x0) + x0^3 * cos(x0)

# 2. Finite differences
h <- 1e-8
finite_diff <- (f_num(x0 + h) - f_num(x0 - h)) / (2 * h)

# 3. Automatic differentiation
ad_result <- D(f, x0)

# Compare
data.frame(
  method = c("Analytical", "Finite Diff", "AD (nabla)"),
  derivative = c(analytical, finite_diff, ad_result),
  error_vs_analytical = c(0, finite_diff - analytical, ad_result - analytical)
)

## ----fig-error-comparison, fig.width=5, fig.height=3.5------------------------
errors <- abs(c(finite_diff - analytical, ad_result - analytical))
# Use log10 scale; clamp AD error to .Machine$double.eps if exactly zero
errors[errors == 0] <- .Machine$double.eps

oldpar <- par(mar = c(4, 5, 2, 1))
bp <- barplot(log10(errors),
              names.arg = c("Finite Diff", "AD (nabla)"),
              col = c("coral", "steelblue"),
              ylab = expression(log[10] ~ "|error|"),
              main = "Absolute error vs analytical derivative",
              ylim = c(-16, 0), border = NA)
abline(h = log10(.Machine$double.eps), lty = 2, col = "grey40")
text(mean(bp), log10(.Machine$double.eps) + 0.8,
     "machine epsilon", cex = 0.8, col = "grey40")
par(oldpar)

## -----------------------------------------------------------------------------
# A dual variable: value = 3, derivative seed = 1
x <- dual_variable(3)
value(x)
deriv(x)

# A dual constant: value = 5, derivative seed = 0
k <- dual_constant(5)
value(k)
deriv(k)

# Explicit constructor
y <- dual(2, 1)
value(y)
deriv(y)

## -----------------------------------------------------------------------------
x <- dual_variable(3)

# Addition: d/dx(x + 2) = 1
r_add <- x + 2
value(r_add)
deriv(r_add)

# Subtraction: d/dx(5 - x) = -1
r_sub <- 5 - x
value(r_sub)
deriv(r_sub)

# Multiplication: d/dx(x * 4) = 4
r_mul <- x * 4
value(r_mul)
deriv(r_mul)

# Division: d/dx(1/x) = -1/x^2 = -1/9
r_div <- 1 / x
value(r_div)
deriv(r_div)

# Power: d/dx(x^3) = 3*x^2 = 27
r_pow <- x^3
value(r_pow)
deriv(r_pow)

## -----------------------------------------------------------------------------
x <- dual_variable(1)

# exp: d/dx exp(x) = exp(x)
r_exp <- exp(x)
value(r_exp)
deriv(r_exp)

# log: d/dx log(x) = 1/x
r_log <- log(x)
value(r_log)
deriv(r_log)

# sin: d/dx sin(x) = cos(x)
x2 <- dual_variable(pi / 4)
r_sin <- sin(x2)
value(r_sin)
deriv(r_sin)  # cos(pi/4)

# sqrt: d/dx sqrt(x) = 1/(2*sqrt(x))
x3 <- dual_variable(4)
r_sqrt <- sqrt(x3)
value(r_sqrt)
deriv(r_sqrt)  # 1/(2*2) = 0.25

# Gamma-related
x4 <- dual_variable(3)
r_lgamma <- lgamma(x4)
value(r_lgamma)     # log(2!) = log(2)
deriv(r_lgamma)     # digamma(3)

## -----------------------------------------------------------------------------
# sum() works
a <- dual_variable(2)
b <- dual_constant(3)
total <- sum(a, b, dual_constant(1))
value(total)  # 6
deriv(total)  # 1 (only a has deriv = 1)

# prod() works
p <- prod(a, dual_constant(3))
value(p)  # 6
deriv(p)  # 3 (product rule: 3*1 + 2*0)

# c() creates a dual_vector
v <- c(a, b)
length(v)

# is.numeric returns TRUE (for compatibility)
is.numeric(dual_variable(1))

