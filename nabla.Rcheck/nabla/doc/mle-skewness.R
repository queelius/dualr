## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ----setup--------------------------------------------------------------------
library(nabla)

## -----------------------------------------------------------------------------
# Simulate data
set.seed(42)
n <- 200
alpha_true <- 3
beta_true <- 2
x_data <- rgamma(n, shape = alpha_true, rate = beta_true)

# Sufficient statistics
S1 <- sum(log(x_data))
S2 <- sum(x_data)

# Log-likelihood using sufficient statistics
ll_gamma <- function(theta) {
  alpha <- theta[1]; beta <- theta[2]
  n * alpha * log(beta) - n * lgamma(alpha) +
    (alpha - 1) * S1 - beta * S2
}

## -----------------------------------------------------------------------------
fit <- optim(
  par = c(2, 1),
  fn = function(theta) -ll_gamma(theta),
  gr = function(theta) -gradient(ll_gamma, theta),
  method = "L-BFGS-B",
  lower = c(1e-4, 1e-4)
)
theta_hat <- fit$par
names(theta_hat) <- c("alpha", "beta")
theta_hat

## -----------------------------------------------------------------------------
# Score at MLE — should be near zero
score <- gradient(ll_gamma, theta_hat)
score

# Observed information
H <- hessian(ll_gamma, theta_hat)
J <- -H
cat("Observed information matrix:\n")
J

## -----------------------------------------------------------------------------
alpha_hat <- theta_hat[1]; beta_hat <- theta_hat[2]
J_analytical <- matrix(c(
  n * trigamma(alpha_hat),   -n / beta_hat,
  -n / beta_hat,              n * alpha_hat / beta_hat^2
), nrow = 2, byrow = TRUE)

cat("Max |AD - analytical| in J:", max(abs(J - J_analytical)), "\n")

## -----------------------------------------------------------------------------
J_inv <- solve(J)
se <- sqrt(diag(J_inv))
cat("SE(alpha_hat) =", se[1], "\n")
cat("SE(beta_hat)  =", se[2], "\n")

## -----------------------------------------------------------------------------
T3 <- D(ll_gamma, theta_hat, order = 3)
T3

## -----------------------------------------------------------------------------
T3_analytical <- array(0, dim = c(2, 2, 2))
T3_analytical[1, 1, 1] <- -n * psigamma(alpha_hat, deriv = 2)
# T3[1,1,2] and permutations = 0 (already zero)
T3_analytical[1, 2, 2] <- n / beta_hat^2
T3_analytical[2, 1, 2] <- n / beta_hat^2
T3_analytical[2, 2, 1] <- n / beta_hat^2
T3_analytical[2, 2, 2] <- -2 * n * alpha_hat / beta_hat^3

cat("Max |AD - analytical| in T3:", max(abs(T3 - T3_analytical)), "\n")

## -----------------------------------------------------------------------------
mle_skewness <- function(J_inv, T3, n) {
  p <- nrow(J_inv)
  skew <- numeric(p)
  for (a in seq_len(p)) {
    s <- 0
    for (r in seq_len(p)) {
      for (s_ in seq_len(p)) {
        for (t in seq_len(p)) {
          s <- s + J_inv[a, r] * J_inv[a, s_] * J_inv[a, t] *
            4 * T3[r, s_, t] / n
        }
      }
    }
    skew[a] <- s / J_inv[a, a]^(3/2)
  }
  skew
}

theoretical_skew <- mle_skewness(J_inv, T3, n)
cat("Theoretical skewness of alpha_hat:", theoretical_skew[1], "\n")
cat("Theoretical skewness of beta_hat: ", theoretical_skew[2], "\n")

## -----------------------------------------------------------------------------
skewness <- function(x) {
  m <- mean(x); s <- sd(x)
  mean(((x - m) / s)^3)
}

set.seed(42)
B <- 5000
alpha_hats <- numeric(B)
beta_hats <- numeric(B)

for (b in seq_len(B)) {
  x_b <- rgamma(n, shape = alpha_true, rate = beta_true)
  S1_b <- sum(log(x_b))
  S2_b <- sum(x_b)

  ll_b <- function(theta) {
    a <- theta[1]; be <- theta[2]
    n * a * log(be) - n * lgamma(a) + (a - 1) * S1_b - be * S2_b
  }

  fit_b <- optim(
    par = c(2, 1),
    fn = function(theta) -ll_b(theta),
    method = "L-BFGS-B",
    lower = c(1e-4, 1e-4)
  )
  alpha_hats[b] <- fit_b$par[1]
  beta_hats[b] <- fit_b$par[2]
}

empirical_skew <- c(skewness(alpha_hats), skewness(beta_hats))

## -----------------------------------------------------------------------------
comparison <- data.frame(
  parameter = c("alpha", "beta"),
  theoretical = round(theoretical_skew, 4),
  empirical = round(empirical_skew, 4)
)
comparison

## ----fig-mle-histograms, fig.width=7, fig.height=3.5--------------------------
oldpar <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

# alpha_hat
hist(alpha_hats, breaks = 50, freq = FALSE, col = "grey85", border = "grey60",
     main = expression(hat(alpha)),
     xlab = expression(hat(alpha)))
curve(dnorm(x, mean = mean(alpha_hats), sd = sd(alpha_hats)),
      col = "steelblue", lwd = 2, add = TRUE)
abline(v = alpha_true, col = "firebrick", lty = 2, lwd = 2)
legend("topright",
       legend = c("Normal approx", expression(alpha[true])),
       col = c("steelblue", "firebrick"), lty = c(1, 2), lwd = 2,
       bty = "n", cex = 0.8)

# beta_hat
hist(beta_hats, breaks = 50, freq = FALSE, col = "grey85", border = "grey60",
     main = expression(hat(beta)),
     xlab = expression(hat(beta)))
curve(dnorm(x, mean = mean(beta_hats), sd = sd(beta_hats)),
      col = "steelblue", lwd = 2, add = TRUE)
abline(v = beta_true, col = "firebrick", lty = 2, lwd = 2)
legend("topright",
       legend = c("Normal approx", expression(beta[true])),
       col = c("steelblue", "firebrick"), lty = c(1, 2), lwd = 2,
       bty = "n", cex = 0.8)
par(oldpar)

