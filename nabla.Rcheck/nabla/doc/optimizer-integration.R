## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ----setup--------------------------------------------------------------------
library(nabla)

## -----------------------------------------------------------------------------
set.seed(42)
data_norm <- rnorm(100, mean = 5, sd = 2)
n <- length(data_norm)
sum_x <- sum(data_norm)
sum_x2 <- sum(data_norm^2)

ll_normal <- function(x) {
  mu <- x[1]
  sigma <- exp(x[2])  # unconstrained: x[2] = log(sigma)
  -n * log(sigma) - (1 / (2 * sigma^2)) * (sum_x2 - 2 * mu * sum_x + n * mu^2)
}

## -----------------------------------------------------------------------------
# Negated versions for minimization
nll <- function(x) -ll_normal(x)
ngr <- function(x) -gradient(ll_normal, x)

# Starting value: mu = 0, log(sigma) = 0 (i.e., sigma = 1)
start <- c(0, 0)

# Without gradient: optim uses internal finite differences
fit_no_gr <- optim(start, fn = nll, method = "BFGS")

# With AD gradient
fit_ad_gr <- optim(start, fn = nll, gr = ngr, method = "BFGS")

# Compare (transform x[2] back to sigma for display)
data.frame(
  method = c("BFGS (no gradient)", "BFGS (AD gradient)"),
  mu = c(fit_no_gr$par[1], fit_ad_gr$par[1]),
  sigma = exp(c(fit_no_gr$par[2], fit_ad_gr$par[2])),
  fn_calls = c(fit_no_gr$counts["function"], fit_ad_gr$counts["function"]),
  gr_calls = c(fit_no_gr$counts["gradient"], fit_ad_gr$counts["gradient"]),
  convergence = c(fit_no_gr$convergence, fit_ad_gr$convergence)
)

## -----------------------------------------------------------------------------
# Verify against analytical MLE
mle_mu <- mean(data_norm)
mle_sigma <- sqrt(mean((data_norm - mle_mu)^2))
cat("Analytical MLE: mu =", round(mle_mu, 6), " sigma =", round(mle_sigma, 6), "\n")
cat("optim+AD MLE:   mu =", round(fit_ad_gr$par[1], 6),
    " sigma =", round(exp(fit_ad_gr$par[2]), 6), "\n")

## -----------------------------------------------------------------------------
nhess <- function(x) -hessian(ll_normal, x)

fit_nlminb <- nlminb(start,
                     objective = nll,
                     gradient = ngr,
                     hessian = nhess)

cat("nlminb MLE: mu =", round(fit_nlminb$par[1], 6),
    " sigma =", round(exp(fit_nlminb$par[2]), 6), "\n")
cat("Converged:", fit_nlminb$convergence == 0, "\n")
cat("Iterations:", fit_nlminb$iterations, "\n")
cat("Function evaluations:", fit_nlminb$evaluations["function"], "\n")
cat("Gradient evaluations:", fit_nlminb$evaluations["gradient"], "\n")

## -----------------------------------------------------------------------------
# Simulate data
set.seed(7)
n_lr <- 200
x1 <- rnorm(n_lr)
x2 <- rnorm(n_lr)
X <- cbind(1, x1, x2)
beta_true <- c(-0.5, 1.2, -0.8)
eta_true <- X %*% beta_true
prob_true <- 1 / (1 + exp(-eta_true))
y <- rbinom(n_lr, 1, prob_true)

# Log-likelihood for nabla
ll_logistic <- function(x) {
  result <- 0
  for (i in seq_len(n_lr)) {
    eta_i <- x[1] * X[i, 1] + x[2] * X[i, 2] + x[3] * X[i, 3]
    result <- result + y[i] * eta_i - log(1 + exp(eta_i))
  }
  result
}

# Numeric version for optim's fn
ll_logistic_num <- function(beta) {
  eta <- X %*% beta
  sum(y * eta - log(1 + exp(eta)))
}

## -----------------------------------------------------------------------------
nll_lr <- function(x) -ll_logistic_num(x)
ngr_lr <- function(x) -gradient(ll_logistic, x)

fit_lr <- optim(c(0, 0, 0), fn = nll_lr, gr = ngr_lr, method = "BFGS")

## -----------------------------------------------------------------------------
fit_glm <- glm(y ~ x1 + x2, family = binomial)

# Coefficient comparison
data.frame(
  parameter = c("Intercept", "x1", "x2"),
  optim_AD = round(fit_lr$par, 6),
  glm = round(coef(fit_glm), 6),
  difference = round(fit_lr$par - coef(fit_glm), 10)
)

## -----------------------------------------------------------------------------
# Observed information at the optimum (negative Hessian)
obs_info <- -hessian(ll_logistic, fit_lr$par)

# Variance-covariance matrix
vcov_ad <- solve(obs_info)

# Standard errors
se_ad <- sqrt(diag(vcov_ad))

# Compare with glm's standard errors
se_glm <- summary(fit_glm)$coefficients[, "Std. Error"]

data.frame(
  parameter = c("Intercept", "x1", "x2"),
  SE_AD = round(se_ad, 6),
  SE_glm = round(se_glm, 6),
  ratio = round(se_ad / se_glm, 8)
)

## ----fig-convergence, fig.width=6, fig.height=4.5-----------------------------
# Track convergence for three methods on the Normal(mu, sigma) model
# We'll use a callback-style approach via optim's trace

# Method 1: Nelder-Mead (no gradient)
fit_nm <- optim(start, fn = nll, method = "Nelder-Mead",
                control = list(maxit = 200))

# Method 2: BFGS with AD gradient
fit_bfgs <- optim(start, fn = nll, gr = ngr, method = "BFGS")

# Method 3: nlminb with AD gradient + Hessian
fit_nlm <- nlminb(start, objective = nll, gradient = ngr, hessian = nhess)

# Summary comparison
data.frame(
  method = c("Nelder-Mead", "BFGS + AD grad", "nlminb + AD grad+hess"),
  fn_evals = c(fit_nm$counts["function"],
               fit_bfgs$counts["function"],
               fit_nlm$evaluations["function"]),
  converged = c(fit_nm$convergence == 0,
                fit_bfgs$convergence == 0,
                fit_nlm$convergence == 0),
  nll = c(fit_nm$value, fit_bfgs$value, fit_nlm$objective)
)

## ----fig-optim-contour, fig.width=6, fig.height=5-----------------------------
# Contour plot with optimizer paths
# Recompute paths by running each optimizer and collecting iterates

# Helper: collect iterates using an environment (avoids <<-)
make_tracer <- function(fn) {
  env <- new.env(parent = emptyenv())
  env$trace <- list()
  list(
    fn = function(x) {
      env$trace[[length(env$trace) + 1L]] <- x
      fn(x)
    },
    path = function() do.call(rbind, env$trace)
  )
}

# Collect iterates via BFGS
bfgs <- make_tracer(function(x) -ll_normal(x))
optim(start, fn = bfgs$fn, gr = ngr, method = "BFGS")
bfgs_path <- bfgs$path()

# Collect iterates via Nelder-Mead
nm <- make_tracer(function(x) -ll_normal(x))
optim(start, fn = nm$fn, method = "Nelder-Mead",
      control = list(maxit = 200))
nm_path <- nm$path()

# Collect iterates via nlminb
nlm <- make_tracer(function(x) -ll_normal(x))
nlminb(start, objective = nlm$fn, gradient = ngr, hessian = nhess)
nlm_path <- nlm$path()

# Build contour grid
mu_grid <- seq(-1, 7, length.out = 100)
sigma_grid <- seq(0.5, 4, length.out = 100)
nll_surface <- outer(mu_grid, sigma_grid, Vectorize(function(m, s) {
  -ll_normal(c(m, log(s)))
}))

# Transform paths from (mu, log_sigma) to (mu, sigma) for plotting
bfgs_path[, 2] <- exp(bfgs_path[, 2])
nm_path[, 2] <- exp(nm_path[, 2])
nlm_path[, 2] <- exp(nlm_path[, 2])

oldpar <- par(mar = c(4, 4, 2, 1))
contour(mu_grid, sigma_grid, nll_surface, nlevels = 30,
        xlab = expression(mu), ylab = expression(sigma),
        main = "Optimizer paths on NLL contours",
        col = "grey75")

# Nelder-Mead path (subsample for clarity)
nm_sub <- nm_path[seq(1, nrow(nm_path), length.out = min(50, nrow(nm_path))), ]
lines(nm_sub[, 1], nm_sub[, 2], col = "coral", lwd = 1.5, lty = 3)
points(nm_sub[1, 1], nm_sub[1, 2], pch = 17, col = "coral", cex = 1.2)

# BFGS path
lines(bfgs_path[, 1], bfgs_path[, 2], col = "steelblue", lwd = 2, type = "o",
      pch = 19, cex = 0.6)

# nlminb path
lines(nlm_path[, 1], nlm_path[, 2], col = "forestgreen", lwd = 2, type = "o",
      pch = 15, cex = 0.6)

# MLE
points(mle_mu, mle_sigma, pch = 3, col = "black", cex = 2, lwd = 2)

legend("topright",
       legend = c("Nelder-Mead", "BFGS + AD grad", "nlminb + AD grad+hess", "MLE"),
       col = c("coral", "steelblue", "forestgreen", "black"),
       lty = c(3, 1, 1, NA), pch = c(17, 19, 15, 3),
       lwd = c(1.5, 2, 2, 2), bty = "n", cex = 0.85)
par(oldpar)

