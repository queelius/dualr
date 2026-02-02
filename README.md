
# nabla <img src="man/figures/logo.png" align="right" height="139" alt="" />

> Exact Derivatives via Automatic Differentiation

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/nabla)](https://CRAN.R-project.org/package=nabla)
[![R-CMD-check](https://github.com/queelius/nabla/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/queelius/nabla/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Why automatic differentiation?

Computing derivatives is fundamental to optimization, statistics, and scientific
computing. The three main approaches are:

| | **Finite Differences** | **Symbolic Diff** | **AD (nabla)** |
|---|---|---|---|
| **Accuracy** | O(h) or O(h²) truncation error | Exact | Exact (machine precision) |
| **Cost** | 2p function evaluations | Expression swell for complex f | p forward passes |
| **Control flow** | Works | Breaks on if/for/while | Works through any code |
| **Implementation** | Easy | Requires CAS | Operator overloading |

Finite differences are inaccurate for ill-conditioned problems. Symbolic
differentiation suffers from expression swell and cannot handle control flow.
**Automatic differentiation** gives exact derivatives at machine precision
through any R code — loops, branches, and all.

## How it works

A dual number extends the reals with an infinitesimal ε where ε² = 0:

$$f(x + \varepsilon) = f(x) + f'(x)\,\varepsilon$$

By propagating ε through every arithmetic operation and math function,
`nabla` computes f(x) and f'(x) simultaneously in a single forward pass.
No symbolic manipulation, no finite step sizes — just exact derivatives
via the algebra of dual numbers.

## Installation

```r
# Install from CRAN
install.packages("nabla")

# Or install development version from GitHub
remotes::install_github("queelius/nabla")
```

## Quick start

**Single variable derivative:**

```r
library(nabla)

# Create a dual variable at x = 2 (seeds derivative tracking)
x <- dual_variable(2)

# Evaluate any expression — derivatives propagate automatically
result <- x^3 + sin(x)
value(result)   # f(2) = 8.909
deriv(result)   # f'(2) = 3*4 + cos(2) = 11.584
```

**Multi-parameter derivatives:**

```r
f <- function(x) x[1]^2 * x[2] + sin(x[2])

# Gradient — exact, no finite differences
gradient(f, c(3, 4))

# Hessian matrix — exact second derivatives
hessian(f, c(3, 4))

# The D operator composes to any order
D(f, c(3, 4))              # gradient
D(f, c(3, 4), order = 2)   # Hessian
D(f, c(3, 4), order = 3)   # third-order tensor
```

## Use cases

- **Optimization** — supply exact gradients to `optim()` and `nlminb()`
- **Maximum likelihood estimation** — gradients, Hessians, and standard errors
- **Sensitivity analysis** — how outputs change with respect to inputs
- **Taylor approximation** — compute polynomial approximations with exact coefficients
- **Curvature analysis** — second-order geometric properties of curves

## Vignettes

<!-- Links work on the pkgdown site; use browseVignettes("nabla") locally -->
- [Introduction to nabla](articles/introduction.html) — dual numbers, arithmetic, composition, comparison with finite differences
- [MLE Workflow](articles/mle-workflow.html) — gradient, Hessian, Newton-Raphson on statistical models
- [Higher-Order Derivatives](articles/higher-order.html) — the `D` operator, curvature, Taylor expansion
- [Higher-Order MLE Analysis](articles/mle-skewness.html) — third-order derivative tensors and asymptotic skewness of MLEs
- [Optimizer Integration](articles/optimizer-integration.html) — using `gradient()` and `hessian()` with `optim()` and `nlminb()`

## License

MIT
