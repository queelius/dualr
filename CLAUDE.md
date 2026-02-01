# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Is

`dualr` is an R package for forward-mode automatic differentiation using dual numbers. It provides exact derivatives at machine precision through any R code (loops, branches, control flow) via operator overloading on S4 classes.

## Build & Test Commands

```bash
# Full R CMD check (matches CI)
R CMD check . --no-manual --compact-vignettes=gs+qpdf

# Run all tests
Rscript -e 'devtools::test()'

# Run a single test file
Rscript -e 'testthat::test_file("tests/testthat/test-arithmetic.R")'

# Regenerate NAMESPACE and man/ pages after changing roxygen comments
Rscript -e 'roxygen2::roxygenise()'

# Test coverage
Rscript -e 'covr::package_coverage()'
Rscript -e 'covr::report()'

# Build vignettes
Rscript -e 'devtools::build_vignettes()'
```

## Architecture

### S4 Class Hierarchy

Two S4 classes defined in `R/dual-class.R`:

- **`dualr`**: Core dual number with `value` (ANY) and `deriv` (ANY) slots. Both slots accept numeric OR another `dualr` object, enabling nesting for higher-order derivatives.
- **`dual_vector`**: List-based container of `dualr` objects for multi-parameter functions. Supports `[` and `[[` indexing so user code can write `theta[1]`, `theta[2]`.

### Source File Roles (Collate Order Matters)

The `Collate:` field in DESCRIPTION defines load order — dependencies flow top to bottom:

1. `dual-class.R` — Class definitions, constructors (`dual()`, `dual_variable()`, `dual_constant()`), accessors (`value()`, `deriv()`), coercion, display
2. `dual-arithmetic.R` — Arithmetic operators (+, -, *, /, ^) with derivative rules; Ops/Summary group generics
3. `dual-math.R` — Math/Math2 group generics (trig, exp, log, gamma, etc.); standalone `atan2()`, `log(x, base)`
4. `dual-special.R` — Non-base-generic functions: `erf()`, `erfc()`, `beta()`, `lbeta()`, `psigamma()`
5. `dual-higher.R` — Second-order via nested duals: `dual2_variable()`, `differentiate2()`
6. `mle-helpers.R` — High-level API: `score()`, `hessian()`, `observed_information()`
7. `dualr-package.R` — Package-level roxygen docs

### Method Dispatch Pattern

- **Hot-path operations** (`+`, `-`, `*`, `/`, `^`, `exp`, `sqrt`) have direct S4 methods for performance — they bypass group generic switch overhead.
- **Everything else** uses group generics (`Ops`, `Math`, `Math2`, `Summary`) with a `switch(.Generic, ...)` fallback.
- Each binary operator typically has three signatures: `(dualr,dualr)`, `(dualr,numeric)`, `(numeric,dualr)`.
- `.as_dual()` promotes plain numerics to dual constants (deriv=0) for mixed operations.

### Derivative Propagation

The fundamental identity: `f(a + b*ε) = f(a) + f'(a)*b*ε` where `ε² = 0`.

- Arithmetic uses product rule, quotient rule, power rule (with `.is_zero()` optimization for constant exponents)
- Math functions apply chain rule: `deriv(f(x)) = f'(value(x)) * deriv(x)`
- Second-order derivatives nest duals: `dual(dual(x, 1), dual(1, 0))` — after evaluation, `deriv(deriv(result))` gives f''(x)

### MLE Workflow

`score()` runs p forward passes (one per parameter, seeding deriv=1 on each in turn). `hessian()` uses nested duals with p*(p+1)/2 passes exploiting symmetry. Both use internal `.make_dual_vector()` / `.make_dual2_vector()` for seeding.

## Testing Conventions

- Framework: testthat 3rd edition
- Test files mirror source structure: `test-arithmetic.R`, `test-math.R`, `test-special.R`, `test-higher-order.R`, `test-mle-helpers.R`
- `test-coverage.R` targets uncovered edge cases specifically
- `test-optimizer-integration.R` tests AD gradients with `optim()` and `nlminb()`
- `tests/testthat/helper-numerical.R` provides `central_difference()`, `numerical_gradient()`, `numerical_hessian()` for verification
- Tolerances: 1e-6 to 1e-10 for first-order, 1e-5 for second-order (numerical stability)
- Standard verification pattern: compare AD result against both analytical formula and numerical finite difference

## CI

GitHub Actions runs R CMD check on macOS (release), Windows (release), Ubuntu (release + devel). Config in `.github/workflows/R-CMD-check.yaml`.

## Dependencies

Hard: `R (>= 4.0.0)`, `methods`. Import: `stats::pnorm` (for erf/erfc). Suggests: `testthat`, `knitr`, `rmarkdown`.
