# dualr 0.5.0

* Generalized to arbitrary-order exact derivatives via recursive nesting.
  New API: `dual_variable_n()`, `dual_constant_n()`, `deriv_n()`,
  `differentiate_n()`.
* Removed Rcpp/C++ fast paths — package is now pure R with no compiled code.
  This simplifies installation and aligns with the package's focus on exact
  arbitrary-order derivatives rather than scalar speed.
* Deprecated second-order-only functions (`dual2_variable()`, `dual2_constant()`,
  `value2()`, `first_deriv()`, `second_deriv()`, `differentiate2()`) as thin
  wrappers around the new generalized API.
* Removed benchmarks (speed is not this package's value proposition).

# dualr 0.4.0

* `score()` now computes the full gradient in 1 forward pass (was p passes)
  using vector-valued derivatives, exploiting the `ANY` slots of the `dualr` class.
* `hessian()` now computes the full Hessian in p forward passes (was p(p+1)/2)
  using vector-gradient inner duals with nested outer duals.
* Internal `.is_scalar_dual()` now also checks `length() == 1L` to correctly
  distinguish scalar duals (C++ fast path) from vector-gradient duals (R path).

# dualr 0.3.0

* Added Rcpp-based C++ fast paths for first-order dual arithmetic (`+`, `-`, `*`, `/`, `^`), math (`exp`, `sqrt`, `log`), and `sum`. Provides 3-10x speedup on scalar dual operations while preserving full R fallback for nested (second-order) duals.
* New internal `.is_scalar_dual()` predicate gates C++ vs R paths using `is.double()` on slot contents.
* Added `Rcpp` to `Imports` and `LinkingTo`; package now requires C++ compilation.

# dualr 0.2.0

* Renamed S4 class from `dual` to `dualr` to avoid conflict with base R's `dual` usage.
* Added dedicated `setMethod` dispatches for hot-path arithmetic (`+`, `-`, `*`, `/`, `^`) and math (`exp`, `sqrt`) operations, bypassing group generic overhead.
* Extracted `.dual_min()` / `.dual_max()` internal helpers, deduplicating 6 inline lambdas across `dual-arithmetic.R` and `dual-math.R`.
* Removed dead `switch` branches (`sqrt`, `exp`, `log`) from `Math` group generic that were shadowed by dedicated methods.
* Standardized `sum()` in `Summary` group generic to use `.as_dual()` promotion, consistent with `prod`, `min`, `max`, and `range`.
* Fixed stale `\code{compositional.mle}` reference in `score()` documentation.

# dualr 0.1.0

* Initial CRAN release.
* S4 dual number class with full arithmetic and math function support.
* Nested duals for exact second-order derivatives (`dual2_variable`, `differentiate2`).
* MLE workflow helpers: `score`, `hessian`, `observed_information`, `score_and_hessian`.
* Special functions: `erf`, `erfc`, `beta`, `lbeta`, `psigamma`.
* Four vignettes: introduction, MLE workflow, higher-order derivatives, optimizer integration.
