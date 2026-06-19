# GFisher 0.3.0

Initial public release.

- Omnibus generalized Fisher's combination tests (`oGFisher`) that adapt the
  weights and degrees of freedom to the data, alongside the base GFisher family
  (Fisher's combination, Good's, Lancaster's, weighted-Z).
- Accurate p-value calculation under dependence via moment-ratio matching and
  joint-distribution approximations (methods `HYB`, `MR`, `GB`).
- High-performance C++ (RcppArmadillo) implementations with a three-tier
  fallback (cached C++ -> non-cached C++ -> pure R) verified to agree.
- Fast `*_ind` paths for independent inputs via the optional `coga` package,
  with `check_coga()` to detect/install it.

<!-- .release-source: built from monorepo commit 442a15893efc173e9d5f9ac71e78bace78076c0f -->
