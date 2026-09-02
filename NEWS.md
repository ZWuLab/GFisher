# GFisher 0.3.1

- The compiled code now links BLAS and LAPACK explicitly, through `src/Makevars`
  and `src/Makevars.win`. The C++ sources use RcppArmadillo, whose matrix
  operations call BLAS and LAPACK routines; R links a BLAS implicitly on Linux
  and macOS, but the Windows (Rtools) toolchain does not, so without this the
  Windows build failed at link time with undefined references to `dgemm_` and
  similar. There is no change to any computation.

  This fix has been on the repository's `main` branch since 2026-06-19, but no
  released tag carried it: `v0.3.0` points at a tree published before the fix.
  Installing `ZWuLab/GFisher@v0.3.1` now gets it.
- The package maintainer is now Zheyang Wu (zheyangwu@wpi.edu). Hong Zhang
  remains an author. Bug reports go to https://github.com/ZWuLab/GFisher/issues.

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
