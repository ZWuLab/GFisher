# GFisher: Generalized Fisher's Combination Tests Under Dependence

## Overview

**GFisher** is an R package that provides accurate and computationally efficient methods for computing
    p-values for a general family of Fisher-type statistics (GFisher), including
    Fisher's combination, Good's statistic, Lancaster's statistic, and weighted
    Z-score combinations. It supports flexible weighting schemes and an omnibus
    procedure that adapts the weights and degrees of freedom to the data. The p-value
    calculation methods used in the testing procedure are based on moment-ratio
    matching and joint distribution approximations. The new version includes both
    pure R implementations and optimized C++ code via RcppArmadillo to improve
    performance on large-scale data. Based on Zhang, H. and Wu, Z. (2022)
    ``The generalized Fisher's combination and accurate p-value calculation under
    dependence", Biometrics, 79(2), 1159-1172.

## Key Features

- **Generalized Fisher's combination tests**: Implements the GFisher family of tests with flexible degrees of freedom and weighting schemes
- **Omnibus testing (oGFisher)**: Adapts test weights and degrees of freedom to the data for increased power
- **Multiple calculation methods**: Supports moment ratio matching (MR), hybrid quadratic approximation (HYB), and Brown's method (GB)
- **Handles dependence**: Accounts for correlation structure among input test statistics through correlation matrices
- **Independent case optimization**: Special fast implementations when input p-values are independent
- **C++ acceleration**: High-performance implementations using RcppArmadillo achieving 2-400x speedups with caching
- **CRAN compliant**: Follows all R package development standards with comprehensive documentation and tests
- **Cross-platform**: Works on Windows, macOS, and Linux

## Installation

### From GitHub (recommended)

```r
# install.packages("remotes")
remotes::install_github("ZWuLab/GFisher")
```

### From source

```r
# requires a C++ compiler
install.packages("devtools")
devtools::install_local("path/to/GFisher")
```

### System Requirements

- R >= 4.3.0
- C++11 compiler
- RcppArmadillo (automatically installed)

### Optional Dependencies

For optimized calculations under independence (`*_ind` functions), install the 'coga' package:

```r
# Check if coga is available
GFisher::check_coga()

# Install coga from CRAN
install.packages("coga")

# Or use the helper function
GFisher::check_coga(install = TRUE)
```

**Note:** If 'coga' is not installed, you can still perform calculations under independence
using the general functions with `M = diag(n)`, though this will be slower than the specialized `*_ind` functions.

## Quick Start

```r
library(GFisher)

# Example: Combine p-values under independence
set.seed(123)
pvals <- c(0.01, 0.03, 0.05, 0.1, 0.2)

# Calculate GFisher statistic (standard Fisher's method: df=2, equal weights)
stat <- stat.GFisher(pvals, df = 2, w = 1)

# Calculate p-value (fast method for independent p-values)
p_value <- p.GFisher_ind_w1(stat, df = rep(2, length(pvals)))
print(p_value)
# Calculate p-value using other functions
p.GFisher_ind(stat, df=rep(2, length(pvals)), w=rep(1, length(pvals)))
p.GFisher(stat, df=2, w=1, M=diag(length(pvals)))

# Example: Account for dependence with correlation matrix
n <- 5
M <- matrix(0.3, n, n) + diag(0.7, n, n)  # Correlation matrix
p_dependent <- p.GFisher(stat, df = 2, w = 1, M = M, method = "HYB")
print(p_dependent)
```

## Package Structure

```
GFisher/
├── R/                  # R function implementations
├── src/                # C++ source code
├── tests/testthat/     # Unit tests
├── vignettes/          # Documentation and examples
└── inst/               # Additional package files
```

## Development Status

**Version 0.3.0**

This version includes fully implemented GFisher methodology with both R and optimized C++ implementations.

### Main Functions

**Basic GFisher Tests:**
- `stat.GFisher()`: Calculate GFisher test statistic
- `p.GFisher()`: Calculate p-value under dependence (methods: HYB, MR, GB)

**Omnibus Tests:**
- `stat.oGFisher()`: Calculate multiple GFisher statistics with different configurations
- `pval.oGFisher()`: Combined p-value using Cauchy combination or minimum p-value

**Independent Case (Optimized, require 'coga' package):**
- `p.GFisher_ind()`: Fast p-value calculation for independent inputs (requires coga)
- `p.GFisher_ind_w1()`: Fastest method for equal weights (no coga needed)
- `stat.oGFisher_ind()`: Omnibus test for independent inputs (requires coga)
- `pval.oGFisher_ind()`: Combined p-value for independent inputs (requires coga)

**Helper Functions:**
- `check_coga()`: Check and optionally install 'coga' package for independent functions
- `getGFisherGM()`: Calculate covariance of GFisher statistic
- `getGFisherlam()`: Calculate eigenvalues for moment matching
- `getGFisherCOR()`: Calculate correlation among multiple GFisher tests
- `getGFishercov()`: Calculate covariance matrix
- `getGFishercoef()`: Calculate coefficients for variance calculation

## Documentation

- **Vignette**: Run `vignette("introduction", package = "GFisher")` for an introduction
- **Function help**: Use `?function_name` for detailed documentation
- **Package help**: Run `?GFisher` for package overview

## Testing

The package includes comprehensive unit tests using testthat:

```r
devtools::test()
```

## Performance Benchmarks

The C++ implementations provide substantial speedups over pure R code:

- **2-10x speedup**: Basic coefficient and covariance calculations
- **50-400x speedup**: With caching enabled for repeated calculations
- **Moment ratio method**: C++ implementation handles large simulations efficiently

**Recent Optimization (2025-10-16):** Refactored `getGFisherGM()` to eliminate code duplication and integrate with the caching system. Both R and C++ versions now share coefficient calculations, providing:
- Eliminated ~80 lines of duplicated code
- 2-4x speedup on repeated calls with same degrees of freedom
- 30-70% time reduction for typical genomic pipelines analyzing thousands of gene sets
- Automatic caching with negligible memory overhead (<2 KB)

## Contributing

For issues or contributions, please file an issue on GitHub or contact the package maintainer.

## License

GPL-3

## Citation

If you use GFisher in your research, please cite:

```
Zhang, H. and Wu, Z. (2022). The generalized Fisher's combination and accurate p-value calculation under dependence. *Biometrics*, 79(2), 1159-1172. doi:10.1111/biom.13634
```

## Contact

Maintainer: Zheyang Wu (zheyangwu@wpi.edu). Please file bugs and questions at https://github.com/ZWuLab/GFisher/issues.

## Authors

- **Hong Zhang** (consistencyzhang@gmail.com) — Author
- **Zheyang Wu** (zheyangwu@wpi.edu) — Author, Maintainer

---

*GFisher is also used as the combination-test backend of the GLOW software family.*

> Developed with AI assistance; see [AI-USE.md](AI-USE.md).
