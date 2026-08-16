# scFPCDE

`scFPCDE` tests for differential gene expression along single-cell pseudotime
using functional principal component analysis (FPCA). The package provides
functions for fitting smooth expression trajectories, tuning smoothing
parameters, running permutation-based D- and F-tests, and visualizing fitted
gene curves.

## Installation

Install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("NLM-DIR/scFPCDE", build_vignettes = TRUE)
```

## Quick start

The bundled simulation contains 1,000 cells and 500 genes. For a quick example,
the code below uses 20 differentially expressed genes and 40 null genes and a
small number of permutations. Increase `n_perm` for an analysis intended for
inference.

```r
library(scFPCDE)

data(scFPCDE_simdata)
gene_index <- c(1:20, 101:140)
cell_index <- seq(1, length(scFPCDE_simdata$tt), length.out = 200)
cell_index <- unique(round(cell_index))

yt <- scFPCDE_simdata$yt[cell_index, gene_index]
tt <- scFPCDE_simdata$tt[cell_index]

set.seed(2026)
res <- scFPCDE_run(
  yt = yt,
  tt = tt,
  n_perm = 25,
  ncores = 1
)

head(res$D_test_result)
hist(
  res$D_test_result$p_value,
  breaks = 20,
  main = "D-test permutation p-values",
  xlab = "p-value"
)

scFPCDE_gene_curves(
  tt = tt,
  yt = scale(yt, center = TRUE, scale = FALSE),
  yt_fit = res$fpca_result$xt_hat,
  cell_cluster = scFPCDE_simdata$clusters[cell_index],
  subset = order(res$D_test_result$D_obs, decreasing = TRUE)[1:12]
)
```

`scFPCDE_run()` returns three named components:

- `fpca_result`: the fitted FPCA model and trajectories;
- `D_test_result`: gene-level D statistics, p-values, and BH-adjusted q-values;
- `F_test_result`: gene-level F-test results when `use_FPC_F = TRUE`, otherwise
  `NULL`.

See `vignette("scFPCDE-overview")` for a guided workflow and the function help
pages for parameter details.

## Citation

Run `citation("scFPCDE")` to obtain the package citation. Publication metadata
can be added when the associated manuscript citation is final.

## License

This package is released under the MIT License.
