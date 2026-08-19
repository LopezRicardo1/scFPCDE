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

## Included example data

`scFPCDE_hb6` contains the two analysis-ready HB6 B-cell trajectories used in
the manuscript: uncentered log-normalized expression, aligned raw counts,
pseudotime, cluster labels, and provenance metadata. The source Seurat and
Monocle objects are not bundled. `scFPCDE_simdata` is a separate known-truth
simulation for method checks and teaching.

```r
data(scFPCDE_hb6)
hb6 <- scFPCDE_hb6$traj1

stopifnot(
  identical(rownames(hb6$yt), names(hb6$tt)),
  identical(rownames(hb6$yt), names(hb6$clusters))
)

genes_kept <- scFPCDE_filter_genes(hb6$yt, qz = 0.05)
hb6$yt[, genes_kept, drop = FALSE]
```

The complete real-data tutorial in `vignette("scFPCDE-overview")` demonstrates
fixed-basis roughness tuning, the GCV curve, 100-permutation D- and F-tests,
FPC score visualization, and fitted curves on the original log-expression
scale.

## Quick simulation check

The bundled simulation contains 1,000 cells and 500 genes. For a quick example,
the code below uses 20 differentially expressed genes, 40 null genes, and 100
permutations. Increase `n_perm` to at least 1,000 for an analysis intended for
inference and examine stability across random seeds.

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
  n_perm = 100,
  ncores = 1
)

head(res$D_test_result)
scores_for_plot <- scFPCDE_fpc_scores(
  res,
  components = 1:2,
  transform = "signed_log10"
)
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
  subset = order(res$D_test_result$D_obs, decreasing = TRUE)[1:12],
  point_size = 0.45,
  point_alpha = 0.55,
  curve_linewidth = 0.9
)
```

`scFPCDE_run()` returns three named components:

- `fpca_result`: the fitted FPCA model and trajectories;
- `D_test_result`: gene-level D statistics, p-values, and BH-adjusted q-values;
- `F_test_result`: gene-level F-test results when `use_FPC_F = TRUE`, otherwise
  `NULL`.

Use `scFPCDE_fpc_scores()` to extract raw gene-level FPC scores or create
signed-log10-compressed and standardized coordinates for visualization. Keep
the raw scores for inference, distance calculations, and null boundaries.

See `vignette("scFPCDE-overview")` for the guided real-data and known-truth
workflows, and the function help pages for parameter details.

## Citation

Run `citation("scFPCDE")` to obtain the package citation. Publication metadata
can be added when the associated manuscript citation is final.

## License

This package is released under the MIT License.
