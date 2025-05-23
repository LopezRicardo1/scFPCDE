README for scFPCDE

Functional PCA-Based Differential Expression for Single-Cell Trajectories
------------------------------------------------------------------------------

The scFPCDE package provides statistical tools to test for differential gene
expression along pseudotime inferred from single-cell RNA-seq data using 
Functional Principal Component Analysis (FPCA).

Installation
------------------------------------------------------------------------------

Install the package directly from GitHub:

    devtools::install_github("LopezRicardo1/scFPCDE")

Example
------------------------------------------------------------------------------

    library(scFPCDE)

    # Load the built-in simulated dataset
    data(scFPCDE_simdata)
    yt <- scale(scFPCDE_simdata$yt)  # Standardize gene expression
    tt <- scFPCDE_simdata$tt

    # Run the differential expression pipeline
    res <- scFPCDE_run(yt, tt)

    # Plot distribution of p-values from D-statistic test
    hist(res$D_test$pval, breaks = 40,
         main = "P-value Distribution (D-test)",
         xlab = "P-value", col = "gray")

    # Plot the smoothed and observed gene expression trajectories
    scFPCDE_gene_curves(
      tt = tt,
      yt = yt,
      yt_fit = res$fpca_result$xt_hat,
      cell_cluster = scFPCDE_simdata$clusters,
      subset = 1:12
    )

License
------------------------------------------------------------------------------

This package is released under the MIT License.
