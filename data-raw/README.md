# Rebuilding the HB6 example data

`scFPCDE_hb6.R` extracts the two exact cell and gene subsets used by the
manuscript analysis from `res_traj1.rds` and `res_traj2.rds`. It reads the
uncentered `logcounts` assay from each saved `cds_sub` object and writes only
plain R matrices, vectors, and provenance metadata to the package data file.
The large Monocle and Seurat objects are not included in the package.

From the package root, run:

```sh
Rscript data-raw/scFPCDE_hb6.R /path/to/bcell_two_traj
```

The source directory is intentionally supplied at build time so that no
machine-specific path is embedded in the repository. The generated object is
validated for cell alignment, unique gene names, finite values, and missing
cluster labels before it is saved.
