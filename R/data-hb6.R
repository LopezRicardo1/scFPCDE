#' HB6 B-Cell Trajectory Inputs from the Manuscript Analysis
#'
#' Log-normalized expression, pseudotime, and B-cell subset labels for the two
#' HB6 trajectories analyzed in the scFPCDE manuscript. These are derived,
#' package-ready inputs extracted from the saved `cds_sub` objects used by the
#' official analysis; Seurat and Monocle objects are not included.
#'
#' @format A named list with components `traj1` and `traj2`. Each trajectory is
#'   a list containing:
#' \describe{
#'   \item{yt}{An uncentered numeric matrix with cells in rows and genes in
#'   columns. Values come from the Monocle `logcounts` assay.}
#'   \item{tt}{A named numeric pseudotime vector aligned to the rows of `yt`.}
#'   \item{clusters}{A named factor of B-cell subset labels aligned to the rows
#'   of `yt`.}
#'   \item{cell_id}{Cell identifiers in matrix-row order.}
#'   \item{gene_id}{Unique gene identifiers in matrix-column order.}
#'   \item{metadata}{Trajectory provenance, the source analysis tag, expression
#'   scale, study identifiers, and the principal manuscript-analysis settings.}
#' }
#'
#' Trajectory 1 contains 2,928 cells and 2,995 unique genes from the Trans,
#' Naive, and C-mem1 subsets. Trajectory 2 contains 1,714 cells and 1,999 genes
#' from the M-mem1, DN2, and DN3 subsets. The duplicated `ACTB` column present
#' in the trajectory-1 paper matrix is represented once so every packaged gene
#' identifier is unique.
#'
#' @usage data(scFPCDE_hb6)
#' @source Stewart A, Ng JC, Wallis G, et al. (2021). Single-Cell
#'   Transcriptomic Analyses Define Distinct Peripheral B Cell Subsets and
#'   Discrete Development Pathways. *Frontiers in Immunology*, 12:602539.
#'   \doi{10.3389/fimmu.2021.602539}. ArrayExpress accession E-MTAB-9544.
#' @references
#' Stewart A, Ng JC, Wallis G, Tsioligka V, Fraternali F, Dunn-Walters DK
#' (2021). Single-Cell Transcriptomic Analyses Define Distinct Peripheral B
#' Cell Subsets and Discrete Development Pathways. *Frontiers in Immunology*,
#' 12:602539. \doi{10.3389/fimmu.2021.602539}.
"scFPCDE_hb6"
