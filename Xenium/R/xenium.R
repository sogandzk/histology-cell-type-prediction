# R/xenium.R
# Xenium loading, preprocessing and query-cell selection.
# Dependencies: Seurat, SeuratObject and Matrix.
# These functions reproduce the original exploratory Xenium preprocessing
# (remove zero-count cells; select genes detected in >= 1 cell;
# LogNormalize; ScaleData; RunPCA).

#' Load one Xenium sample and exclude zero-transcript cells.
#'
#' @param data_dir Folder containing standard Xenium output.
#' @param assay Assay to create (default: "Xenium").
#' @param fov Field-of-view name (default: "fov").
#' @param remove_zero Whether to exclude cells with zero counts.
#' @return A Seurat object (not normalized by this function).
load_xenium <- function(data_dir, assay = "Xenium", fov = "fov",
                        remove_zero = TRUE) {
  if (!dir.exists(data_dir)) {
    stop("Xenium output directory does not exist: ", data_dir)
  }

  obj <- Seurat::LoadXenium(
    data.dir = data_dir,
    fov = fov,
    assay = assay
  )
  SeuratObject::DefaultAssay(obj) <- assay

  if (remove_zero) {
    count_column <- paste0("nCount_", assay)
    meta <- obj[[]]
    if (!count_column %in% colnames(meta)) {
      stop("Expected cell metadata column missing: ", count_column)
    }
    ncounts <- meta[colnames(obj), count_column]
    keep <- colnames(obj)[is.finite(ncounts) & ncounts > 0]
    if (!length(keep)) stop("No cells with positive transcript counts.")
    obj <- subset(obj, cells = keep)
  }

  obj
}

#' Prepare a Xenium Seurat object for the CCA label-transfer workflow.
#'
#' @param obj Seurat object containing a raw Xenium counts layer.
#' @param assay Xenium assay name.
#' @param npcs Number of principal components (default: 50).
#' @param scale_factor Library-size normalization factor (default: 10000).
#' @param pca_seed Seed passed to RunPCA (default: 42).
#' @param verbose Show preprocessing output.
#' @return A Seurat object with counts, data, scale.data and PCA.
prepare_xenium_query <- function(obj, assay = "Xenium", npcs = 50L,
                                 scale_factor = 10000, pca_seed = 42L,
                                 verbose = FALSE) {
  if (!inherits(obj, "Seurat")) stop("obj must be a Seurat object.")
  if (length(npcs) != 1L || is.na(npcs) || npcs < 1L || npcs != as.integer(npcs)) {
    stop("npcs must be a positive integer.")
  }
  SeuratObject::DefaultAssay(obj) <- assay

  counts <- SeuratObject::LayerData(obj, assay = assay, layer = "counts")
  if (is.null(counts) || !nrow(counts) || !ncol(counts)) {
    stop("A nonempty raw counts layer is required in assay '", assay, "'.")
  }

  # This matches the original script: treat all genes expressed in >= 1 cell
  # as variable features. We are NOT running FindVariableFeatures here.
  expressed_genes <- rownames(counts)[Matrix::rowSums(counts) > 0]
  if (length(expressed_genes) < 2L) stop("Too few detected genes for PCA.")
  if (ncol(obj) < 2L) stop("At least two query cells are required for PCA.")
  if (npcs >= min(length(expressed_genes), ncol(obj))) {
    stop("npcs must be smaller than both the gene and cell counts.")
  }
  SeuratObject::VariableFeatures(obj) <- expressed_genes

  obj <- Seurat::NormalizeData(
    obj, assay = assay, normalization.method = "LogNormalize",
    scale.factor = scale_factor, verbose = verbose
  )
  obj <- Seurat::ScaleData(
    obj, assay = assay, features = expressed_genes, verbose = verbose
  )
  obj <- Seurat::RunPCA(
    obj, assay = assay, features = expressed_genes,
    npcs = as.integer(npcs), seed.use = pca_seed, verbose = verbose
  )
  obj
}

#' Select a pilot query, or return the complete query unmodified.
#'
#' IMPORTANT: This function does NOT recalculate normalization, scaling, or PCA.
#' If supplied a preprocessed object, it carries over that preprocessing.
#'
#' @param obj Seurat query (prepared or unprepared).
#' @param cell_ids Explicit cell IDs to retain; use to reproduce a prior pilot.
#' @param n_cells Number of cells to sample at random instead.
#' @param seed Sampling seed used only when n_cells is specified.
#' @return A Seurat object containing the chosen cells.
select_xenium_query <- function(obj, cell_ids = NULL, n_cells = NULL,
                                seed = 42L) {
  if (!inherits(obj, "Seurat")) stop("obj must be a Seurat object.")
  if (!is.null(cell_ids) && !is.null(n_cells)) {
    stop("Specify cell_ids OR n_cells, not both.")
  }
  if (is.null(cell_ids) && is.null(n_cells)) return(obj)

  available <- colnames(obj)
  if (!is.null(n_cells)) {
    if (length(n_cells) != 1L || is.na(n_cells) ||
        n_cells < 1L || n_cells != as.integer(n_cells) ||
        n_cells > length(available)) {
      stop("n_cells must be an integer between 1 and ", length(available), ".")
    }
    set.seed(seed)
    cell_ids <- sample(available, size = n_cells)
  }

  if (!is.character(cell_ids) || !length(cell_ids) ||
      anyNA(cell_ids) || anyDuplicated(cell_ids)) {
    stop("cell_ids must be a nonempty vector of unique, nonmissing cell IDs.")
  }
  if (!all(cell_ids %in% available)) {
    stop("Some requested cell IDs are absent from the query.")
  }
  selected <- subset(obj, cells = cell_ids)
  stopifnot(setequal(colnames(selected), cell_ids))
  selected
}
