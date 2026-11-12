# R/reference.R
# CELLxGENE SingleCellExperiment loading, reference-cell selection, and
# conversion of an existing expression matrix to a Seurat reference.
#
# Dependencies: zellkonverter, SingleCellExperiment,
# SummarizedExperiment, Matrix, Seurat, SeuratObject.
#
# Important: this module does NOT normalize the reference. It copies the
# requested .h5ad assay (normally "X") into the Seurat RNA/data layer.
# The user must confirm that the source assay is appropriate for transfer.

#' Load a CELLxGENE .h5ad reference as a SingleCellExperiment.
#'
#' @param path Path to an .h5ad file.
#' @param assay_name Name of the expression assay to use, usually "X".
#' @param label_col Column holding the known broad cell-type labels.
#' @param gene_col Row-data column holding gene symbols.
#' @param use_hdf5 Keep expression HDF5-backed, instead of fully loading it.
#' @return SingleCellExperiment, with original metadata and assays intact.
load_scrna_reference <- function(path, assay_name = "X",
                                 label_col = "author_cell_type",
                                 gene_col = "gene_symbols",
                                 use_hdf5 = TRUE) {
  if (!file.exists(path)) stop("Reference file does not exist: ", path)
  sce <- zellkonverter::readH5AD(
    path, reader = "R", use_hdf5 = use_hdf5
  )
  if (!inherits(sce, "SingleCellExperiment")) {
    stop("readH5AD did not return a SingleCellExperiment.")
  }
  if (!assay_name %in% SummarizedExperiment::assayNames(sce)) {
    stop("Expression assay not found: ", assay_name)
  }
  if (!label_col %in% colnames(SummarizedExperiment::colData(sce))) {
    stop("Reference label column not found: ", label_col)
  }
  if (!gene_col %in% colnames(SummarizedExperiment::rowData(sce))) {
    stop("Reference gene column not found: ", gene_col)
  }
  sce
}

#' Select reference cell indices without accessing the expression matrix.
#'
#' @param sce Full SingleCellExperiment.
#' @param method One of "per_type", "per_state", "balanced_within_type".
#' @param n_per_type Number of cells per broad type for per_type or
#'   balanced_within_type methods.
#' @param max_per_state Maximum cells in each (broad type, state) group for
#'   per_state method.
#' @param celltypes Allowed broad types, in desired order.
#' @param label_col Broad cell-type column in colData(sce).
#' @param state_col Detailed state column in colData(sce).
#' @param seed Sampling seed.
#' @param preserve Optional existing numeric column indices for nested per_state
#'   sampling. All provided indices must fit the new cap and eligible groups.
#' @return Integer column positions into the ORIGINAL sce object.
select_reference_cells <- function(
    sce,
    method = c("per_type", "per_state", "balanced_within_type"),
    n_per_type = NULL,
    max_per_state = NULL,
    celltypes = c("T", "B", "Mye", "Fibro", "Endo", "Peri", "Tumor"),
    label_col = "author_cell_type",
    state_col = "cell_state",
    seed = 42L,
    preserve = NULL) {

  method <- match.arg(method)
  if (!inherits(sce, "SingleCellExperiment")) {
    stop("sce must be a SingleCellExperiment.")
  }
  cd <- SummarizedExperiment::colData(sce)
  if (!label_col %in% colnames(cd)) stop("Missing label column: ", label_col)
  if (anyNA(celltypes) || !length(celltypes) || anyDuplicated(celltypes)) {
    stop("celltypes must contain distinct, nonmissing types.")
  }
  types <- as.character(cd[[label_col]])
  if (length(seed) != 1L || is.na(seed)) stop("seed must be a single value.")
  set.seed(seed)

  choose <- function(ids, n) ids[sample.int(length(ids), size = n)]
  is_positive_integer <- function(x) {
    length(x) == 1L && !is.na(x) && is.numeric(x) &&
      is.finite(x) && x >= 1 && x == floor(x)
  }

  if (method == "per_type") {
    if (!is_positive_integer(n_per_type)) {
      stop("n_per_type must be a positive integer for per_type.")
    }
    if (!is.null(preserve)) stop("preserve is only supported for per_state.")
    ans <- unlist(lapply(celltypes, function(ct) {
      ids <- which(!is.na(types) & types == ct)
      if (length(ids) < n_per_type) {
        stop("Not enough cells of type ", ct, " for n_per_type=", n_per_type)
      }
      choose(ids, n_per_type)
    }), use.names = FALSE)
    return(as.integer(ans))
  }

  if (!state_col %in% colnames(cd)) stop("Missing state column: ", state_col)
  states <- as.character(cd[[state_col]])
  valid <- !is.na(types) & types %in% celltypes &
    !is.na(states) & nzchar(states)

  if (method == "per_state") {
    if (!is_positive_integer(max_per_state)) {
      stop("max_per_state must be a positive integer for per_state.")
    }
    valid_ids <- which(valid)
    group_ids <- paste(types[valid_ids], states[valid_ids], sep = "::")
    groups <- split(valid_ids, group_ids)
    if (!length(groups)) stop("No eligible (type, state) groups.")

    if (!is.null(preserve)) {
      if (!is.numeric(preserve) || anyNA(preserve) ||
          any(preserve != floor(preserve)) || anyDuplicated(preserve) ||
          !all(preserve %in% valid_ids)) {
        stop("preserve must contain unique eligible column indices into sce.")
      }
    }

    ans <- unlist(lapply(groups, function(ids) {
      target <- min(length(ids), max_per_state)
      existing <- if (is.null(preserve)) integer(0) else intersect(ids, preserve)
      if (length(existing) > target) {
        stop("A preserved state group already exceeds max_per_state.")
      }
      remaining <- setdiff(ids, existing)
      n_add <- target - length(existing)
      additional <- if (n_add > 0L) choose(remaining, n_add) else integer(0)
      c(existing, additional)
    }), use.names = FALSE)
    if (!is.null(preserve) && !all(preserve %in% ans)) {
      stop("Some preserved cells were not retained.")
    }
    return(as.integer(ans))
  }

  # balanced_within_type: approximately equal allocation across detailed
  # states within EACH broad type, keeping n_per_type identical by type.
  if (!is_positive_integer(n_per_type)) {
    stop("n_per_type must be a positive integer for balanced_within_type.")
  }
  if (!is.null(preserve)) stop("preserve is only supported for per_state.")

  ans <- unlist(lapply(celltypes, function(ct) {
    ids <- which(valid & types == ct)
    if (length(ids) < n_per_type) {
      stop("Not enough state-annotated cells of type ", ct)
    }
    state_groups <- split(ids, states[ids])
    available <- lengths(state_groups)
    quota <- pmin(available, floor(n_per_type / length(state_groups)))
    while (sum(quota) < n_per_type) {
      eligible <- which(quota < available)
      n_add <- min(n_per_type - sum(quota), length(eligible))
      chosen <- eligible[sample.int(length(eligible), n_add)]
      quota[chosen] <- quota[chosen] + 1L
    }
    unlist(Map(choose, state_groups, as.integer(quota)), use.names = FALSE)
  }), use.names = FALSE)
  as.integer(ans)
}

#' Construct a Seurat RNA/data-only reference from an existing assay.
#'
#' @param sce Full SingleCellExperiment with existing expression data.
#' @param selected_cells Numeric column positions, or cell IDs, to retain.
#' @param genes Unique gene SYMBOLS to keep (e.g. shared_genes).
#' @param expression_assay Existing SCE assay containing normalized values.
#' @param gene_col SCE rowData column with gene symbols.
#' @param label_col SCE colData column with known broad types.
#' @param seurat_assay Desired Seurat assay name.
#' @param project Seurat project label.
#' @param verify_values Compare a small sample of values to the original SCE.
#' @return Seurat object with existing values in RNA/data; no raw counts layer.
build_seurat_reference <- function(
    sce, selected_cells, genes,
    expression_assay = "X", gene_col = "gene_symbols",
    label_col = "author_cell_type", seurat_assay = "RNA",
    project = "TNBC_reference", verify_values = TRUE) {

  if (!inherits(sce, "SingleCellExperiment")) {
    stop("sce must be a SingleCellExperiment.")
  }
  if (!expression_assay %in% SummarizedExperiment::assayNames(sce)) {
    stop("Expression assay not found: ", expression_assay)
  }
  if (!gene_col %in% colnames(SummarizedExperiment::rowData(sce))) {
    stop("Gene column not found: ", gene_col)
  }
  if (!label_col %in% colnames(SummarizedExperiment::colData(sce))) {
    stop("Label column not found: ", label_col)
  }
  if (!is.character(genes) || !length(genes) ||
      anyNA(genes) || any(!nzchar(genes)) || anyDuplicated(genes)) {
    stop("genes must be unique, nonempty, nonmissing gene symbols.")
  }
  if (is.character(selected_cells)) {
    if (anyNA(selected_cells) || anyDuplicated(selected_cells)) {
      stop("selected_cells must be unique and nonmissing.")
    }
    selected_cells <- match(selected_cells, colnames(sce))
  }
  if (!is.numeric(selected_cells) || !length(selected_cells) ||
      anyNA(selected_cells) || anyDuplicated(selected_cells) ||
      any(selected_cells != floor(selected_cells)) ||
      any(selected_cells < 1L | selected_cells > ncol(sce))) {
    stop("selected_cells must be valid unique column positions or cell IDs.")
  }
  selected_cells <- as.integer(selected_cells)

  ref_symbols <- as.character(SummarizedExperiment::rowData(sce)[[gene_col]])
  gene_idx <- match(genes, ref_symbols)
  if (anyNA(gene_idx)) stop("Some requested gene symbols are absent.")
  if (sum(ref_symbols %in% genes, na.rm = TRUE) != length(genes)) {
    stop("Some requested gene symbols are duplicated in reference metadata.")
  }

  subset_sce <- sce[gene_idx, selected_cells, drop = FALSE]
  rownames(subset_sce) <- genes
  subset_sce$celltype <- as.character(
    SummarizedExperiment::colData(subset_sce)[[label_col]]
  )
  if (anyNA(subset_sce$celltype) || any(!nzchar(subset_sce$celltype))) {
    stop("Selected reference cells include missing broad-type labels.")
  }
  SingleCellExperiment::mainExpName(subset_sce) <- seurat_assay

  # Avoid an explicit dense as.matrix() of the full HDF5-backed assay.
  sparse_expr <- as(
    SummarizedExperiment::assay(subset_sce, expression_assay),
    "dgCMatrix"
  )
  dimnames(sparse_expr) <- dimnames(subset_sce)
  SummarizedExperiment::assay(subset_sce, expression_assay) <- sparse_expr
  rm(sparse_expr)

  sc <- Seurat::as.Seurat(
    subset_sce, counts = NULL, data = expression_assay, project = project
  )
  SeuratObject::DefaultAssay(sc) <- seurat_assay
  if (!identical(rownames(sc), genes) ||
      !identical(colnames(sc), colnames(sce)[selected_cells]) ||
      !identical(as.character(sc$celltype), as.character(subset_sce$celltype))) {
    stop("Conversion changed gene/cell order or cell-type labels.")
  }
  if (verify_values) {
    # This comparison is sampled, not exhaustive.
    test_gene_pos <- unique(round(seq(1, length(genes), length.out = min(20L, length(genes)))))
    test_cell_pos <- unique(round(seq(1, length(selected_cells),
                                      length.out = min(25L, length(selected_cells)))))
    original_test <- as.matrix(
      SummarizedExperiment::assay(sce, expression_assay)[
        gene_idx[test_gene_pos], selected_cells[test_cell_pos], drop = FALSE
      ]
    )
    converted_test <- as.matrix(
      SeuratObject::LayerData(sc, assay = seurat_assay, layer = "data")[
        genes[test_gene_pos], colnames(sc)[test_cell_pos], drop = FALSE
      ]
    )
    if (!isTRUE(all.equal(unname(original_test), unname(converted_test),
                          tolerance = 1e-12))) {
      stop("Expression values did not match the original SCE assay.")
    }
  }
  sc
}
