# R/label_transfer.R
# Reusable Seurat anchor-based broad cell-type label transfer.
#
# Dependencies: Seurat, SeuratObject. For run_kfold_transfer(), also source
# R/gene_portions.R to make get_transfer_gene_split() available.
#
# The defaults mirror the verified Xenium P10 pilot:
#   FindTransferAnchors(reduction = 'cca', normalization.method = 'LogNormalize',
#                       k.filter = NA, dims = 1:30)
#   TransferData(weight.reduction = 'pca', dims = 1:30)
#
# This module never normalizes expression, never computes a consensus and
# never modifies the caller's reference/query variables. Those are separate
# stages of the pipeline.

#' Check one standardized prediction data frame against the query.
#' @return Invisibly TRUE; error when cell IDs, labels or scores are invalid.
validate_transfer_predictions <- function(predictions, query_cells, celltypes) {
  if (!is.data.frame(predictions)) {
    stop("predictions must be a data.frame.")
  }
  required <- c("cell_id", "predicted.id")
  if (!all(required %in% names(predictions))) {
    stop("Predictions must contain: ", paste(required, collapse = ", "))
  }
  if (!is.character(query_cells) || !length(query_cells) ||
      anyNA(query_cells) || anyDuplicated(query_cells)) {
    stop("query_cells must be unique, nonmissing cell IDs.")
  }
  if (!is.character(celltypes) || !length(celltypes) ||
      anyNA(celltypes) || anyDuplicated(celltypes)) {
    stop("celltypes must be a unique character vector.")
  }
  if (nrow(predictions) != length(query_cells) ||
      !identical(as.character(predictions$cell_id), query_cells)) {
    stop("Prediction cell IDs are missing, reordered or duplicated.")
  }
  labels <- as.character(predictions$predicted.id)
  if (anyNA(labels) || !all(labels %in% celltypes)) {
    stop("Predictions contain missing or unexpected cell-type labels.")
  }
  if ("prediction.score.max" %in% names(predictions)) {
    scores <- predictions$prediction.score.max
    if (!is.numeric(scores) || any(!is.finite(scores)) ||
        any(scores < -1e-8 | scores > 1 + 1e-8)) {
      stop("prediction.score.max must be finite and between 0 and 1.")
    }
  }
  invisible(TRUE)
}

#' Convert TransferData(query=...) output to a compact, ordered data.frame.
#'
#' The output format is compatible with the verified Experiment 4 outputs:
#'   cell_id, predicted.id, prediction.score.max,
#'   prediction.score.<celltype> for each cell type.
#'
#' @param transferred Seurat object returned by TransferData(query = ...).
#' @param celltypes Reference cell types in the desired column order.
#' @param score_assay Score assay returned by transferring named 'id' labels.
extract_transfer_predictions <- function(
    transferred, celltypes, score_assay = "prediction.score.id") {

  if (!inherits(transferred, "Seurat")) {
    stop("transferred must be a Seurat object from TransferData(query=...).")
  }
  cell_ids <- colnames(transferred)
  meta <- transferred[[]]
  required <- c("predicted.id", "predicted.id.score")
  if (!all(required %in% colnames(meta))) {
    stop("TransferData did not produce predicted.id and predicted.id.score metadata.")
  }
  if (!score_assay %in% names(transferred@assays)) {
    stop("Prediction-score assay is missing: ", score_assay)
  }
  scores <- SeuratObject::LayerData(
    transferred, assay = score_assay, layer = "data"
  )
  if (!setequal(rownames(scores), celltypes) ||
      !setequal(colnames(scores), cell_ids)) {
    stop("Prediction-score matrix does not match cell types or query cells.")
  }
  scores <- scores[celltypes, cell_ids, drop = FALSE]
  score_df <- as.data.frame(t(as.matrix(scores)), check.names = FALSE)
  colnames(score_df) <- paste0("prediction.score.", celltypes)

  ordered_meta <- meta[cell_ids, required, drop = FALSE]
  out <- data.frame(
    cell_id = cell_ids,
    predicted.id = as.character(ordered_meta$predicted.id),
    prediction.score.max = as.numeric(ordered_meta$predicted.id.score),
    score_df,
    row.names = NULL,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  validate_transfer_predictions(out, cell_ids, celltypes)
  if (any(!is.finite(as.matrix(score_df)))) {
    stop("Non-finite per-class prediction scores detected.")
  }
  # Maximum class score should match the metadata score for each prediction.
  top_scores <- apply(as.matrix(score_df), 1L, max)
  if (any(abs(top_scores - out$prediction.score.max) > 1e-6)) {
    stop("Prediction maximum scores disagree with the score assay.")
  }
  out
}

#' Find CCA transfer anchors and transfer broad cell-type labels ONCE.
#'
#' @return Compact prediction data.frame, not a Seurat object. An n_anchors
#'   attribute is added for diagnostics. Full Seurat/anchor intermediates
#'   stay local to this function to reduce memory retention.
run_label_transfer <- function(
    reference, query, genes_train,
    label_col = "celltype",
    reference_assay = "RNA", query_assay = "Xenium",
    celltypes = NULL,
    normalization_method = "LogNormalize",
    reduction = "cca", weight_reduction = "pca",
    dims = 1:30, k_filter = NA,
    seed = NULL, verbose = FALSE) {

  if (!inherits(reference, "Seurat") || !inherits(query, "Seurat")) {
    stop("reference and query must both be Seurat objects.")
  }
  if (!reference_assay %in% names(reference@assays) ||
      !query_assay %in% names(query@assays)) {
    stop("Requested reference/query assay is missing.")
  }
  if (!is.character(genes_train) || !length(genes_train) ||
      anyNA(genes_train) || anyDuplicated(genes_train) ||
      !all(genes_train %in% rownames(reference[[reference_assay]])) ||
      !all(genes_train %in% rownames(query[[query_assay]]))) {
    stop("genes_train must be distinct shared genes in both assay objects.")
  }
  if (length(dims) < 1L || !is.numeric(dims) || anyNA(dims) ||
      any(dims < 1L) || any(dims != as.integer(dims)) || anyDuplicated(dims)) {
    stop("dims must be a nonempty vector of unique positive integers.")
  }
  if (!is.null(seed)) {
    if (length(seed) != 1L || is.na(seed) || !is.finite(seed)) {
      stop("seed must be NULL or one finite number.")
    }
    set.seed(seed)
  }
  meta <- reference[[]]
  if (!label_col %in% names(meta)) {
    stop("Reference label field is absent: ", label_col)
  }
  labels <- as.character(meta[colnames(reference), label_col])
  if (anyNA(labels) || any(!nzchar(labels))) {
    stop("Reference label field contains missing or empty labels.")
  }
  names(labels) <- colnames(reference)
  if (is.null(celltypes)) celltypes <- sort(unique(labels))
  if (anyNA(celltypes) || anyDuplicated(celltypes) ||
      !setequal(unique(labels), celltypes)) {
    stop("celltypes must match the unique reference labels.")
  }

  anchors <- Seurat::FindTransferAnchors(
    reference = reference,
    query = query,
    reference.assay = reference_assay,
    query.assay = query_assay,
    features = genes_train,
    normalization.method = normalization_method,
    reduction = reduction,
    k.filter = k_filter,
    dims = dims,
    verbose = verbose
  )
  if (!inherits(anchors, "TransferAnchorSet")) {
    stop("FindTransferAnchors returned an unexpected object type.")
  }
  n_anchors <- if ("anchors" %in% methods::slotNames(anchors)) {
    nrow(methods::slot(anchors, "anchors"))
  } else {
    NA_integer_
  }

  transferred <- Seurat::TransferData(
    anchorset = anchors,
    query = query,
    refdata = labels,
    weight.reduction = weight_reduction,
    dims = dims,
    verbose = verbose
  )

  out <- extract_transfer_predictions(transferred, celltypes = celltypes)
  # Do not assume that Seurat will preserve cell order across future versions.
  # Reorder returned data.frame by query ID if needed, after validating IDs.
  query_cells <- colnames(query)
  if (!setequal(out$cell_id, query_cells)) {
    stop("Transferred predictions are missing query cells.")
  }
  out <- out[match(query_cells, out$cell_id), , drop = FALSE]
  rownames(out) <- NULL
  validate_transfer_predictions(out, query_cells, celltypes)
  attr(out, "n_anchors") <- n_anchors
  out
}

#' Run leave-one-portion-out label transfer for all K gene portions.
#'
#' @param output_dir Optional directory where each run is saved to an RDS.
#' @param resume FALSE by default; TRUE skips previously verified RDS files
#'   only when the experiment manifest matches exactly.
#' @param overwrite FALSE by default; TRUE removes old run RDS files and
#'   starts this directory from scratch (never combine with resume=TRUE).
#' @return Named list of K prediction data.frames (Run_1, ..., Run_K).
#'
#' IMPORTANT: The manifest checks IDs, labels, genes, settings, and Seurat
#' version, but does NOT hash every expression-matrix value. Use a different
#' output directory whenever you change normalization or input expression.
run_kfold_transfer <- function(
    reference, query, gene_portions,
    output_dir = NULL, label_col = "celltype",
    reference_assay = "RNA", query_assay = "Xenium",
    celltypes = NULL,
    normalization_method = "LogNormalize",
    reduction = "cca", weight_reduction = "pca",
    dims = 1:30, k_filter = NA,
    seed = 42L, verbose = FALSE,
    resume = FALSE, overwrite = FALSE) {

  if (!exists("get_transfer_gene_split", mode = "function")) {
    stop("Source R/gene_portions.R before run_kfold_transfer().")
  }
  if (!is.list(gene_portions) || length(gene_portions) < 2L ||
      !identical(names(gene_portions),
                 as.character(seq_along(gene_portions)))) {
    stop("gene_portions must be a named list '1', '2', ... 'K'.")
  }
  if (isTRUE(resume) && isTRUE(overwrite)) {
    stop("resume and overwrite cannot both be TRUE.")
  }
  if (length(seed) != 1L || !is.numeric(seed) ||
      is.na(seed) || !is.finite(seed)) {
    stop("seed must be one finite number.")
  }
  ref_meta <- reference[[]]
  if (!label_col %in% names(ref_meta)) stop("Reference label column is missing.")
  labels <- as.character(ref_meta[colnames(reference), label_col])
  if (is.null(celltypes)) celltypes <- sort(unique(labels))
  if (anyNA(labels) || anyNA(celltypes) ||
      !setequal(unique(labels), celltypes)) {
    stop("Reference cell labels are incomplete or inconsistent with celltypes.")
  }

  # Captures inputs relevant to reproducibility and safe checkpoint reuse.
  manifest <- list(
    reference_cells = colnames(reference),
    reference_labels = labels,
    reference_genes = rownames(reference[[reference_assay]]),
    query_cells = colnames(query),
    query_genes = rownames(query[[query_assay]]),
    gene_portions = gene_portions,
    label_col = label_col,
    reference_assay = reference_assay,
    query_assay = query_assay,
    celltypes = celltypes,
    normalization_method = normalization_method,
    reduction = reduction,
    weight_reduction = weight_reduction,
    dims = dims,
    k_filter = k_filter,
    seed = seed,
    seurat_version = as.character(utils::packageVersion("Seurat"))
  )
  k <- length(gene_portions)
  result <- setNames(vector("list", k), paste0("Run_", seq_len(k)))
  if (!is.null(output_dir)) {
    if (!is.character(output_dir) || length(output_dir) != 1L ||
        is.na(output_dir) || !nzchar(output_dir)) {
      stop("output_dir must be NULL or a single nonempty path.")
    }
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    files <- file.path(output_dir, sprintf("portion_%02d.rds", seq_len(k)))
    mf <- file.path(output_dir, "experiment_manifest.rds")
    has_files <- any(file.exists(files))
    if (isTRUE(resume)) {
      if (file.exists(mf)) {
        previous <- readRDS(mf)
        if (!identical(previous, manifest)) {
          stop("Saved manifest differs from inputs/settings; use a new output_dir.")
        }
      } else if (has_files) {
        stop("Saved predictions exist without a manifest; cannot resume safely.")
      } else {
        saveRDS(manifest, mf)
      }
    } else {
      if ((has_files || file.exists(mf)) && !isTRUE(overwrite)) {
        stop("Output files already exist. Choose a new output_dir, ",
             "overwrite=TRUE, or resume=TRUE with matching manifest.")
      }
      if (isTRUE(overwrite)) unlink(files[file.exists(files)])
      saveRDS(manifest, mf)
    }
  }
  for (i in seq_len(k)) {
    if (!is.null(output_dir) && isTRUE(resume) && file.exists(files[i])) {
      d <- readRDS(files[i])
      validate_transfer_predictions(d, colnames(query), celltypes)
      result[[i]] <- d
      if (verbose) message("Loaded validated Run ", i, " from disk.")
      next
    }
    split <- get_transfer_gene_split(gene_portions, held_out = i)
    if (verbose) {
      message("Run ", i, "/", k, ": ", length(split$genes_train),
              " training genes; ", length(split$genes_test), " held out.")
    }
    d <- run_label_transfer(
      reference = reference, query = query,
      genes_train = split$genes_train,
      label_col = label_col,
      reference_assay = reference_assay,
      query_assay = query_assay,
      celltypes = celltypes,
      normalization_method = normalization_method,
      reduction = reduction,
      weight_reduction = weight_reduction,
      dims = dims, k_filter = k_filter,
      seed = seed + i, verbose = verbose
    )
    validate_transfer_predictions(d, colnames(query), celltypes)
    result[[i]] <- d
    if (!is.null(output_dir)) saveRDS(d, files[i])
    rm(d)
    gc(verbose = FALSE)
  }
  result
}
