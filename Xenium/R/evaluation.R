# R/evaluation.R
# Reusable evaluation of Xenium broad-cell-type consensus annotations against
# deposited author annotations. No Seurat or other package dependency.
#
# Important: agreement with author-provided labels is a reproduction metric,
# not independent biological ground-truth accuracy. Author-unlabeled cells
# may have been removed by additional processing or confidence filters.
#
# Inputs: consensus data.frames from calculate_consensus(), and author_meta
# with unique cell identifiers (default: cellname) and labels (cell_type).
#
# Outputs align by cell ID, never by row position. Predictions outside the
# allowed celltypes (including LOWCONF) are rejected, except LOWCONF itself.

.check_eval_celltypes <- function(celltypes) {
  if (!is.character(celltypes) || length(celltypes) == 0L ||
      anyNA(celltypes) || any(!nzchar(celltypes)) ||
      anyDuplicated(celltypes) || "LOWCONF" %in% celltypes) {
    stop("celltypes must be unique, nonempty class names excluding LOWCONF.")
  }
  invisible(TRUE)
}

.check_eval_ids <- function(ids, label) {
  if (!is.character(ids) || length(ids) == 0L ||
      anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop(label, " must be unique, nonempty, nonmissing character IDs.")
  }
  invisible(TRUE)
}

.check_eval_comparison <- function(comparison, celltypes, prediction_cols) {
  .check_eval_celltypes(celltypes)
  if (!is.data.frame(comparison) ||
      !all(c("cell_id", "author") %in% names(comparison)) ||
      nrow(comparison) == 0L) {
    stop("comparison must have nonempty cell_id and author columns.")
  }
  .check_eval_ids(as.character(comparison$cell_id), "comparison$cell_id")
  if (!is.character(prediction_cols) || length(prediction_cols) == 0L ||
      anyNA(prediction_cols) || anyDuplicated(prediction_cols) ||
      !all(prediction_cols %in% names(comparison)) ||
      any(prediction_cols %in% c("cell_id", "author"))) {
    stop("prediction_cols must name distinct prediction columns in comparison.")
  }
  for (col in prediction_cols) {
    values <- as.character(comparison[[col]])
    if (anyNA(values) || !all(values %in% c(celltypes, "LOWCONF"))) {
      stop("Unexpected or missing prediction labels in column: ", col)
    }
  }
  invisible(TRUE)
}

#' Match multiple consensus results to the same author-labeled query cells.
#'
#' @param consensus_list Named list of consensus data.frames, one per reference
#'   experiment. Each must have cell_id and css_id_full columns.
#' @param author_meta Data.frame with author cell IDs and annotations.
#' @param celltypes Vector of permitted confident broad cell types.
#' @param cell_ids Optional cell IDs in desired output order; defaults to the
#'   first consensus table's cell_id order. All consensus datasets must contain
#'   exactly these cells (order within each table may differ).
#' @param author_id_col Column in author_meta containing cell IDs.
#' @param author_label_col Column containing author labels. NA is allowed here
#'   and means unlabeled during evaluation.
#' @param consensus_label_col Column with final predictions (default css_id_full).
#' @return data.frame with cell_id, author, and one column per named experiment.
make_annotation_comparison <- function(
    consensus_list, author_meta, celltypes, cell_ids = NULL,
    author_id_col = "cellname", author_label_col = "cell_type",
    consensus_label_col = "css_id_full") {

  .check_eval_celltypes(celltypes)
  if (!is.list(consensus_list) || length(consensus_list) == 0L ||
      is.null(names(consensus_list)) ||
      anyNA(names(consensus_list)) ||
      any(!nzchar(names(consensus_list))) ||
      anyDuplicated(names(consensus_list)) ||
      any(names(consensus_list) %in% c("cell_id", "author"))) {
    stop("consensus_list must be a named list with unique experiment names.")
  }
  if (!is.data.frame(author_meta) ||
      !all(c(author_id_col, author_label_col) %in% names(author_meta))) {
    stop("author_meta is missing its cell-ID or annotation column.")
  }
  author_ids <- as.character(author_meta[[author_id_col]])
  .check_eval_ids(author_ids, paste0("author_meta$", author_id_col))

  if (is.null(cell_ids)) {
    first <- consensus_list[[1L]]
    if (!is.data.frame(first) || !"cell_id" %in% names(first)) {
      stop("First consensus table must contain cell_id.")
    }
    cell_ids <- as.character(first$cell_id)
  }
  .check_eval_ids(cell_ids, "cell_ids")
  author_idx <- match(cell_ids, author_ids)
  if (anyNA(author_idx)) {
    stop(sum(is.na(author_idx)), " query cell IDs are missing in author_meta.")
  }

  result <- data.frame(
    cell_id = cell_ids,
    author = as.character(author_meta[[author_label_col]][author_idx]),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  for (nm in names(consensus_list)) {
    d <- consensus_list[[nm]]
    if (!is.data.frame(d) ||
        !all(c("cell_id", consensus_label_col) %in% names(d))) {
      stop("Consensus dataset ", nm, " must contain cell_id and ", consensus_label_col, ".")
    }
    ids <- as.character(d$cell_id)
    .check_eval_ids(ids, paste0("consensus_list$", nm, "$cell_id"))
    if (length(ids) != length(cell_ids) || !setequal(ids, cell_ids)) {
      stop("Consensus dataset ", nm, " does not have exactly the query cell IDs.")
    }
    labels <- as.character(d[[consensus_label_col]][match(cell_ids, ids)])
    if (anyNA(labels) || !all(labels %in% c(celltypes, "LOWCONF"))) {
      stop("Consensus dataset ", nm, " contains missing or unexpected labels.")
    }
    result[[nm]] <- labels
  }
  rownames(result) <- NULL
  result
}

#' Evaluate each experiment against the author-annotated cells.
#'
#' - confident_cells: number of predictions that are one of celltypes.
#' - coverage_pct: confident_cells / total query cells * 100.
#' - jointly_labeled: cells labeled by both author and experiment.
#' - matching_labels: correct matches among jointly labeled cells.
#' - agreement_pct: matching_labels / jointly_labeled * 100.
#' - author_recovery_pct: matching_labels / all author-labeled cells * 100.
#'
#' Author labels that are NA or not in celltypes count as unlabeled. They do
#' not count as either correct or incorrect in author_recovery_pct.
#' @return data.frame with one row per experiment, preserving input names.
evaluate_annotation_comparison <- function(
    comparison, celltypes,
    prediction_cols = setdiff(names(comparison), c("cell_id", "author"))) {

  .check_eval_comparison(comparison, celltypes, prediction_cols)
  author <- as.character(comparison$author)
  author_labeled <- !is.na(author) & author %in% celltypes
  n_author <- sum(author_labeled)

  rows <- lapply(prediction_cols, function(nm) {
    labels <- as.character(comparison[[nm]])
    ours_labeled <- labels %in% celltypes
    both <- ours_labeled & author_labeled
    correct <- sum(labels[both] == author[both])
    data.frame(
      confident_cells = sum(ours_labeled),
      coverage_pct = 100 * mean(ours_labeled),
      jointly_labeled = sum(both),
      matching_labels = correct,
      agreement_pct = if (sum(both)) 100 * correct / sum(both) else NA_real_,
      author_recovery_pct = if (n_author) 100 * correct / n_author else NA_real_
    )
  })
  result <- do.call(rbind, rows)
  rownames(result) <- prediction_cols
  result
}

#' Two-by-two confidence coverage comparison against author labels.
#' @return Table with rows Labeled, LOWCONF and columns Labeled, Unlabeled.
annotation_coverage_table <- function(comparison, prediction_col, celltypes) {
  .check_eval_comparison(comparison, celltypes, prediction_col)
  author <- as.character(comparison$author)
  ours <- as.character(comparison[[prediction_col]])
  table(
    Ours = factor(
      ifelse(ours %in% celltypes, "Labeled", "LOWCONF"),
      levels = c("Labeled", "LOWCONF")
    ),
    Authors = factor(
      ifelse(!is.na(author) & author %in% celltypes,
             "Labeled", "Unlabeled"),
      levels = c("Labeled", "Unlabeled")
    )
  )
}

#' Per-author-cell-type recovery, for one or multiple reference strategies.
#' @return Long-form data.frame with experiment, cell_type, author_cells,
#'   matching_labels, recovery_pct. Small author groups warrant caution.
annotation_recovery_by_type <- function(
    comparison, celltypes,
    prediction_cols = setdiff(names(comparison), c("cell_id", "author"))) {

  .check_eval_comparison(comparison, celltypes, prediction_cols)
  author <- as.character(comparison$author)
  pieces <- lapply(prediction_cols, function(nm) {
    predictions <- as.character(comparison[[nm]])
    do.call(rbind, lapply(celltypes, function(ct) {
      idx <- which(!is.na(author) & author == ct)
      n <- length(idx)
      correct <- sum(predictions[idx] == ct)
      data.frame(
        experiment = nm,
        cell_type = ct,
        author_cells = n,
        matching_labels = correct,
        recovery_pct = if (n) 100 * correct / n else NA_real_,
        stringsAsFactors = FALSE
      )
    }))
  })
  result <- do.call(rbind, pieces)
  rownames(result) <- NULL
  result
}

#' Compare confident-label overlap between any two experiments.
#' @return List with 2-by-2 status table, count confidently labeled by both,
#'   and label agreement percentage among jointly confident cells.
compare_confident_annotations <- function(
    comparison, first, second, celltypes) {

  .check_eval_comparison(comparison, celltypes, c(first, second))
  a <- as.character(comparison[[first]])
  b <- as.character(comparison[[second]])
  a_ok <- a %in% celltypes
  b_ok <- b %in% celltypes
  both <- a_ok & b_ok
  list(
    coverage_overlap = table(
      First = factor(ifelse(a_ok, "Labeled", "LOWCONF"),
                     levels = c("Labeled", "LOWCONF")),
      Second = factor(ifelse(b_ok, "Labeled", "LOWCONF"),
                      levels = c("Labeled", "LOWCONF"))
    ),
    jointly_confident = sum(both),
    label_agreement_pct = if (any(both)) 100 * mean(a[both] == b[both]) else NA_real_
  )
}
