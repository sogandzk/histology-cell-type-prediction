# R/consensus.R
# Reusable consensus aggregation for leave-one-gene-portion-out label transfer.
# No Seurat dependency: takes data.frames returned by run_kfold_transfer().
#
# Defaults reproduce the existing P10 pilot analysis:
# - K runs, one prediction per query cell per run.
# - css_id: most frequent predicted class (ties resolved by celltypes order).
# - nportion.css_id: number of votes for css_id.
# - css_id_full: css_id only if all K runs agree; otherwise LOWCONF.
# - nportion.<celltype>: vote counts for each reference class.
#
# At relaxed thresholds, an exact tie is never considered confident.

.check_consensus_args <- function(celltypes, min_agreement) {
  if (!is.character(celltypes) || !length(celltypes) ||
      anyNA(celltypes) || any(!nzchar(celltypes)) ||
      anyDuplicated(celltypes) || "LOWCONF" %in% celltypes) {
    stop("celltypes must contain unique, nonempty class names, excluding LOWCONF.")
  }
  if (!is.numeric(min_agreement) || length(min_agreement) != 1L ||
      !is.finite(min_agreement) || min_agreement <= 0 ||
      min_agreement > 1) {
    stop("min_agreement must be a single number in (0, 1].")
  }
  invisible(TRUE)
}

#' Calculate a consensus label across K prediction tables.
#'
#' @param predictions_by_portion List of at least two data.frames, each with
#'   cell_id and predicted.id columns (as produced by run_kfold_transfer()).
#' @param celltypes Vector of permitted cell-type names. Order determines
#'   tie-breaking for the informational css_id majority label.
#' @param min_agreement Minimum fraction of matching votes for confidence.
#'   Default 1 means unanimity, matching our authors-reproduction experiments.
#' @param cell_ids Optional query cell IDs, in required output order.
#'   If omitted, use the first prediction table's order. Later tables may have
#'   different row orders, but must contain exactly the same cell IDs.
#' @return data.frame with css_id, nportion.css_id, css_id_full and vote counts.
calculate_consensus <- function(
    predictions_by_portion, celltypes, min_agreement = 1,
    cell_ids = NULL) {

  .check_consensus_args(celltypes, min_agreement)
  if (!is.list(predictions_by_portion) ||
      length(predictions_by_portion) < 2L ||
      any(vapply(predictions_by_portion, is.null, logical(1)))) {
    stop("predictions_by_portion must be a list of at least two prediction tables.")
  }

  if (is.null(cell_ids)) {
    if (!is.data.frame(predictions_by_portion[[1L]]) ||
        !"cell_id" %in% names(predictions_by_portion[[1L]])) {
      stop("First prediction table must contain cell_id.")
    }
    cell_ids <- as.character(predictions_by_portion[[1L]]$cell_id)
  }
  if (!is.character(cell_ids) || !length(cell_ids) ||
      anyNA(cell_ids) || any(!nzchar(cell_ids)) || anyDuplicated(cell_ids)) {
    stop("cell_ids must be unique, nonempty, nonmissing character IDs.")
  }

  k <- length(predictions_by_portion)
  n <- length(cell_ids)
  vote_matrix <- matrix(
    NA_character_, nrow = n, ncol = k,
    dimnames = list(cell_ids, paste0("Run_", seq_len(k)))
  )

  for (i in seq_len(k)) {
    d <- predictions_by_portion[[i]]
    if (!is.data.frame(d) ||
        !all(c("cell_id", "predicted.id") %in% names(d))) {
      stop("Run ", i, " must be a data.frame with cell_id and predicted.id.")
    }
    ids <- as.character(d$cell_id)
    labels <- as.character(d$predicted.id)
    if (length(ids) != n || anyNA(ids) || anyDuplicated(ids) ||
        !setequal(ids, cell_ids)) {
      stop("Run ", i, " does not contain exactly the expected query cells.")
    }
    if (anyNA(labels) || !all(labels %in% celltypes)) {
      stop("Run ", i, " has missing or invalid predicted cell types.")
    }
    vote_matrix[, i] <- labels[match(cell_ids, ids)]
  }

  vote_counts <- matrix(
    0L, nrow = n, ncol = length(celltypes),
    dimnames = list(cell_ids, celltypes)
  )
  for (j in seq_along(celltypes)) {
    vote_counts[, j] <- rowSums(vote_matrix == celltypes[[j]])
  }
  stopifnot(all(rowSums(vote_counts) == k))

  # Like the original script, a tied informational majority uses the first
  # cell type in the provided vector. Ties never count as confident.
  winning_column <- max.col(vote_counts, ties.method = "first")
  max_votes <- as.integer(vote_counts[cbind(seq_len(n), winning_column)])
  tied <- rowSums(vote_counts == max_votes) > 1L
  min_votes <- ceiling(min_agreement * k)
  confident <- (max_votes >= min_votes) & !tied
  majority <- celltypes[winning_column]

  result <- data.frame(
    cell_id = cell_ids,
    css_id = majority,
    nportion.css_id = max_votes,
    css_id_full = ifelse(confident, majority, "LOWCONF"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  vote_df <- as.data.frame(vote_counts, check.names = FALSE)
  names(vote_df) <- paste0("nportion.", celltypes)
  result <- cbind(result, vote_df)
  rownames(result) <- NULL

  validate_consensus(
    result, celltypes = celltypes, n_runs = k,
    min_agreement = min_agreement, expected_cells = cell_ids
  )
  result
}

#' Validate vote counts, ordering, majority label and confidence rule.
#' @return Invisibly TRUE; raises an informative error on failure.
validate_consensus <- function(
    consensus, celltypes, n_runs, min_agreement = 1,
    expected_cells = NULL) {

  .check_consensus_args(celltypes, min_agreement)
  if (!is.numeric(n_runs) || length(n_runs) != 1L ||
      is.na(n_runs) || n_runs < 2 || n_runs != as.integer(n_runs)) {
    stop("n_runs must be an integer >= 2.")
  }
  cols <- paste0("nportion.", celltypes)
  required <- c("cell_id", "css_id", "nportion.css_id", "css_id_full", cols)
  if (!is.data.frame(consensus) || !all(required %in% names(consensus))) {
    stop("Consensus table is missing required output columns.")
  }
  ids <- as.character(consensus$cell_id)
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop("Consensus cell IDs must be unique and nonmissing.")
  }
  if (!is.null(expected_cells) && !identical(ids, expected_cells)) {
    stop("Consensus cell IDs do not match expected_cells in order.")
  }
  counts <- as.matrix(consensus[, cols, drop = FALSE])
  if (!is.numeric(counts) || any(!is.finite(counts)) ||
      any(counts < 0) || any(counts != floor(counts)) ||
      any(rowSums(counts) != n_runs)) {
    stop("Cell-type vote counts must be nonnegative integers totaling n_runs.")
  }
  winner <- max.col(counts, ties.method = "first")
  max_votes <- as.integer(counts[cbind(seq_len(nrow(counts)), winner)])
  tied <- rowSums(counts == max_votes) > 1L
  confident <- (max_votes >= ceiling(min_agreement * n_runs)) & !tied
  expected_majority <- celltypes[winner]
  expected_full <- ifelse(confident, expected_majority, "LOWCONF")

  if (!identical(as.character(consensus$css_id), expected_majority) ||
      !identical(as.integer(consensus$nportion.css_id), max_votes) ||
      !identical(as.character(consensus$css_id_full), expected_full)) {
    stop("Consensus majority or confidence labels do not match the vote counts.")
  }
  invisible(TRUE)
}

#' Summarize consensus coverage and distribution.
#' @return Named list: counts, coverage_pct, labels, max_votes.
summarize_consensus <- function(consensus, celltypes) {
  if (!is.data.frame(consensus) ||
      !all(c("css_id_full", "nportion.css_id") %in% names(consensus))) {
    stop("consensus must contain css_id_full and nportion.css_id.")
  }
  final <- as.character(consensus$css_id_full)
  if (!all(final %in% c(celltypes, "LOWCONF"))) {
    stop("Consensus contains unexpected final labels.")
  }
  n <- length(final)
  if (n == 0L) stop("consensus cannot be empty.")
  n_confident <- sum(final != "LOWCONF")
  list(
    n_cells = n,
    n_confident = n_confident,
    n_lowconf = n - n_confident,
    coverage_pct = 100 * n_confident / n,
    labels = table(factor(final, levels = c(celltypes, "LOWCONF"))),
    max_votes = table(consensus$nportion.css_id)
  )
}
