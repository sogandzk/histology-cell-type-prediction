# R/gene_portions.R
# Shared-gene selection, CELLxGENE gene-statistic filtering, and
# reproducible held-out gene portions for Xenium label transfer.
#
# Input expression is NOT read by these functions. All gene statistics come
# from rowData(sce), so cell sampling does not change the gene selection.
#
# Important: The optional mean-stratified generator below is a NEW, explicit
# algorithm. It is not claimed to reproduce the original authors' or the
# earlier pilot's exact portion assignments. To reuse an existing experiment
# exactly, supply its original gene_portions to make_gene_portions(existing=).

#' Find shared Xenium/reference genes and filter by reference gene statistics.
#'
#' Keeps the order of gene symbols in the Xenium assay, as in
#' intersect(rownames(xmo[["Xenium"]]), ref_genes).
#'
#' @return A list with shared_genes, eligible_genes, excluded_genes,
#'   gene_stats (all shared genes and statistics), and counts.
prepare_transfer_genes <- function(
    sce, query, query_assay = "Xenium",
    gene_col = "gene_symbols",
    mean_col = "vf_vst_counts_mean",
    variance_col = "vf_vst_counts_variance.standardized") {

  if (!inherits(sce, "SingleCellExperiment")) {
    stop("sce must be a SingleCellExperiment.")
  }
  if (!inherits(query, "Seurat")) stop("query must be a Seurat object.")
  if (!query_assay %in% names(query@assays)) {
    stop("Query assay is missing: ", query_assay)
  }

  rd <- SummarizedExperiment::rowData(sce)
  required <- c(gene_col, mean_col, variance_col)
  missing_cols <- setdiff(required, colnames(rd))
  if (length(missing_cols)) {
    stop("Reference rowData is missing: ", paste(missing_cols, collapse = ", "))
  }

  ref_symbols <- as.character(rd[[gene_col]])
  query_symbols <- rownames(query[[query_assay]])
  query_symbols <- query_symbols[!is.na(query_symbols) & nzchar(query_symbols)]

  if (anyDuplicated(query_symbols)) {
    stop("Query assay has duplicated gene symbols.")
  }

  shared <- intersect(query_symbols, ref_symbols)
  if (!length(shared)) stop("No shared gene symbols found.")

  # Reject ambiguity in the full reference metadata; match() alone would
  # silently choose the first of multiple rows with an identical symbol.
  shared_ref_symbols <- ref_symbols[!is.na(ref_symbols) & ref_symbols %in% shared]
  if (anyDuplicated(shared_ref_symbols)) {
    stop("Reference contains duplicated shared gene symbols.")
  }

  ref_idx <- match(shared, ref_symbols)
  if (anyNA(ref_idx)) stop("Failed to match some shared gene symbols.")

  means <- suppressWarnings(as.numeric(rd[[mean_col]][ref_idx]))
  variances <- suppressWarnings(as.numeric(rd[[variance_col]][ref_idx]))

  # The legacy pilot excluded genes with zero mean or zero standardized
  # variance; treat missing/nonfinite or negative statistics as ineligible.
  eligible <- is.finite(means) & is.finite(variances) &
    means > 0 & variances > 0

  stats <- data.frame(
    gene = shared,
    ref_row_index = as.integer(ref_idx),
    mean = means,
    variance_standardized = variances,
    eligible = eligible,
    stringsAsFactors = FALSE
  )

  list(
    shared_genes = shared,
    eligible_genes = shared[eligible],
    excluded_genes = shared[!eligible],
    gene_stats = stats,
    counts = c(
      reference_genes = length(ref_symbols),
      xenium_genes = length(query_symbols),
      shared_genes = length(shared),
      excluded_genes = sum(!eligible),
      eligible_genes = sum(eligible)
    )
  )
}

#' Verify that gene portions form a complete, nonoverlapping partition.
#' @return Invisibly, a table containing the number of genes per portion.
validate_gene_portions <- function(portions, eligible_genes, k = length(portions)) {
  if (!is.list(portions) || !length(portions) ||
      length(k) != 1L || !is.finite(k) || k != as.integer(k) ||
      k < 2L || length(portions) != k) {
    stop("portions must be a list with exactly k >= 2 components.")
  }
  if (!identical(names(portions), as.character(seq_len(k)))) {
    stop("Portions must be named '1', '2', ..., as.character(k), in order.")
  }
  if (!is.character(eligible_genes) || !length(eligible_genes) ||
      anyNA(eligible_genes) || anyDuplicated(eligible_genes)) {
    stop("eligible_genes must be a unique, nonmissing character vector.")
  }
  valid <- vapply(portions, function(x) {
    is.character(x) && length(x) > 0 && !anyNA(x) &&
      all(nzchar(x)) && !anyDuplicated(x)
  }, logical(1))
  if (!all(valid)) stop("Each portion must contain distinct, valid gene symbols.")

  assigned <- unlist(portions, use.names = FALSE)
  if (anyDuplicated(assigned)) stop("Some genes occur in more than one portion.")
  if (!setequal(assigned, eligible_genes)) {
    stop("Portions do not contain exactly the eligible genes.")
  }
  invisible(data.frame(
    portion = seq_len(k),
    n_genes = as.integer(lengths(portions))
  ))
}

#' Return training genes (all but held_out) and test genes (held_out).
get_transfer_gene_split <- function(portions, held_out, eligible_genes = NULL) {
  k <- length(portions)
  if (!is.null(eligible_genes)) {
    validate_gene_portions(portions, eligible_genes, k)
  }
  if (length(held_out) != 1L || is.na(held_out) ||
      !is.numeric(held_out) || held_out != floor(held_out) ||
      held_out < 1L || held_out > k) {
    stop("held_out must be one integer between 1 and length(portions).")
  }
  train <- unlist(portions[setdiff(as.character(seq_len(k)),
                                   as.character(held_out))], use.names = FALSE)
  test <- portions[[as.character(held_out)]]
  if (length(intersect(train, test))) stop("Training and test genes overlap.")
  list(genes_train = train, genes_test = test)
}

#' Make or reuse K gene portions.
#'
#' IMPORTANT: For exact reproduction of the current pilot's experiments,
#' pass existing = gene_portions. The generator is a NEW mean-stratified
#' option and does not promise the exact same bin assignments.
#'
#' @param gene_info Output of prepare_transfer_genes().
#' @param existing Optional previously built and verified portion list.
#' @param method "mean_stratified" or "random" (for NEW experiments only).
#' @param block_size Number of adjacent mean-ranked genes per block, used by
#'   mean_stratified allocation. Each block is shuffled before round-robin
#'   allocation, approximately preserving the mean distribution in each bin.
make_gene_portions <- function(
    gene_info, k = 10L, seed = 42L,
    method = c("mean_stratified", "random"),
    block_size = NULL, existing = NULL) {

  if (!is.list(gene_info) || is.null(gene_info$eligible_genes) ||
      is.null(gene_info$gene_stats)) {
    stop("gene_info must be an output of prepare_transfer_genes().")
  }
  eligible <- gene_info$eligible_genes
  if (length(k) != 1L || is.na(k) || !is.numeric(k) ||
      k != floor(k) || k < 2L || k > length(eligible)) {
    stop("k must be an integer between 2 and the eligible gene count.")
  }
  k <- as.integer(k)

  if (!is.null(existing)) {
    validate_gene_portions(existing, eligible, k)
    return(existing)
  }

  method <- match.arg(method)
  if (length(seed) != 1L || is.na(seed) || !is.numeric(seed)) {
    stop("seed must be a single numeric value.")
  }
  if (is.null(block_size)) block_size <- 5L * k
  if (length(block_size) != 1L || is.na(block_size) ||
      !is.numeric(block_size) || !is.finite(block_size) ||
      block_size != floor(block_size) || block_size < k) {
    stop("block_size must be a positive integer >= k.")
  }

  set.seed(seed)
  portions <- setNames(replicate(k, character(0), simplify = FALSE),
                       as.character(seq_len(k)))

  if (method == "random") {
    ordered_genes <- sample(eligible, length(eligible))
  } else {
    stats <- gene_info$gene_stats
    pos <- match(eligible, stats$gene)
    if (anyNA(pos) || any(!is.finite(stats$mean[pos]))) {
      stop("Missing finite mean statistics for some eligible genes.")
    }
    ordered_genes <- eligible[order(stats$mean[pos], eligible)]
  }

  # Round-robin allocation keeps total portion sizes within one gene of
  # each other. Mean-stratified mode randomizes within nearby mean ranks.
  blocks <- split(seq_along(ordered_genes),
                  ceiling(seq_along(ordered_genes) / block_size))
  offset <- 0L
  for (ids in blocks) {
    genes <- ordered_genes[ids]
    if (method == "mean_stratified") genes <- sample(genes, length(genes))
    bucket <- ((seq_along(genes) - 1L + offset) %% k) + 1L
    for (j in seq_len(k)) {
      portions[[j]] <- c(portions[[j]], genes[bucket == j])
    }
    offset <- offset + length(genes)
  }
  validate_gene_portions(portions, eligible, k)
  portions
}
