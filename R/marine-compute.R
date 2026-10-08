# Pure operators scoped to the frozen Marine/S2 recipe. No C3 or oracle reads.

.lt_marine_effect <- function(values, margin) {
  fun <- function(x) {
    x <- x[is.finite(x)]
    if (length(x)) c(mean = mean(x), median = stats::median(x)) else
      c(mean = NA_real_, median = NA_real_)
  }
  t(apply(values, margin, fun))
}

.lt_marine_baseline <- function(raw) {
  values <- raw$values
  cutoffs <- vapply(seq_len(nrow(values)), function(i) {
    x <- values[i, ]
    if (all(is.na(x))) return(NA_real_)
    as.numeric(stats::quantile(x, probs = 0.975, type = 7,
                              na.rm = TRUE, names = FALSE))
  }, numeric(1L))
  trimmed <- values
  mask <- matrix(FALSE, nrow(values), ncol(values), dimnames = dimnames(values))
  for (i in seq_len(nrow(values))) {
    selected <- is.finite(values[i, ]) & !is.na(cutoffs[i]) &
      values[i, ] >= cutoffs[i]
    mask[i, selected] <- TRUE
    trimmed[i, selected] <- NA_real_
  }
  ge <- .lt_marine_effect(trimmed, 1L)
  be <- .lt_marine_effect(trimmed, 2L)
  # Preserve the frozen source's floating-point operation order.
  gbi <- sweep(sweep(trimmed, 2L, be[, "mean"], "/"), 1L, ge[, "mean"], "/")
  gbi[!is.finite(gbi)] <- NA_real_
  effect_table <- function(margin, effect) {
    count <- function(x, f) if (margin == 1L) rowSums(f(x)) else colSums(f(x))
    ids <- if (margin == 1L) raw$gene_ids else raw$branch_ids
    tab <- data.frame(id = ids, median = effect[, "median"], mean = effect[, "mean"],
      n_present_before_trim = count(values, is.finite),
      n_present_after_trim = count(trimmed, is.finite),
      n_trimmed = if (margin == 1L) rowSums(mask) else colSums(mask),
      zero_count_before_trim = count(values, function(x) is.finite(x) & x == 0),
      zero_count_after_trim = count(trimmed, function(x) is.finite(x) & x == 0),
      row.names = NULL, check.names = FALSE)
    names(tab)[1:3] <- if (margin == 1L) c("gene", "GE_median", "GE_mean") else
      c("branch", "BE_median", "BE_mean")
    tab$ratio_median_over_mean <- tab[[2L]] / tab[[3L]]
    tab
  }
  list(gbi = gbi, trimmed_values = trimmed, trim_mask = mask,
       trim_cutoff = stats::setNames(cutoffs, raw$gene_ids),
       gene_effect = effect_table(1L, ge), branch_effect = effect_table(2L, be))
}

.lt_marine_branch_join <- function(tree, old_keys, branch_ids) {
  .lt_assert_unique_ids(old_keys$branch_label, "old branch labels")
  .lt_assert_unique_ids(old_keys$canonical_split_key, "supplied scientific split keys")
  if (!identical(old_keys$branch_label, branch_ids)) {
    .lt_marine_abort("Old split-key ledger and matrix branch order differ.", "axis")
  }
  children <- split(tree$edge[, 2L], tree$edge[, 1L])
  ntip <- length(tree$tip.label)
  descendants <- function(node) {
    if (node <= ntip) return(tree$tip.label[node])
    unlist(lapply(children[[as.character(node)]], descendants), use.names = FALSE)
  }
  # Join tree annotation edges to existing split sides. No new coordinate truth,
  # cross-gene matching, reclassification or branch-label conversion is performed.
  side_key <- function(x) paste(sort(x, method = "radix"), collapse = ";")
  normalize_side <- function(x) vapply(strsplit(x, ";", fixed = TRUE), side_key, "")
  a <- normalize_side(old_keys$side_A_taxa)
  b <- normalize_side(old_keys$side_B_taxa)
  edge_side <- vapply(tree$edge[, 2L], function(node) side_key(descendants(node)), "")
  matches <- lapply(edge_side, function(s) which(a == s | b == s))
  if (any(lengths(matches) != 1L)) {
    .lt_marine_abort("Tree edge does not join uniquely to a supplied old split key.", "branch_join")
  }
  index <- unlist(matches, use.names = FALSE)
  if (anyDuplicated(index) || length(index) != length(branch_ids)) {
    .lt_marine_abort("Tree/ledger edge join is not one-to-one and exhaustive.", "branch_join")
  }
  edge_order <- match(seq_along(branch_ids), index)
  data.frame(branch_id = branch_ids,
    canonical_split_key = old_keys$canonical_split_key,
    ancestor = tree$edge[edge_order, 1L], offspring = tree$edge[edge_order, 2L],
    tree_edge_index = edge_order, branch_type = old_keys$branch_type,
    terminal_taxon = old_keys$terminal_taxon_if_terminal, stringsAsFactors = FALSE)
}

.lt_marine_asr <- function(tree, y) {
  # castor's compiled entry initializes .Random.seed via its Rcpp wrapper even
  # for this deterministic routine. Preserve caller state, including absence.
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (had_seed) assign(".Random.seed", seed, envir = .GlobalEnv) else
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  fit <- castor::asr_max_parsimony(tree, y + 1, length(unique(y)))
  if (had_seed && !identical(seed, get0(".Random.seed", envir = .GlobalEnv,
                                       inherits = FALSE))) {
    .lt_marine_abort("Deterministic ASR unexpectedly consumed the caller RNG stream.", "asr_rng")
  }
  fit
}

.lt_marine_annotate <- function(tree, traits, branch_join, column, trait_id) {
  .lt_assert_unique_ids(traits$species, "trait species")
  .lt_assert_unique_ids(tree$tip.label, "tree tips")
  if (!setequal(traits$species, tree$tip.label) || !column %in% names(traits)) {
    .lt_marine_abort("Trait/tree species or frozen trait column mismatch.", "trait")
  }
  y <- traits[[column]][match(tree$tip.label, traits$species)]
  allowed <- if (trait_id == "marine_binary") c(0, 1) else c(0, 0.5, 1)
  if (!is.numeric(y) || anyNA(y) || any(!y %in% allowed)) {
    .lt_marine_abort("Trait contains a state outside its frozen vocabulary.", "trait")
  }
  # Deliberately mirror the historical call, including fractional tip coding.
  # Internal output uses column_index - 1, not a newly designed state remapping.
  fit <- .lt_marine_asr(tree, y)
  likelihood <- fit$ancestral_likelihoods
  if (is.null(likelihood) || anyNA(likelihood) || any(!is.finite(likelihood))) {
    .lt_marine_abort("Frozen ASR failed; no alternative reconstruction is used.", "asr")
  }
  state <- c(y, max.col(likelihood, ties.method = "first") - 1)
  branch_state <- state[branch_join$offspring]
  if (anyNA(branch_state) || any(!branch_state %in% allowed)) {
    .lt_marine_abort("ASR generated an unsupported frozen branch state.", "asr")
  }
  data.frame(run_id = paste0("fix_", trait_id), trait = trait_id,
    branch_id = branch_join$branch_id, state = branch_state,
    screened = branch_state %in% c(0, 1),
    reason = ifelse(branch_state == 0.5, "intermediate_excluded", NA_character_),
    ancestor = branch_join$ancestor, offspring = branch_join$offspring,
    stringsAsFactors = FALSE)
}

.lt_marine_historical_fdr <- function(sorted_p) {
  n <- length(sorted_p)
  if (!n) return(integer())
  pass <- sorted_p < seq_len(n) / n * 0.01
  failed <- which(!pass | is.na(pass))
  end <- if (length(failed)) failed[1L] - 1L else n
  seq_len(end)
}

.lt_marine_screen <- function(gbi, branch_states) {
  if (!identical(colnames(gbi), branch_states$branch_id)) {
    .lt_marine_abort("Screen branch states and GBI axes differ.", "axis")
  }
  y <- branch_states$state
  groups <- list(which(y == 0), which(y == 1))
  rows <- vector("list", nrow(gbi)); reasons <- rep(NA_character_, nrow(gbi))
  nf <- nr <- integer(nrow(gbi))
  for (i in seq_len(nrow(gbi))) {
    ref <- gbi[i, groups[[1L]]]; ref <- ref[is.finite(ref)]
    focal <- gbi[i, groups[[2L]]]; focal <- focal[is.finite(focal)]
    nr[i] <- length(ref); nf[i] <- length(focal)
    if (min(nr[i], nf[i]) < 2L) {
      reasons[i] <- "insufficient_finite_group_size"
      next
    }
    ans <- tryCatch(stats::t.test(ref, focal, alternative = "two.sided",
                                 var.equal = FALSE), error = identity)
    if (inherits(ans, "error")) {
      .lt_marine_abort(paste0("Welch failed for gene ", rownames(gbi)[i],
        ": ", conditionMessage(ans), ". No silent skip or alternate test."), "welch")
    }
    if (!is.finite(ans$p.value) || !is.finite(ans$statistic)) {
      .lt_marine_abort(paste0("Non-finite Welch result for ", rownames(gbi)[i]), "welch")
    }
    rows[[i]] <- data.frame(gene = rownames(gbi)[i],
      tvalue = unname(ans$statistic), pvalue = ans$p.value,
      mean1 = unname(ans$estimate[1L]), mean2 = unname(ans$estimate[2L]),
      marine = nf[i], non_marine = nr[i], stringsAsFactors = FALSE)
  }
  tested <- do.call(rbind, rows)
  if (is.null(tested)) tested <- data.frame(gene = character(), tvalue = numeric(),
    pvalue = numeric(), mean1 = numeric(), mean2 = numeric(), marine = integer(),
    non_marine = integer())
  tested <- tested[order(tested$pvalue), , drop = FALSE]; rownames(tested) <- NULL
  selected <- .lt_marine_historical_fdr(tested$pvalue)
  significant <- tested[selected, , drop = FALSE]
  ledger <- data.frame(gene = rownames(gbi), gene_order = seq_len(nrow(gbi)),
    n_focal = nf, n_reference = nr, tested = is.na(reasons), reason = reasons,
    FDR_selected = rownames(gbi) %in% significant$gene,
    stringsAsFactors = FALSE)
  ti <- match(ledger$gene, tested$gene)
  ledger$beta_background_minus_focal <- tested$mean1[ti] - tested$mean2[ti]
  ledger$direction <- ifelse(is.na(ti), NA_character_,
    ifelse(tested$tvalue[ti] > 0, "slow", ifelse(tested$tvalue[ti] < 0, "fast", "zero")))
  summary <- data.frame(run_id = unique(branch_states$run_id),
    trait_name = unique(branch_states$trait), n_branches_input = ncol(gbi),
    n_branches_screened = sum(branch_states$screened), n_tested_genes = nrow(tested),
    n_sig_fdr01 = nrow(significant), positive_t_all = sum(tested$tvalue > 0),
    negative_t_all = sum(tested$tvalue < 0), positive_t_sig = sum(significant$tvalue > 0),
    negative_t_sig = sum(significant$tvalue < 0), stringsAsFactors = FALSE)
  list(tested = tested, significant = significant, gene_ledger = ledger, summary = summary)
}
