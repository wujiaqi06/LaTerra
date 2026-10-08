# Literal ADDENDUM002 historical direction, deliberately isolated from the
# general S2 screen. A positive t here retains the historical name n_slow even
# though the numerator is focal-minus-background. Never a generic default.
.lt_marine_fig5b_screen <- function(mat, states) {
  x0 <- mat[, names(states)[states == 0], drop = FALSE]
  x1 <- mat[, names(states)[states == 1], drop = FALSE]
  n0 <- rowSums(!is.na(x0)); n1 <- rowSums(!is.na(x1))
  s0 <- rowSums(x0, na.rm = TRUE); s1 <- rowSums(x1, na.rm = TRUE)
  ss0 <- rowSums(x0 * x0, na.rm = TRUE); ss1 <- rowSums(x1 * x1, na.rm = TRUE)
  m0 <- s0 / n0; m1 <- s1 / n1
  v0 <- (ss0 - s0 * s0 / n0) / (n0 - 1)
  v1 <- (ss1 - s1 * s1 / n1) / (n1 - 1)
  v0[v0 < 0 & v0 > -1e-10] <- 0
  v1[v1 < 0 & v1 > -1e-10] <- 0
  se2 <- v0 / n0 + v1 / n1
  tv <- (m1 - m0) / sqrt(se2)
  df <- se2 * se2 / ((v0 / n0)^2 / (n0 - 1) + (v1 / n1)^2 / (n1 - 1))
  pv <- 2 * stats::pt(abs(tv), df = df, lower.tail = FALSE)
  ok <- n0 > 1 & n1 > 1 & is.finite(tv) & is.finite(pv)
  p <- pv[ok]; tv <- tv[ok]; genes <- rownames(mat)[ok]
  ord <- order(p, na.last = NA)
  selected <- ord[.lt_marine_historical_fdr(p[ord])]
  sig_t <- tv[selected]
  list(n_tested = length(p), n_sig = length(selected),
    n_slow = sum(sig_t > 0, na.rm = TRUE), n_fast = sum(sig_t < 0, na.rm = TRUE),
    slow_proportion = if (length(selected)) sum(sig_t > 0, na.rm = TRUE) / length(selected) else NA_real_,
    sig_genes = genes[selected], positive_t_all = sum(tv > 0, na.rm = TRUE),
    negative_t_all = sum(tv < 0, na.rm = TRUE))
}

.lt_marine_fig5b <- function(gbi, tree, branch_join, drop_traits, worlds, observed_screen,
                             progress = NULL) {
  required <- c("perm_id", "seed", "branch_id", "species_id", "perm_state")
  if (!identical(names(worlds), required) || nrow(worlds) != 53600L ||
      !identical(unique(worlds$perm_id), 1:200) || anyNA(worlds) ||
      any(worlds$seed != 20260520L) || !all(worlds$perm_state %in% 0:1))
    .lt_marine_abort("Fig5B requires the complete original 200-world terminal ledger.", "worlds")
  terminal <- .lt_marine_terminal_frame(branch_join, drop_traits, "drop_whale")
  eligible <- terminal$response != 0.5
  if (sum(!eligible) != 34L || sum(terminal$response == 1) != 17L || sum(terminal$response == 0) != 251L)
    .lt_marine_abort("Fig5B trait domain differs from the admitted drop-whale input.", "worlds")
  rows <- vector("list", 200L)
  branch_worlds <- matrix(NA_real_, nrow = ncol(gbi), ncol = 200L,
    dimnames = list(colnames(gbi), as.character(1:200)))
  for (i in 1:200) {
    world <- worlds[worlds$perm_id == i, , drop = FALSE]
    if (nrow(world) != 268L || anyDuplicated(world$branch_id) ||
        !identical(world$branch_id, terminal$branch[eligible]) ||
        !identical(world$species_id, terminal$species[eligible]) || sum(world$perm_state == 1) != 17L)
      .lt_marine_abort(paste0("Fig5B world identity/count/order mismatch: ", i), "worlds")
    traits <- drop_traits
    traits$drop_whale[match(world$species_id, traits$species)] <- world$perm_state
    bs <- .lt_marine_annotate(tree, traits, branch_join, "drop_whale", "drop_whale")
    branch_worlds[, i] <- bs$state
    states <- stats::setNames(bs$state, bs$branch_id)
    states <- states[states != 0.5 & !is.na(states)]
    r <- .lt_marine_fig5b_screen(gbi, states)
    rows[[i]] <- data.frame(perm_id = i, seed = 20260520L,
      n_positive = sum(states == 1), n_negative = sum(states == 0),
      n_tested_genes = r$n_tested, n_sig_FDR_0_01 = r$n_sig,
      n_slow = r$n_slow, n_fast = r$n_fast, slow_proportion = r$slow_proportion,
      slow_proportion_defined = !is.na(r$slow_proportion),
      positive_t_all = r$positive_t_all, negative_t_all = r$negative_t_all,
      notes = "Matched-positive eligible-terminal-label permutation followed by deterministic ASR; old BH/FDR while-loop rule.",
      stringsAsFactors = FALSE)
    if (!is.null(progress)) progress(i, 200L, rows[[i]])
  }
  summary <- do.call(rbind, rows)
  observed <- observed_screen$significant
  ns <- nrow(observed); slow <- sum(observed$tvalue > 0); fast <- sum(observed$tvalue < 0)
  prop <- if (ns) slow / ns else NA_real_
  defined <- summary$slow_proportion[!is.na(summary$slow_proportion)]
  nd <- length(defined)
  comparison <- data.frame(observed_run_id = "fix_drop_whale", observed_n_sig_FDR_0_01 = ns,
    observed_n_slow = slow, observed_n_fast = fast, observed_slow_proportion = prop,
    n_permutations = 200L, null_median_sig = stats::median(summary$n_sig_FDR_0_01),
    null_max_sig = max(summary$n_sig_FDR_0_01), null_min_sig = min(summary$n_sig_FDR_0_01),
    null_defined_slow_prop_count = nd, null_undefined_slow_prop_count = sum(is.na(summary$slow_proportion)),
    null_median_slow_prop = if (nd) stats::median(defined) else NA_real_,
    null_max_slow_prop = if (nd) max(defined) else NA_real_,
    null_min_slow_prop = if (nd) min(defined) else NA_real_,
    empirical_p_sig_count = (sum(summary$n_sig_FDR_0_01 >= ns) + 1) / 201,
    empirical_p_slow_prop = if (nd) (sum(defined >= prop) + 1) / (nd + 1) else NA_real_,
    notes = "Permutation tests sample-size/null-label behavior only; no biological signal interpretation is made here.",
    stringsAsFactors = FALSE)
  list(screening = summary, observed_vs_null = comparison, branch_worlds = branch_worlds,
    provenance = list(qualifier = "historical_direction_mismatch_preserved",
      terminal_worlds = "original_200_persisted_worlds_no_new_sampling",
      observed_direction = "background_minus_focal", null_direction = "focal_minus_background",
      positive_t_historical_name = "n_slow", correction_authorized = FALSE))
}
