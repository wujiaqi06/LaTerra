# Deterministic terminal reporting from fresh S2 outputs. These helpers neither
# screen genes nor compute a baseline, trait, fit, random world, or new null.
.lt_marine_table_s1 <- function(traits) {
  fields <- c("species", "aquaticity_score_sum_0_18", "marine_binary", "aquatic_v1_status",
    "aquatic_v2_status", "aquatic_v3_status", "legacy_marine2", "legacy_aquatic2",
    "aquatic_trait_source", "aquatic_scoring_scheme", "aquatic_score_range")
  if (!is.data.frame(traits) || !all(fields %in% names(traits)))
    .lt_marine_abort("Missing consumed trait fields for S1 reporting.", "reporting_fields")
  .lt_assert_unique_ids(traits$species, "S1 consumed species")
  if (anyNA(traits[c("marine_binary", "aquatic_v2_status")]) ||
      any(!traits$marine_binary %in% c(0, 1)) || any(!traits$aquatic_v2_status %in% c(0, .5, 1)))
    .lt_marine_abort("Invalid existing trait coding; S1 does not rescore traits.", "reporting_payload")
  out <- traits[fields]; names(out)[1L] <- "species_id"
  out$endpointfix_trait_category <- ifelse(out$marine_binary == 1, "marine",
    ifelse(out$aquatic_v2_status == 1, "non_marine_aquatic",
      ifelse(out$aquatic_v2_status == .5, "semi_aquatic", "terrestrial")))
  out
}

.lt_marine_table_bh <- function(p) {
  if (!is.numeric(p) || anyNA(p) || any(!is.finite(p)) || any(p < 0 | p > 1))
    .lt_marine_abort("Displayed S2 BH requires the complete finite tested p-values.", "reporting_payload")
  if (!length(p)) return(numeric())
  # Literal May22 source order: ranked_p * n / (i+1), reverse cumulative min.
  idx <- order(p, seq_along(p), method = "radix")
  q <- rev(cummin(rev(p[idx] * length(p) / seq_along(p))))
  out <- numeric(length(p)); out[idx] <- pmin(q, 1)
  out
}

.lt_marine_table_s2 <- function(tested, significant, trait) {
  if (length(trait) != 1L || is.na(trait) || !trait %in% c("marine_binary", "aquatic_v2"))
    .lt_marine_abort("Unknown historical S2 reporting trait.", "reporting_keys")
  required <- c("gene", "tvalue", "pvalue", "mean1", "mean2", "marine", "non_marine")
  if (!all(required %in% names(tested)) || !all(required %in% names(significant)))
    .lt_marine_abort("Incomplete this-run screen output for S2.", "reporting_fields")
  .lt_assert_unique_ids(tested$gene, "S2 tested genes")
  .lt_assert_unique_ids(significant$gene, "S2 already-selected genes")
  idx <- match(significant$gene, tested$gene)
  if (anyNA(idx) || !isTRUE(all.equal(significant[required], tested[idx, required],
      check.attributes = FALSE, tolerance = 0)))
    .lt_marine_abort("S2 significant rows must be the unchanged selected subset of this-run tested rows.", "reporting_keys")
  q <- .lt_marine_table_bh(tested$pvalue)
  direction <- ifelse(significant$tvalue > 0, "slow", "fast")
  # The new q-values are presentation only: never reselect the strict prefix.
  data.frame(gene_id = significant$gene, trait = trait, run_id = paste0("fix_", trait),
    t_value = significant$tvalue, p_value = significant$pvalue, FDR = q[idx],
    direction = direction, slow_or_fast = direction,
    mean_trait_state_1 = significant$mean2, mean_trait_state_0 = significant$mean1,
    n_state_1 = significant$marine, n_state_0 = significant$non_marine,
    stringsAsFactors = FALSE)
}

.lt_marine_table_overlap <- function(marine_tested, aquatic_tested, marine_selected, aquatic_selected) {
  # Final source ac80fb64..., lines 217-281. Counts are built before looking
  # up fixed metric labels, and no stored count/null result is used below.
  for (x in list(marine_tested, aquatic_tested, marine_selected, aquatic_selected))
    .lt_assert_unique_ids(x$gene, "baseline overlap gene keys")
  if (!all(marine_selected$gene %in% marine_tested$gene) ||
      !all(aquatic_selected$gene %in% aquatic_tested$gene))
    .lt_marine_abort("Overlap selected sets are not within their current tested universes.", "reporting_keys")
  ms <- marine_selected$gene[marine_selected$tvalue > 0]
  mf <- marine_selected$gene[marine_selected$tvalue <= 0]
  as <- aquatic_selected$gene[aquatic_selected$tvalue > 0]
  af <- aquatic_selected$gene[aquatic_selected$tvalue <= 0]
  sig_m <- marine_selected$gene; sig_a <- aquatic_selected$gene
  n <- function(x) length(x)
  both <- intersect(ms, as); union_slow <- union(ms, as)
  both_sig <- intersect(sig_m, sig_a); union_sig <- union(sig_m, sig_a)
  cross1 <- intersect(ms, af); cross2 <- intersect(mf, as); fast_both <- intersect(mf, af)
  concordant <- n(both) + n(fast_both); discordant <- n(cross1) + n(cross2)
  fmt <- function(x, digits) as.numeric(sprintf(paste0("%.", digits, "f"), x))
  pct <- function(a, b) {
    if (b == 0L) .lt_marine_abort("Historical overlap percentage has a zero denominator.", "reporting_payload")
    fmt(100 * a / b, 2L)
  }
  ratio <- function(a, b) {
    if (b == 0L) .lt_marine_abort("Historical overlap Jaccard has a zero denominator.", "reporting_payload")
    fmt(a / b, 4L)
  }
  metric <- c("marine_tested_genes", "aquatic_tested_genes", "marine_FDR_significant_genes",
    "marine_slow_genes", "marine_fast_genes", "aquatic_FDR_significant_genes", "aquatic_slow_genes",
    "aquatic_fast_genes", "shared_slow", "marine_only_slow", "aquatic_only_slow", "all_slow_union",
    "marine_slow_shared_with_aquatic_slow_percent", "aquatic_slow_shared_with_marine_slow_percent",
    "slow_jaccard_index", "marine_slow_intersect_aquatic_fast", "marine_fast_intersect_aquatic_slow",
    "marine_fast_intersect_aquatic_fast", "marine_significant_intersect_aquatic_significant",
    "significant_gene_jaccard_index", "direction_concordant_among_both_significant",
    "direction_discordant_among_both_significant")
  value <- c(nrow(marine_tested), nrow(aquatic_tested), n(sig_m), n(ms), n(mf), n(sig_a), n(as), n(af),
    n(both), n(setdiff(ms, as)), n(setdiff(as, ms)), n(union_slow), pct(n(both), n(ms)),
    pct(n(both), n(as)), ratio(n(both), n(union_slow)), n(cross1), n(cross2), n(fast_both),
    n(both_sig), ratio(n(both_sig), n(union_sig)), concordant, discordant)
  denominator <- rep(NA_real_, length(metric))
  denominator[c(3, 4, 5, 6, 7, 8, 13, 14, 21, 22)] <-
    c(nrow(marine_tested), n(sig_m), n(sig_m), nrow(aquatic_tested), n(sig_a), n(sig_a),
      n(ms), n(as), n(both_sig), n(both_sig))
  percent <- rep(NA_real_, length(metric))
  for (i in c(3, 4, 5, 6, 7, 8, 21, 22)) percent[i] <- pct(value[i], denominator[i])
  # Explicit alias join between the recovered count-table metric IDs and the
  # ADDENDUM007 final display schema; no positional gene/result selection.
  labels <- c("marine tested genes", "aquatic tested genes", "marine FDR-significant genes",
    "marine slow genes", "marine fast genes", "aquatic FDR-significant genes", "aquatic slow genes",
    "aquatic fast genes", "shared slow genes", "marine-only slow genes", "aquatic-only slow genes",
    "all slow union", "percent of marine slow genes shared with aquatic slow",
    "percent of aquatic slow genes shared with marine slow", "slow-gene Jaccard index",
    "marine slow \u2229 aquatic fast", "marine fast \u2229 aquatic slow", "marine fast \u2229 aquatic fast",
    "all significant gene overlap regardless of direction", "significant-gene Jaccard index",
    "direction concordant among genes significant in both screens",
    "direction discordant among genes significant in both screens")
  alias <- stats::setNames(labels, metric)
  out <- data.frame(metric = unname(alias[metric]), value = value, denominator = denominator, percent_or_rate = percent,
    stringsAsFactors = FALSE)
  layout <- .lt_marine_reporting_layout("OUT_TABLE_S4", "Sheet1")
  keys <- layout$schema$ordered_keys[1:22]
  original_metric <- vapply(keys, `[[`, "", 2L)
  if (!setequal(out$metric, original_metric))
    .lt_marine_abort("Computed overlap metric coverage differs from the original schema.", "reporting_keys")
  group <- vapply(keys, `[[`, "", 1L)
  out$metric_group <- group[match(out$metric, original_metric)]
  # Check the context of the separately imported six-row historical reference.
  # Its historical null outcomes remain untouched and are never recalculated.
  historical_metrics <- vapply(layout$schema$ordered_keys[23:28], `[[`, "", 2L)
  context_rows <- match(c("null universe size", "observed shared slow genes"), historical_metrics) + 23L
  observed_context <- c(n(intersect(marine_tested$gene, aquatic_tested$gene)), n(both))
  archived_context <- vapply(context_rows, function(r) as.numeric(layout$values[[r, 3L]]), 0)
  if (anyNA(context_rows) || !identical(as.numeric(observed_context), archived_context))
    .lt_marine_abort("Fresh baseline universe/overlap conflicts with the context of the fixed historical reference; neither is repaired.",
      "reporting_historical_context")
  attr(out, "historical_context_check") <- list(
    fresh_tested_intersection = n(intersect(marine_tested$gene, aquatic_tested$gene)),
    fresh_slow_overlap = n(both), historical_metric_keys = historical_metrics,
    source = "this_run_baseline_sets_not_null_outcomes")
  out
}
