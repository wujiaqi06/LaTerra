.lt_s2_module_of <- function(genes, lookup) {
  result <- unname(lookup[genes])
  result[is.na(result)] <- "REVIEW_UNANNOTATED"
  result
}

.lt_s2_set_metrics <- function(A, B, lookup) {
  A <- unique(A); B <- unique(B)
  Am <- A[.lt_s2_module_of(A, lookup) != "REVIEW_UNANNOTATED"]
  Bm <- B[.lt_s2_module_of(B, lookup) != "REVIEW_UNANNOTATED"]
  ma <- unique(.lt_s2_module_of(Am, lookup)); mb <- unique(.lt_s2_module_of(Bm, lookup))
  modules <- sort(unique(c(ma, mb)))
  ca <- table(factor(.lt_s2_module_of(Am, lookup), levels = modules))
  cb <- table(factor(.lt_s2_module_of(Bm, lookup), levels = modules))
  cosine <- if (!length(modules) || sqrt(sum(ca^2)) * sqrt(sum(cb^2)) == 0) NA_real_ else
    sum(ca * cb) / (sqrt(sum(ca^2)) * sqrt(sum(cb^2)))
  data.frame(set_A_size_total = length(A), set_B_size_total = length(B),
    set_A_size_module_annotated = length(Am), set_B_size_module_annotated = length(Bm),
    gene_overlap_Jaccard = ifelse(length(union(A, B)) > 0, length(intersect(A, B)) / length(union(A, B)), NA_real_),
    module_presence_Jaccard = ifelse(length(union(ma, mb)) > 0, length(intersect(ma, mb)) / length(union(ma, mb)), NA_real_),
    module_count_cosine_similarity = as.numeric(cosine), stringsAsFactors = FALSE)
}

.lt_marine_turnover_draws <- function(U, Um, nA, nB, nAm, nBm, comparison, lookup) {
  .lt_s2_with_preserved_rng({
    seed <- 20260524L + nchar(comparison) + length(U)
    set.seed(seed)
    rows <- vector("list", 10000L)
    draw_hashes <- character(10000L)
    for (i in seq_len(10000L)) {
      # Exact four-call historical stream, not an equivalent distribution.
      A <- sample(U, nA, replace = FALSE)
      B <- sample(U, nB, replace = FALSE)
      Am <- if (nAm > 0) sample(Um, nAm, replace = FALSE) else character()
      Bm <- if (nBm > 0) sample(Um, nBm, replace = FALSE) else character()
      ma <- unique(.lt_s2_module_of(Am, lookup)); mb <- unique(.lt_s2_module_of(Bm, lookup))
      modules <- sort(unique(c(ma, mb)))
      ca <- table(factor(.lt_s2_module_of(Am, lookup), levels = modules))
      cb <- table(factor(.lt_s2_module_of(Bm, lookup), levels = modules))
      cosine <- if (!length(modules) || sqrt(sum(ca^2)) * sqrt(sum(cb^2)) == 0) NA_real_ else
        sum(ca * cb) / (sqrt(sum(ca^2)) * sqrt(sum(cb^2)))
      rows[[i]] <- data.frame(comparison = comparison, permutation = i,
        gene_overlap_Jaccard = ifelse(length(union(A, B)) > 0, length(intersect(A, B)) / length(union(A, B)), NA_real_),
        module_presence_Jaccard = ifelse(length(union(ma, mb)) > 0, length(intersect(ma, mb)) / length(union(ma, mb)), NA_real_),
        module_count_cosine_similarity = as.numeric(cosine), stringsAsFactors = FALSE)
      draw_hashes[i] <- digest::digest(serialize(list(A = A, B = B, Am = Am, Bm = Bm),
        NULL, version = 3), algo = "sha256", serialize = FALSE)
    }
    list(metrics = do.call(rbind, rows), draw_hashes = draw_hashes, seed = seed,
      final_rng = get(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
  })
}

.lt_marine_turnover <- function(models, full_data, annotation, progress = NULL) {
  lookup <- stats::setNames(annotation$candidate_module, annotation$gene_id)
  .lt_assert_unique_ids(names(lookup), "turnover annotation keys")
  selected <- function(id) {
    beta <- full_data$models[[id]]$beta
    if (is.null(beta)) .lt_marine_abort("Missing this-run full-data predictor set.", "dependency")
    sort(unique(names(beta)[beta != 0]))
  }
  features <- function(id) models[[id]]$global_screen$significant$gene
  slow <- models$fix_drop_whale$global_screen$significant
  comparison <- list(
    list(name = "marine baseline vs whale-only", A_name = "marine baseline", B_name = "whale-only",
      A = selected("fix_marine_binary"), B = selected("fix_whale_only"),
      U = union(features("fix_marine_binary"), features("fix_whale_only"))),
    list(name = "marine baseline vs drop-cetaceans slow genes", A_name = "marine baseline", B_name = "drop-cetaceans slow genes",
      A = selected("fix_marine_binary"), B = slow$gene[slow$tvalue > 0],
      U = union(features("fix_marine_binary"), slow$gene)),
    list(name = "binary aquatic baseline vs aquatic no Cetacea", A_name = "binary aquatic baseline", B_name = "aquatic no Cetacea",
      A = selected("fix_aquatic_v2"), B = selected("fix_aquatic_v2_noCetacea"),
      U = union(features("fix_aquatic_v2"), features("fix_aquatic_v2_noCetacea"))))
  observed <- universes <- worlds <- provenance <- list()
  for (j in seq_along(comparison)) {
    cmp <- comparison[[j]]
    obs <- .lt_s2_set_metrics(cmp$A, cmp$B, lookup)
    U <- unique(cmp$U)
    Um <- U[.lt_s2_module_of(U, lookup) != "REVIEW_UNANNOTATED"]
    universes[[cmp$name]] <- data.frame(comparison = cmp$name, set_A_name = cmp$A_name,
      set_B_name = cmp$B_name,
      universe_definition = "comparison-specific candidate-gene union; module null restricted to module-annotated genes in that union",
      U_gene_size = length(U), U_module_size = length(Um), set_A_total = length(cmp$A),
      set_B_total = length(cmp$B), set_A_module_annotated = obs$set_A_size_module_annotated,
      set_B_module_annotated = obs$set_B_size_module_annotated,
      module_annotation_coverage_A = obs$set_A_size_module_annotated / length(cmp$A),
      module_annotation_coverage_B = obs$set_B_size_module_annotated / length(cmp$B),
      stringsAsFactors = FALSE)
    draws <- .lt_marine_turnover_draws(U, Um, length(cmp$A), length(cmp$B),
      obs$set_A_size_module_annotated, obs$set_B_size_module_annotated, cmp$name, lookup)
    worlds[[cmp$name]] <- draws$metrics
    observed[[cmp$name]] <- cbind(data.frame(comparison = cmp$name, set_A_name = cmp$A_name,
      set_B_name = cmp$B_name, stringsAsFactors = FALSE), obs)
    provenance[[cmp$name]] <- list(ordered_U_gene = U, ordered_U_module = Um,
      ordered_A = cmp$A, ordered_B = cmp$B, seed = draws$seed,
      ordered_draw_sha256 = draws$draw_hashes, final_rng = draws$final_rng,
      source_rule = "20260524 + nchar(comparison) + length(U_gene); sequential A/B/Am/Bm")
    if (!is.null(progress)) progress(j, length(comparison), observed[[cmp$name]])
  }
  turnover <- do.call(rbind, observed)
  perms <- do.call(rbind, worlds)
  summaries <- lapply(split(perms, perms$comparison), function(df) {
    obs <- turnover[turnover$comparison == unique(df$comparison), , drop = FALSE]
    do.call(rbind, lapply(c("gene_overlap_Jaccard", "module_presence_Jaccard", "module_count_cosine_similarity"), function(metric) {
      vals <- df[[metric]]; value <- obs[[metric]]
      data.frame(comparison = unique(df$comparison), metric = metric, observed_value = value,
        null_median = stats::median(vals, na.rm = TRUE),
        null_95_interval = paste(sprintf("%.4f", stats::quantile(vals, c(0.025, 0.975), na.rm = TRUE)), collapse = "-"),
        empirical_p_value = (sum(vals >= value, na.rm = TRUE) + 1) / (sum(is.finite(vals)) + 1),
        z_score = ifelse(stats::sd(vals, na.rm = TRUE) > 0,
          (value - mean(vals, na.rm = TRUE)) / stats::sd(vals, na.rm = TRUE), NA_real_),
        enrichment_ratio = value / stats::median(vals, na.rm = TRUE),
        n_permutations = sum(is.finite(vals)), stringsAsFactors = FALSE)
    }))
  })
  list(observed = turnover, universes = do.call(rbind, universes),
    worlds = perms, summary = do.call(rbind, summaries), provenance = provenance)
}
