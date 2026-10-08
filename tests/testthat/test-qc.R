qc_matrix_fixture <- function(values, coord_state = NULL,
                              gene_ids = paste0("g", seq_len(nrow(values))),
                              branch_ids = paste0("b", seq_len(ncol(values)))) {
  values <- as.matrix(values); storage.mode(values) <- "double"
  dimnames(values) <- list(gene_ids, branch_ids)
  if (is.null(coord_state)) {
    coord_state <- matrix("observed", nrow(values), ncol(values),
                          dimnames = dimnames(values))
  } else dimnames(coord_state) <- dimnames(values)
  reason <- matrix(NA_character_, nrow(values), ncol(values),
                   dimnames = dimnames(values))
  reason[coord_state == "observed" & is.na(values)] <- "current_payload_unavailable"
  lt_matrix(
    values, coord_state, fixture_payload(), value_reason = reason,
    coordinate_provenance = fixture_coordinate_provenance()
  )
}

qc_rule <- function(method, ..., rule_id = paste0("test_", method)) {
  list(list(rule_id = rule_id, method = method, ...))
}

test_that("QC A diagnose-only leaves every raw matrix field identical", {
  x <- small_lt_matrix()
  before <- unserialize(serialize(x, NULL, version = 3))
  qc <- lt_qc_matrix(x)
  expect_identical(x$values, before$values)
  expect_identical(x$coord_state, before$coord_state)
  expect_identical(x$value_reason, before$value_reason)
  expect_identical(x$gene_ids, before$gene_ids)
  expect_identical(x$branch_ids, before$branch_ids)
  expect_false(qc$diagnostic_provenance$raw_matrix_modified)
})

test_that("QC B diagnostic flags do not alter eligibility", {
  x <- qc_matrix_fixture(matrix(c(1, 1, 1, 1, 100), 1))
  before <- lt_measurement_eligibility(x)
  qc <- lt_qc_matrix(x, lt_qc_spec(
    profile = "custom", outlier_rules = qc_rule(
      "iqr", multiplier = 1.5, scope = "whole_matrix", direction = "upper"
    )
  ))
  after <- lt_measurement_eligibility(x)
  expect_gt(nrow(qc$flags), 0)
  expect_identical(before$status, after$status)
  expect_identical(before$reason, after$reason)
  expect_false(qc$diagnostic_provenance$eligibility_modified)
})

test_that("QC C explicit apply marks selected observed cells ineligible", {
  x <- qc_matrix_fixture(matrix(c(1, 1, 1, 1, 100), 1))
  qc <- lt_qc_matrix(x, lt_qc_spec(
    profile = "custom", outlier_rules = qc_rule(
      "iqr", multiplier = 1.5, scope = "whole_matrix", direction = "upper"
    )
  ))
  out <- lt_apply_qc(
    x, qc, list(flag_ids = qc$flags$flag_id),
    decision_id = "explicit_extreme_exclusion", authority = "unit_test"
  )
  expect_identical(out$status[1, 5], "ineligible")
  expect_identical(out$reason[1, 5], "measurement_qc_excluded")
  expect_true(out$provenance$qc_application$explicit_application)
})

test_that("QC D non-observed cells remain scoped not_applicable after apply", {
  x <- small_lt_matrix()
  qc <- lt_qc_matrix(x)
  mask <- matrix(FALSE, nrow(x$values), ncol(x$values), dimnames = dimnames(x$values))
  mask[1, 1] <- TRUE
  out <- lt_apply_qc(x, qc, list(mask = mask),
                     decision_id = "observed_only", authority = "unit_test")
  expect_true(all(out$status[x$coord_state != "observed"] == "not_applicable"))
  expect_true(all(out$reason[x$coord_state != "observed"] == "coordinate_not_observed"))
  mask[2, 1] <- TRUE
  expect_error(lt_apply_qc(x, qc, list(mask = mask),
    decision_id = "invalid", authority = "unit_test"), "non-observed")
})

test_that("QC E exact zero survives unless explicitly selected independently", {
  x <- qc_matrix_fixture(matrix(c(0, 1, 1, 1, 100), 1))
  qc <- lt_qc_matrix(x, lt_qc_spec(
    profile = "custom", outlier_rules = qc_rule(
      "iqr", multiplier = 1.5, scope = "whole_matrix", direction = "upper"
    )
  ))
  flagged <- lt_apply_qc(x, qc, list(flag_ids = qc$flags$flag_id),
    decision_id = "exclude_flags", authority = "unit_test")
  expect_identical(flagged$status[1, 1], "eligible")
  expect_identical(x$values[1, 1], 0)
  zero_mask <- matrix(FALSE, 1, 5, dimnames = dimnames(x$values)); zero_mask[1, 1] <- TRUE
  explicit <- lt_apply_qc(x, qc, list(mask = zero_mask),
    decision_id = "independent_zero_decision", authority = "unit_test")
  expect_identical(explicit$status[1, 1], "ineligible")
  expect_identical(x$values[1, 1], 0)
})

test_that("QC F extreme finite value is warning-flagged but eligible by default", {
  x <- qc_matrix_fixture(matrix(c(1, 1, 1, 1, 1e150), 1))
  qc <- lt_qc_matrix(x, lt_qc_spec(
    profile = "custom", outlier_rules = qc_rule(
      "iqr", multiplier = 1.5, scope = "whole_matrix", direction = "upper"
    )
  ))
  eligibility <- lt_measurement_eligibility(x)
  expect_true(any(qc$flags$branch_id == "b5"))
  expect_true(all(qc$flags$severity == "WARNING"))
  expect_identical(eligibility$status[1, 5], "eligible")
  expect_equal(x$values[1, 5], 1e150)
})

test_that("QC G exact point mass is reported without inferred upstream cause", {
  x <- qc_matrix_fixture(matrix(c(2, 2, 2, 3, 4), 1))
  qc <- lt_qc_matrix(x, lt_qc_spec(outlier_rules = list()))
  point <- qc$matrix_summary$point_masses
  expect_true(any(point$value == 2 & point$count == 3))
  expect_true(all(grepl("no_upstream_cause_inferred", point$interpretation, fixed = TRUE)))
})

test_that("QC H quantile flags are keyed-permutation invariant", {
  y <- matrix(c(1, 2, 100, 3, 4, 200), 2, 3, byrow = TRUE)
  x <- qc_matrix_fixture(y, gene_ids = c("ga", "gb"), branch_ids = c("ba", "bb", "bc"))
  spec <- lt_qc_spec(profile = "custom", outlier_rules = qc_rule(
    "quantile", scope = "whole_matrix", direction = "upper",
    lower_probability = 0, upper_probability = .7, tie_semantics = "strict"
  ))
  a <- lt_qc_matrix(x, spec)
  xp <- qc_matrix_fixture(y[c(2, 1), c(3, 1, 2)],
    gene_ids = c("gb", "ga"), branch_ids = c("bc", "ba", "bb"))
  b <- lt_qc_matrix(xp, spec)
  key <- function(z) z$flags[order(z$flags$gene_id, z$flags$branch_id),
    c("flag_id", "gene_id", "branch_id", "diagnostic_value", "threshold_rule")]
  rownames_a <- key(a); rownames_b <- key(b)
  rownames(rownames_a) <- NULL; rownames(rownames_b) <- NULL
  expect_identical(rownames_a, rownames_b)
  expect_identical(a$output_hashes$flag_ledger_sha256,
                   b$output_hashes$flag_ledger_sha256)
})

test_that("QC I zero-MAD is undefined safely without Inf or NaN", {
  x <- qc_matrix_fixture(matrix(c(1, 1, 1, 100), 1))
  spec <- lt_qc_spec(profile = "custom", outlier_rules = qc_rule(
    "mad", scope = "whole_matrix", direction = "two_sided", multiplier = 3
  ))
  qc <- lt_qc_matrix(x, spec)
  expect_equal(nrow(qc$flags), 0)
  expect_true(any(qc$warnings$warning_code == "MAD_ZERO_UNDEFINED"))
  expect_false(any(is.infinite(qc$flags$robust_score)))
  expect_false(any(is.nan(qc$flags$robust_score)))
})

test_that("QC J signed diagnostics name absolute concentration explicitly", {
  x <- qc_matrix_fixture(matrix(c(-4, -1, 0, 2, 3), 1))
  qc <- lt_qc_matrix(x)
  expect_identical(qc$matrix_summary$overall$resolved_value_domain, "signed")
  expect_identical(qc$matrix_summary$overall$concentration_measure, "absolute_value")
  expect_true(all(grepl("absolute", qc$gene_summary$concentration_definition)))
  expect_equal(qc$matrix_summary$overall$negative_count, 2)
})

test_that("QC K positive-only metrics record zero/negative domain loss", {
  x <- qc_matrix_fixture(matrix(c(-1, 0, 1, 2), 1))
  qc <- lt_qc_matrix(x)
  domain <- qc$matrix_summary$metric_domains
  positive <- domain[domain$metric == "positive_log_distribution", ]
  expect_equal(positive$admitted_count, 2)
  expect_equal(positive$domain_loss_count, 2)
  expect_identical(positive$zero_handling, "outside_metric_domain_no_pseudocount")
  expect_identical(x$values, matrix(c(-1, 0, 1, 2), 1,
    dimnames = list("g1", paste0("b", 1:4))))
})

test_that("QC L custom user mask and diagnostic round-trip", {
  x <- qc_matrix_fixture(matrix(1:4, 2, 2))
  mask <- matrix(c(TRUE, FALSE, FALSE, TRUE), 2, 2, dimnames = dimnames(x$values))
  spec <- lt_qc_spec(profile = "custom", outlier_rules = qc_rule(
    "custom_mask", mask = mask, scope = "whole_matrix", direction = "two_sided"
  ))
  qc <- lt_qc_matrix(x, spec)
  expect_equal(nrow(qc$flags), 2)
  path_spec <- tempfile(fileext = ".rds"); path_qc <- tempfile(fileext = ".rds")
  lt_write_object(spec, path_spec); lt_write_object(qc, path_qc)
  expect_identical(lt_read_object(path_spec), spec)
  expect_identical(lt_read_object(path_qc), qc)
})

test_that("QC M threshold changes flags but not raw matrix identity", {
  x <- qc_matrix_fixture(matrix(c(1, 1, 1, 1, 100), 1))
  a <- lt_qc_matrix(x, lt_qc_spec(profile = "custom", outlier_rules = qc_rule(
    "iqr", multiplier = 1.5, scope = "whole_matrix", direction = "upper"
  )))
  b <- lt_qc_matrix(x, lt_qc_spec(profile = "custom", outlier_rules = qc_rule(
    "iqr", multiplier = 200, scope = "whole_matrix", direction = "upper"
  )))
  expect_false(identical(a$diagnostic_spec_id, b$diagnostic_spec_id))
  expect_false(identical(a$output_hashes$flag_ledger_sha256,
                         b$output_hashes$flag_ledger_sha256))
  expect_identical(a$matrix_identity, b$matrix_identity)
})

test_that("QC N same matrix and spec replay identically", {
  x <- qc_matrix_fixture(matrix(c(0, 1, 1, 2, 100, -3), 2, 3))
  spec <- lt_qc_spec()
  a <- lt_qc_matrix(x, spec); b <- lt_qc_matrix(x, spec)
  expect_identical(a$output_hashes, b$output_hashes)
  expect_identical(a$flags, b$flags)
  expect_identical(a$warnings, b$warnings)
})

test_that("QC O warning cannot mutate supplied measurement eligibility", {
  x <- qc_matrix_fixture(matrix(c(1, 1, 1, 1, 100), 1))
  e <- lt_measurement_eligibility(x)
  before <- unserialize(serialize(e, NULL, version = 3))
  qc <- lt_qc_matrix(x)
  expect_gt(nrow(qc$warnings), 0)
  expect_identical(e, before)
  expect_identical(e$status_sha256, before$status_sha256)
})

test_that("QC P HOLD is limited to established declared-domain failure", {
  x <- qc_matrix_fixture(matrix(c(-1, 1, 100), 1))
  ordinary <- lt_qc_matrix(x)
  expect_equal(nrow(ordinary$holds), 0)
  held <- lt_qc_matrix(x, lt_qc_spec(value_domain = "nonnegative"))
  expect_equal(nrow(held$holds), 1)
  expect_identical(held$holds$hold_code,
                   "DECLARED_NONNEGATIVE_DOMAIN_VIOLATION")
  expect_error(lt_apply_qc(x, held, list(flag_ids = character()),
    decision_id = "held", authority = "unit_test"), "held")
})

test_that("QC Q applied Layer-2 object satisfies T1A-FIX001 scope", {
  x <- small_lt_matrix()
  qc <- lt_qc_matrix(x)
  mask <- matrix(FALSE, nrow(x$values), ncol(x$values), dimnames = dimnames(x$values))
  mask[1, 1] <- TRUE
  out <- lt_apply_qc(x, qc, list(mask = mask),
    decision_id = "scoped", authority = "unit_test")
  expect_invisible(validate_lt_measurement_eligibility(out, x))
  expect_true(all(out$status[x$coord_state != "observed"] == "not_applicable"))
  expect_false(any(out$status[x$coord_state == "observed"] == "not_applicable"))
})

test_that("QC R standard diagnose-only profile is threshold-free", {
  x <- qc_matrix_fixture(matrix(c(0, 1, 2, 100), 1))
  spec <- lt_qc_spec(profile = "standard")
  qc <- lt_qc_matrix(x, spec)
  expect_length(spec$outlier_rules, 0)
  expect_equal(nrow(qc$threshold_config), 0)
  expect_equal(nrow(qc$flags), 0)
  expect_gt(nrow(qc$matrix_summary$tail_concentration), 0)
  expect_gt(nrow(qc$matrix_summary$leverage), 0)
  expect_true(qc$diagnostic_provenance$standard_profile_threshold_free)
})

test_that("QC S policy distinctions are frozen in the spec", {
  spec <- lt_qc_spec()
  expect_identical(spec$semantic_contract, c(
    "outlier_not_error", "long_branch_not_anomalous_cell", "trim_not_diagnosis"
  ))
  expect_identical(spec$standard_profile_contract, "threshold_free_diagnose_only")
  expect_identical(spec$threshold_authority, "explicit_optional_rule_not_universal_truth")
})

test_that("QC T explicit quantile tie semantics distinguish strict and inclusive", {
  x <- qc_matrix_fixture(matrix(c(1, 2, 2, 2), 1))
  make_spec <- function(ties) lt_qc_spec(
    profile = "custom", outlier_rules = qc_rule(
      "quantile", scope = "whole_matrix", direction = "upper",
      lower_probability = 0, upper_probability = .75, tie_semantics = ties
    )
  )
  strict <- lt_qc_matrix(x, make_spec("strict"))
  inclusive <- lt_qc_matrix(x, make_spec("inclusive"))
  expect_equal(nrow(strict$flags), 0)
  expect_equal(nrow(inclusive$flags), 3)
  expect_identical(strict$threshold_config$tie_semantics, "strict")
  expect_identical(inclusive$threshold_config$tie_semantics, "inclusive")
  expect_true(all(grepl("tie_semantics=inclusive", inclusive$flags$threshold_rule,
                        fixed = TRUE)))
})

test_that("QC U quantile exact-rank semantics are parked, not silently approximated", {
  expect_error(lt_qc_spec(
    profile = "custom", outlier_rules = qc_rule(
      "quantile", scope = "within_gene", direction = "upper",
      lower_probability = 0, upper_probability = .975,
      tie_semantics = "exact_rank"
    )
  ), "parked and not implemented")
})

test_that("QC V two-axis empirical extremeness is a diagnostic-only placeholder", {
  x <- qc_matrix_fixture(matrix(c(1, 2, 3, 4), 2, 2))
  qc <- lt_qc_matrix(x)
  future <- qc$cell_summary$future_diagnostics$two_axis_empirical_extremeness
  expect_identical(future$status, "PARKED_NOT_IMPLEMENTED")
  expect_true(future$diagnostic_only)
  expect_identical(future$fields$field,
                   c("within_gene_percentile", "within_branch_percentile"))
  expect_true(all(future$fields$status == "not_computed"))
  expect_null(future$threshold)
  expect_identical(future$admission_effect, "none")
})

test_that("QC W standard profile rejects embedded optional threshold rules", {
  expect_error(lt_qc_spec(
    profile = "standard", outlier_rules = qc_rule(
      "iqr", scope = "whole_matrix", direction = "upper", multiplier = 1.5
    )
  ), "standard.*threshold-free")
})

test_that("QC FIX002 A top-N point-mass ranking remains functional", {
  x <- qc_matrix_fixture(matrix(c(0, 0, rep(5, 5), 10, 10), 1))
  qc <- lt_qc_matrix(x, lt_qc_spec(point_mass_top_n = 1L))
  expect_equal(nrow(qc$matrix_summary$point_masses), 1L)
  expect_equal(qc$matrix_summary$point_masses$value, 5)
  expect_equal(qc$matrix_summary$point_masses$count, 5L)
})

test_that("QC FIX002 B global maximum is retained outside top-N frequency", {
  x <- qc_matrix_fixture(matrix(c(0, 0, rep(5, 5), 10, 10), 1))
  qc <- lt_qc_matrix(x, lt_qc_spec(point_mass_top_n = 1L))
  extreme <- qc$matrix_summary$extreme_exact_repetition
  expect_false(any(qc$matrix_summary$point_masses$value == 10))
  expect_equal(extreme$maximum_finite_value, 10)
  expect_equal(extreme$maximum_exact_count, 2L)
  expect_equal(extreme$maximum_exact_fraction_available, 2 / 9)
})

test_that("QC FIX002 C global minimum is retained outside top-N frequency", {
  x <- qc_matrix_fixture(matrix(c(0, 0, rep(5, 5), 10, 10), 1))
  qc <- lt_qc_matrix(x, lt_qc_spec(point_mass_top_n = 1L))
  extreme <- qc$matrix_summary$extreme_exact_repetition
  expect_false(any(qc$matrix_summary$point_masses$value == 0))
  expect_equal(extreme$minimum_finite_value, 0)
  expect_equal(extreme$minimum_exact_count, 2L)
  expect_equal(extreme$minimum_exact_fraction_available, 2 / 9)
})

test_that("QC FIX002 D extreme repetition makes no upstream-cause inference", {
  x <- qc_matrix_fixture(matrix(c(1, 1, 2, 3, 3), 1))
  extreme <- lt_qc_matrix(x)$matrix_summary$extreme_exact_repetition
  expect_identical(
    extreme$interpretation,
    "exact_extreme_repetition_no_upstream_cause_inferred"
  )
})

test_that("QC FIX002 E marine maximum target can record 50 exactly 274 times", {
  x <- qc_matrix_fixture(matrix(c(rep(1, 300), rep(50, 274)), 1))
  extreme <- lt_qc_matrix(
    x, lt_qc_spec(point_mass_top_n = 1L)
  )$matrix_summary$extreme_exact_repetition
  expect_equal(extreme$maximum_finite_value, 50)
  expect_equal(extreme$maximum_exact_count, 274L)
})

test_that("QC FIX002 F selected eligible becomes measurement-QC excluded", {
  x <- qc_matrix_fixture(matrix(1, 1, 1))
  qc <- lt_qc_matrix(x)
  mask <- matrix(TRUE, 1, 1, dimnames = dimnames(x$values))
  out <- lt_apply_qc(x, qc, list(mask = mask),
    decision_id = "eligible_transition", authority = "unit_test")
  expect_identical(out$status[1, 1], "ineligible")
  expect_identical(out$reason[1, 1], "measurement_qc_excluded")
  expect_equal(out$provenance$qc_application$newly_excluded_count, 1L)
})

test_that("QC FIX002 G selected already-ineligible preserves exact reason", {
  x <- qc_matrix_fixture(matrix(1, 1, 1))
  base <- lt_measurement_eligibility(
    x, matrix("ineligible", 1, 1, dimnames = dimnames(x$values)),
    matrix("user_predeclared_excluded", 1, 1, dimnames = dimnames(x$values))
  )
  qc <- lt_qc_matrix(x)
  mask <- matrix(TRUE, 1, 1, dimnames = dimnames(x$values))
  out <- lt_apply_qc(x, qc, list(mask = mask), base,
    decision_id = "preserve_upstream", authority = "unit_test")
  expect_identical(out$status[1, 1], "ineligible")
  expect_identical(out$reason[1, 1], "user_predeclared_excluded")
  expect_equal(out$provenance$qc_application$already_ineligible_selected_count, 1L)
})

test_that("QC FIX002 H selected unresolved records prior state when resolved", {
  x <- qc_matrix_fixture(matrix(1, 1, 1))
  base <- lt_measurement_eligibility(
    x, matrix("unresolved", 1, 1, dimnames = dimnames(x$values)),
    matrix("authority_unresolved", 1, 1, dimnames = dimnames(x$values))
  )
  qc <- lt_qc_matrix(x)
  mask <- matrix(TRUE, 1, 1, dimnames = dimnames(x$values))
  out <- lt_apply_qc(x, qc, list(mask = mask), base,
    decision_id = "resolve_unresolved", authority = "unit_test")
  application <- out$provenance$qc_application
  expect_identical(out$status[1, 1], "ineligible")
  expect_identical(out$reason[1, 1], "measurement_qc_excluded")
  expect_equal(application$unresolved_resolved_to_ineligible_count, 1L)
  expect_identical(application$selected_cell_prior_state$prior_status, "unresolved")
  expect_identical(application$selected_cell_prior_state$prior_reason,
                   "authority_unresolved")
})

test_that("QC FIX002 I selecting not-applicable fails closed", {
  x <- small_lt_matrix()
  qc <- lt_qc_matrix(x)
  mask <- matrix(FALSE, nrow(x$values), ncol(x$values),
                 dimnames = dimnames(x$values))
  mask[2, 1] <- TRUE
  expect_error(lt_apply_qc(x, qc, list(mask = mask),
    decision_id = "forbidden_not_applicable", authority = "unit_test"),
    "not_applicable.*non-observed")
})

test_that("QC FIX002 J application provenance transition counts are exact", {
  x <- qc_matrix_fixture(matrix(1:3, 1, 3))
  base <- lt_measurement_eligibility(
    x,
    matrix(c("eligible", "ineligible", "unresolved"), 1, 3,
           dimnames = dimnames(x$values)),
    matrix(c("eligible_declared", "analysis_domain_excluded", "authority_unresolved"),
           1, 3, dimnames = dimnames(x$values)),
    provenance = list(authority = "upstream_fixture")
  )
  qc <- lt_qc_matrix(x)
  mask <- matrix(TRUE, 1, 3, dimnames = dimnames(x$values))
  out <- lt_apply_qc(x, qc, list(mask = mask), base,
    decision_id = "three_way_partition", authority = "unit_test")
  application <- out$provenance$qc_application
  expect_identical(unname(unlist(application[c(
    "selected_cell_count", "newly_excluded_count",
    "already_ineligible_selected_count",
    "unresolved_resolved_to_ineligible_count"
  )])), c(3L, 1L, 1L, 1L))
  expect_equal(nrow(application$selected_cell_prior_state), 3L)
  expect_identical(application$prior_status_sha256, base$status_sha256)
  expect_identical(application$prior_reason_sha256, base$reason_sha256)
  expect_length(application$selected_cell_keys, 3L)
})

test_that("QC FIX002 K diagnose-only raw matrix identity remains unchanged", {
  x <- small_lt_matrix()
  before <- unserialize(serialize(x, NULL, version = 3))
  qc <- lt_qc_matrix(x)
  expect_identical(x, before)
  expect_identical(qc$matrix_identity$combined_sha256,
                   lt_qc_matrix(before)$matrix_identity$combined_sha256)
})

test_that("QC FIX002 L diagnose-only Layer-2 identity remains unchanged", {
  x <- qc_matrix_fixture(matrix(1:3, 1, 3))
  base <- lt_measurement_eligibility(
    x,
    matrix(c("eligible", "ineligible", "unresolved"), 1, 3,
           dimnames = dimnames(x$values)),
    matrix(c("eligible_declared", "user_predeclared_excluded", "authority_unresolved"),
           1, 3, dimnames = dimnames(x$values))
  )
  before <- unserialize(serialize(base, NULL, version = 3))
  invisible(lt_qc_matrix(x))
  expect_identical(base, before)
  expect_identical(base$status_sha256, before$status_sha256)
  expect_identical(base$reason_sha256, before$reason_sha256)
})
