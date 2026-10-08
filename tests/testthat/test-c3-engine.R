c3_fixture <- function(values, mask = matrix(TRUE, nrow(values), ncol(values)),
                       gene_ids = paste0("g", seq_len(nrow(values))),
                       branch_ids = paste0("b", seq_len(ncol(values)))) {
  stopifnot(identical(dim(values), dim(mask)))
  values <- as.matrix(values)
  storage.mode(values) <- "double"
  values[!mask] <- NA_real_
  dimnames(values) <- list(gene_ids, branch_ids)
  state <- matrix("observed", nrow(values), ncol(values), dimnames = dimnames(values))
  state[!mask] <- "NA_struct"
  lt_matrix(
    values = values, coord_state = state,
    payload = fixture_payload(),
    coordinate_provenance = fixture_coordinate_provenance()
  )
}

keyed_mu <- function(fit) {
  out <- data.frame(
    gene = fit$domain$matrix$gene_ids[fit$edge_gene_index],
    branch = fit$domain$matrix$branch_ids[fit$edge_branch_index],
    mu = fit$mu, stringsAsFactors = FALSE
  ) |>
    transform(key = paste(gene, branch, sep = "\r")) |>
    (function(x) x[order(x$key), c("key", "mu")])()
  rownames(out) <- NULL
  out
}

layer2_coordinate_fixture <- function(coord_state) {
  value <- if (identical(coord_state, "observed")) 1 else NA_real_
  values <- matrix(value, 1, 1, dimnames = list("g1", "b1"))
  state <- matrix(coord_state, 1, 1, dimnames = dimnames(values))
  lt_matrix(
    values, state, fixture_payload(),
    coordinate_provenance = fixture_coordinate_provenance()
  )
}

test_that("A complete positive table has the closed-form C3 solution", {
  y <- matrix(c(1, 2, 3, 4, 5, 6), 2, 3, byrow = TRUE)
  fit <- lt_c3_fit(lt_c3_domain(c3_fixture(y)), run_id = "A")
  expected <- outer(rowSums(y), colSums(y)) / sum(y)
  expect_identical(fit$state, "FINITE_INTERIOR")
  expect_equal(fit$mu, expected[fit$edge_linear], tolerance = 1e-11)
  expect_true(all(vapply(fit$certificates, function(x) {
    validate_lt_c3_certificate(x); TRUE
  }, logical(1))))
})

test_that("B admitted exact zeros retain positive fitted means when interior", {
  y <- matrix(c(1, 0, 0, 1), 2, 2, byrow = TRUE)
  domain <- lt_c3_domain(c3_fixture(y))
  fit <- lt_c3_fit(domain, run_id = "B")
  expect_equal(sum(domain$edge_value == 0), 2)
  expect_identical(fit$state, "FINITE_INTERIOR")
  expect_equal(fit$mu, rep(0.5, 4), tolerance = 1e-12)
  expect_true(all(fit$mu > 0))
})

test_that("C zero-margin gene is an explicit boundary with zero means", {
  y <- matrix(c(1, 2, 0, 0), 2, 2, byrow = TRUE)
  fit <- lt_c3_fit(lt_c3_domain(c3_fixture(y)), run_id = "C")
  expect_identical(fit$support$zero_margin_vertices$vertex_id, "g2")
  ii <- fit$edge_gene_index == 2L
  expect_true(all(fit$baseline_support_state[ii] == "ZERO_MARGIN_BOUNDARY"))
  expect_true(all(fit$mu[ii] == 0))
})

test_that("D zero-margin branch is an explicit boundary with zero means", {
  y <- matrix(c(1, 0, 2, 0), 2, 2, byrow = TRUE)
  fit <- lt_c3_fit(lt_c3_domain(c3_fixture(y)), run_id = "D")
  expect_identical(fit$support$zero_margin_vertices$vertex_id, "b2")
  ii <- fit$edge_branch_index == 2L
  expect_true(all(fit$baseline_support_state[ii] == "ZERO_MARGIN_BOUNDARY"))
  expect_true(all(fit$mu[ii] == 0))
})

test_that("E zero-total component bypasses the interior solver", {
  fit <- lt_c3_fit(lt_c3_domain(c3_fixture(matrix(0, 2, 2))), run_id = "E")
  expect_identical(fit$state, "ZERO_TOTAL_COMPONENT")
  expect_length(fit$component_fits, 0)
  expect_true(all(fit$mu == 0))
  expect_true(all(fit$baseline_support_state == "ZERO_TOTAL_COMPONENT"))
})

test_that("F unsupported gene receives no factor or baseline edge", {
  y <- matrix(c(1, 2, 3, 4), 2, 2, byrow = TRUE)
  mask <- matrix(c(TRUE, TRUE, FALSE, FALSE), 2, 2, byrow = TRUE)
  fit <- lt_c3_fit(lt_c3_domain(c3_fixture(y, mask)), run_id = "F")
  expect_true(any(fit$support$unsupported_vertices$vertex_id == "g2"))
  expect_false(any(fit$edge_gene_index == 2L))
  expect_false(any(vapply(fit$component_fits, function(z) 2L %in% z$gene_index, logical(1))))
})

test_that("G unsupported branch receives no factor or baseline edge", {
  y <- matrix(c(1, 2, 3, 4), 2, 2, byrow = TRUE)
  mask <- matrix(c(TRUE, FALSE, TRUE, FALSE), 2, 2, byrow = TRUE)
  fit <- lt_c3_fit(lt_c3_domain(c3_fixture(y, mask)), run_id = "G")
  expect_true(any(fit$support$unsupported_vertices$vertex_id == "b2"))
  expect_false(any(fit$edge_branch_index == 2L))
})

test_that("H disconnected positive components receive separate gauges", {
  y <- matrix(c(2, 0, 0, 3), 2, 2, byrow = TRUE)
  mask <- matrix(c(TRUE, FALSE, FALSE, TRUE), 2, 2, byrow = TRUE)
  fit <- lt_c3_fit(lt_c3_domain(c3_fixture(y, mask)), run_id = "H")
  expect_equal(nrow(fit$support$final_components$summary), 2)
  expect_length(fit$certificates, 2)
  expect_equal(fit$mu, c(2, 3))
})

test_that("I connected incomplete mask with positive support is interior", {
  y <- matrix(c(1, 2, 0, 3), 2, 2, byrow = TRUE)
  mask <- matrix(c(TRUE, TRUE, FALSE, TRUE), 2, 2, byrow = TRUE)
  fit <- lt_c3_fit(lt_c3_domain(c3_fixture(y, mask)), run_id = "I")
  expect_equal(nrow(fit$support$initial_components$summary), 1)
  expect_identical(fit$state, "FINITE_INTERIOR")
  expect_true(all(fit$mu > 0))
})

test_that("J connected proper face returns minimal forced-zero support", {
  y <- matrix(c(1, 0, 0, 1), 2, 2, byrow = TRUE)
  mask <- matrix(c(TRUE, TRUE, FALSE, TRUE), 2, 2, byrow = TRUE)
  fit <- lt_c3_fit(lt_c3_domain(c3_fixture(y, mask)), run_id = "J")
  expect_identical(fit$state, "FACIAL_BOUNDARY")
  expect_equal(sum(fit$support$facial_forced_edge_mask), 1)
  forced <- fit$baseline_support_state == "FACIAL_BOUNDARY"
  expect_equal(sum(forced), 1)
  expect_equal(fit$mu[forced], 0)
  face <- fit$support$existence_certificates[[1]]$certificate
  expect_true(face$minimal_face_certified)
  expect_identical(face$dual_certificate_type, "EXACT_INTEGER_CONDENSATION_EXPOSER")
})

test_that("K connectivity alone never upgrades a proper face to interior", {
  y <- matrix(c(1, 0, 0, 1), 2, 2, byrow = TRUE)
  mask <- matrix(c(TRUE, TRUE, FALSE, TRUE), 2, 2, byrow = TRUE)
  support <- lt_c3_support(lt_c3_domain(c3_fixture(y, mask)))
  expect_equal(nrow(support$initial_components$summary), 1)
  expect_identical(support$state, "FACIAL_BOUNDARY")
})

test_that("L deterministic certificate resource ceiling fails closed", {
  y <- matrix(c(1, 0, 0, 1), 2, 2, byrow = TRUE)
  mask <- matrix(c(TRUE, TRUE, FALSE, TRUE), 2, 2, byrow = TRUE)
  fit <- lt_c3_fit(
    lt_c3_domain(c3_fixture(y, mask)),
    lt_c3_control(certificate_max_edges = 2), run_id = "L"
  )
  expect_identical(fit$state, "NUMERICALLY_INDETERMINATE")
  expect_null(fit$mu)
})

test_that("M global scaling scales mu and preserves support", {
  y <- matrix(c(1, 0, 2, 3), 2, 2, byrow = TRUE)
  base <- lt_c3_fit(lt_c3_domain(c3_fixture(y)), run_id = "M1")
  for (scale in c(0.1, 10)) {
    scaled <- lt_c3_fit(lt_c3_domain(c3_fixture(y * scale)), run_id = "M2")
    expect_equal(scaled$mu, base$mu * scale, tolerance = 1e-10)
    expect_identical(scaled$baseline_support_state, base$baseline_support_state)
  }
})

test_that("N gene and branch permutations preserve keyed fitted means", {
  y <- matrix(c(1, 0, 2, 3, 4, 5), 2, 3, byrow = TRUE)
  base <- lt_c3_fit(lt_c3_domain(c3_fixture(y)), run_id = "N1")
  perm <- lt_c3_fit(lt_c3_domain(c3_fixture(
    y[c(2, 1), c(3, 1, 2)], gene_ids = c("g2", "g1"),
    branch_ids = c("b3", "b1", "b2")
  )), run_id = "N2")
  expect_equal(keyed_mu(base), keyed_mu(perm), tolerance = 1e-10)
})

test_that("O deterministic replay produces identical hashes", {
  fit <- lt_c3_fit(lt_c3_domain(c3_fixture(matrix(c(1, 2, 3, 4), 2, 2))),
                   run_id = "O")
  expect_true(all(vapply(fit$certificates, `[[`, logical(1),
                         "deterministic_replay_pass")))
  expect_true(all(vapply(fit$certificates, function(z) {
    identical(z$mu_sha256, z$replay_mu_sha256)
  }, logical(1))))
})

test_that("P off-domain products are never materialized", {
  y <- matrix(c(1, 2, 3, 4), 2, 2, byrow = TRUE)
  mask <- matrix(c(TRUE, TRUE, FALSE, TRUE), 2, 2, byrow = TRUE)
  fit <- lt_c3_fit(lt_c3_domain(c3_fixture(y, mask)), run_id = "P")
  expect_length(fit$mu, sum(mask))
  expect_false(any(fit$edge_gene_index == 2L & fit$edge_branch_index == 1L))
})

test_that("Q diagnostic warnings cannot mutate eligibility", {
  x <- c3_fixture(matrix(c(1, 0, 2, 3), 2, 2))
  e <- lt_measurement_eligibility(x)
  a <- lt_c3_domain(x, e)
  b <- lt_c3_domain(x, e, diagnostic_warnings = c("extreme finite value"))
  expect_identical(a$eligibility$status_sha256, b$eligibility$status_sha256)
  expect_identical(a$fit_mask, b$fit_mask)
  expect_identical(a$edge_value, b$edge_value)
})

test_that("R finite extreme magnitudes remain admitted numeric inputs", {
  y <- matrix(c(1e150, 2e150, 3e150, 4e150), 2, 2, byrow = TRUE)
  domain <- lt_c3_domain(c3_fixture(y))
  fit <- lt_c3_fit(domain, run_id = "R")
  expect_equal(length(domain$edge_linear), 4)
  expect_identical(fit$state, "FINITE_INTERIOR")
  expect_true(all(is.finite(fit$mu)))
})

test_that("layer-2 and payload contract violations fail explicitly", {
  x <- c3_fixture(matrix(c(1, 2, 3, 4), 2, 2))
  bad <- matrix("eligible", 1, 2)
  expect_error(lt_measurement_eligibility(x, bad, bad), "dimensions")
  status <- matrix("eligible", 2, 2, dimnames = dimnames(x$values))
  reason <- matrix("eligible_declared", 2, 2, dimnames = dimnames(x$values))
  status[1, 1] <- "maybe"
  expect_error(lt_measurement_eligibility(x, status, reason), "Unknown")

  negative <- c3_fixture(matrix(c(-1, 2, 3, 4), 2, 2))
  status[1, 1] <- "eligible"
  eligibility <- lt_measurement_eligibility(negative, status, reason)
  expect_error(lt_c3_domain(negative, eligibility), "negative")
})

test_that("C3 objects round-trip with all support layers", {
  fit <- lt_c3_fit(lt_c3_domain(c3_fixture(matrix(c(1, 0, 0, 1), 2, 2))),
                   run_id = "round_trip")
  path <- tempfile(fileext = ".rds")
  lt_write_object(fit, path)
  expect_identical(lt_read_object(path), fit)
})

test_that("FIX001 A NA_struct is Layer-2 not_applicable", {
  e <- lt_measurement_eligibility(layer2_coordinate_fixture("NA_struct"))
  expect_identical(e$status[[1]], "not_applicable")
  expect_identical(e$reason[[1]], "coordinate_not_observed")
})

test_that("FIX001 B NA_fuse is Layer-2 not_applicable", {
  e <- lt_measurement_eligibility(layer2_coordinate_fixture("NA_fuse"))
  expect_identical(e$status[[1]], "not_applicable")
  expect_identical(e$reason[[1]], "coordinate_not_observed")
})

test_that("FIX001 C NA_topo is Layer-2 not_applicable", {
  e <- lt_measurement_eligibility(layer2_coordinate_fixture("NA_topo"))
  expect_identical(e$status[[1]], "not_applicable")
  expect_identical(e$reason[[1]], "coordinate_not_observed")
})

test_that("FIX001 D residual_NA is Layer-2 not_applicable", {
  e <- lt_measurement_eligibility(layer2_coordinate_fixture("residual_NA"))
  expect_identical(e$status[[1]], "not_applicable")
  expect_identical(e$reason[[1]], "coordinate_not_observed")
})

test_that("FIX001 E non-observed coordinate cannot be explicitly eligible", {
  x <- layer2_coordinate_fixture("NA_struct")
  expect_error(lt_measurement_eligibility(
    x, matrix("eligible", 1, 1), matrix("eligible_declared", 1, 1)
  ), "Non-observed")
})

test_that("FIX001 F non-observed coordinate cannot be explicitly ineligible", {
  x <- layer2_coordinate_fixture("NA_fuse")
  expect_error(lt_measurement_eligibility(
    x, matrix("ineligible", 1, 1), matrix("measurement_qc_excluded", 1, 1)
  ), "Non-observed")
})

test_that("FIX001 G non-observed coordinate cannot be explicitly unresolved", {
  x <- layer2_coordinate_fixture("NA_topo")
  expect_error(lt_measurement_eligibility(
    x, matrix("unresolved", 1, 1), matrix("authority_unresolved", 1, 1)
  ), "Non-observed")
})

test_that("FIX001 H observed coordinate cannot be explicitly not_applicable", {
  x <- layer2_coordinate_fixture("observed")
  expect_error(lt_measurement_eligibility(
    x, matrix("not_applicable", 1, 1), matrix("coordinate_not_observed", 1, 1)
  ), "Observed coordinates")
})

test_that("FIX001 I scoped state change preserves the C3 fit mask", {
  x <- small_lt_matrix()
  observed <- x$coord_state == "observed"
  legacy_status <- matrix("unresolved", nrow(x$values), ncol(x$values),
                          dimnames = dimnames(x$values))
  legacy_status[observed & is.finite(x$values) & x$values >= 0] <- "eligible"
  legacy_fit_mask <- observed & legacy_status == "eligible" &
    is.finite(x$values) & x$values >= 0
  domain <- lt_c3_domain(x)
  expect_identical(domain$fit_mask, legacy_fit_mask)
  expect_true(all(domain$eligibility$status[!observed] == "not_applicable"))
})

test_that("FIX001 observed coordinate without admissible payload remains unresolved", {
  values <- matrix(NA_real_, 1, 1, dimnames = list("g1", "b1"))
  state <- matrix("observed", 1, 1, dimnames = dimnames(values))
  reason <- matrix("current_payload_missing", 1, 1, dimnames = dimnames(values))
  x <- lt_matrix(
    values, state, fixture_payload(), value_reason = reason,
    coordinate_provenance = fixture_coordinate_provenance()
  )
  e <- lt_measurement_eligibility(x)
  expect_identical(e$status[[1]], "unresolved")
  expect_identical(e$reason[[1]], "authority_unresolved")
})
