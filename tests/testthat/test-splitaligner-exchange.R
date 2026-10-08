test_that("old producer file is immutable and always explicitly rejected", {
  f <- test_path("fixtures", "splitaligner_exchange_v010.rds")
  expect_identical(digest::digest(file = f, algo = "sha256", serialize = FALSE),
    "1855466b4add699d7ae0a56a30acfa7921372d6dce9a589cf51f93ff4180de5e")
  before <- Sys.getlocale()
  expect_error(lt_read_splitaligner_exchange(f), "H01", class = "lt_error_exchange_version")
  expect_identical(Sys.getlocale(), before)
  expect_false("verify" %in% names(formals(lt_read_splitaligner_exchange)))
})

test_that("real successor passes mandatory upstream validation without locale mutation", {
  skip_if_not_installed("SplitAlignerR")
  f <- test_path("fixtures", "splitaligner_exchange_v020.rds")
  expect_identical(digest::digest(file = f, algo = "sha256", serialize = FALSE),
    "b557d98e93d1388a87587b22b9e11a19df9b2ad6ded3927453e2eea712a1b7d5")
  before <- Sys.getlocale(); out <- lt_read_splitaligner_exchange(f)
  expect_identical(out$exchange, sx_fixture_original())
  expect_identical(out$provenance$validation_backend$public_validator, "SplitAlignerR::validate_splitaligner_exchange")
  expect_identical(out$provenance$validation_backend$version, "0.1.0.9002")
  expect_identical(Sys.getlocale(), before)
})

test_that("local test envelope returns exact usable inputs and preserves evidence", {
  x <- sx_local_test(); out <- sx_import_test(x)
  expect_s3_class(out$matrix, "lt_matrix")
  expect_identical(out$matrix$values, x$data$values)
  expect_identical(out$matrix$values[1, 1], unname(x$data$values[1, 1]))
  expect_identical(out$matrix$coord_state, x$data$coord_state)
  expect_identical(out$matrix$value_reason, x$data$value_reason)
  expect_identical(out$matrix$branch_ids, x$axes$coordinate_ids)
  expect_identical(out$matrix$branch_display_labels, x$labels$active_ledger$alias)
  expect_identical(out$matrix$payload$origin, "imported")
  expect_identical(out$provenance$upstream_payload, x$contract$payload)
  expect_identical(out$exchange, x)
  expect_identical(out$tree, x$tree$phylo)
  expect_identical(out$coordinates$branch_id, x$coordinates$coordinate_id)
  expect_identical(LaTerra:::.lt_s2_user_coordinates(out$matrix, out$tree, out$coordinates)$branch_id, x$axes$coordinate_ids)
  expect_false(out$provenance$scientific_analysis_run)
  out$matrix$branch_display_labels <- paste0("editable ", 1:7)
  out$matrix$metadata$user_note <- "Ordinary labels remain editable, not a frozen hash validity barrier."
  expect_silent(validate_lt_matrix(out$matrix))
  f <- tempfile(); on.exit(unlink(f)); saveRDS(out, f, version = 3L)
  expect_identical(readRDS(f), out)
})

test_that("every original array rejects independent axis edits before construction", {
  x <- sx_local_test()
  for (n in names(x$data)) {
    y <- x; rownames(y$data[[n]]) <- rev(rownames(y$data[[n]]))
    expect_error(sx_import_test(sx_seal_test(y)), "original ordered axes", class = "lt_error_exchange_axis")
    y <- x; colnames(y$data[[n]]) <- rev(colnames(y$data[[n]]))
    expect_error(sx_import_test(sx_seal_test(y)), "original ordered axes", class = "lt_error_exchange_axis")
    y <- x; dimnames(y$data[[n]]) <- NULL
    expect_error(sx_import_test(sx_seal_test(y)), class = "lt_error_exchange_axis")
    y <- x; colnames(y$data[[n]])[2] <- colnames(y$data[[n]])[1]
    expect_error(sx_import_test(sx_seal_test(y)), class = "lt_error_exchange_axis")
    y <- x; y$data[[n]] <- y$data[[n]][, -1, drop = FALSE]
    expect_error(sx_import_test(sx_seal_test(y)), class = "lt_error_exchange_axis")
  }
  y <- x; y$axes$gene_ids[2] <- y$axes$gene_ids[1]
  expect_error(sx_import_test(sx_seal_test(y)), class = "lt_error_exchange_axis")
  y <- x; y$axes$coordinate_ids[2] <- y$axes$coordinate_ids[1]
  expect_error(sx_import_test(sx_seal_test(y)), class = "lt_error_exchange_axis")
  y <- x; y$data$values <- matrix(as.integer(y$data$values), 2, 7, dimnames = dimnames(y$data$values))
  expect_error(sx_import_test(sx_seal_test(y)), "storage type")
})

test_that("only explicit common transpose is normalized, preserving original envelope", {
  x <- sx_local_test(); y <- x
  y$data <- lapply(y$data, t); y$axes$branch_margin <- 1L; y$axes$orientation <- "branch_by_gene"
  y <- sx_seal_test(y); a <- sx_import_test(x); b <- sx_import_test(y)
  for (n in c("values", "coord_state", "value_reason")) expect_identical(a$matrix[[n]], b$matrix[[n]])
  expect_identical(b$exchange, y)
  expect_identical(b$provenance$normalization, "transpose_all_six_arrays")
  y$data$value_reason <- t(y$data$value_reason)
  expect_error(sx_import_test(sx_seal_test(y)), class = "lt_error_exchange_axis")
  y <- x; y$axes$branch_margin <- NULL
  expect_error(sx_import_test(sx_seal_test(y)), "axes")
})

test_that("unavailable payloads preserve state reasons rather than fabricate absence", {
  x <- sx_local_test()
  for (status in c("missing_branch_length", "software_failure_marker")) {
    y <- sx_cell_test(x, 1, 2, "mapped", NA_real_, status, "PROJECTED_SPLIT_RECOVERED")
    a <- sx_import_test(y)
    expect_identical(a$matrix$coord_state[1, 2], "observed")
    expect_true(is.na(a$matrix$values[1, 2]))
    expect_identical(a$matrix$value_reason[1, 2], status)
    expect_identical(a$exchange$preserved, y$preserved)
  }
  for (state in c("NA_struct", "NA_topo")) {
    y <- sx_cell_test(x, 1, 2, state, NA_real_, "not_applicable", "EXPLICIT_TEST_SOURCE_REASON")
    a <- sx_import_test(y)
    expect_identical(a$matrix$coord_state[1, 2], state)
    expect_true(is.na(a$matrix$value_reason[1, 2]))
    expect_identical(a$exchange$data$source_reason[1, 2], "EXPLICIT_TEST_SOURCE_REASON")
  }
  y <- sx_cell_test(x, 1, 2, "mapped", NA_real_, "not_applicable", "EXPLICIT_TEST_SOURCE_REASON")
  expect_error(sx_import_test(y), "Mapped-unavailable")
  y <- x; y$data$coord_state[1, 2] <- "unknown"
  expect_error(sx_import_test(sx_seal_test(y)), "state mapping")
  y <- x; y$data$source_state[1, 2] <- NA_character_
  expect_error(sx_import_test(sx_seal_test(y)), "source state")
  y <- x; y$data$value_reason[1, 2] <- "invented"
  expect_error(sx_import_test(sx_seal_test(y)), "value reason")
  for (v in c(-1, NaN, Inf, -Inf)) {
    y <- x; y$data$values[1, 2] <- v
    expect_error(sx_import_test(sx_seal_test(y)), class = "lt_error_exchange_payload")
  }
})

test_that("resealing does not authorize mismatched preserved values and ledgers", {
  x <- sx_local_test()
  y <- x; y$data$values[1, 2] <- 999
  expect_error(sx_import_test(sx_seal_test(y)), "primitive numeric payload")
  y <- x; y$data$source_reason[1, 2] <- "altered"
  expect_error(sx_import_test(sx_seal_test(y)), "source reason ledger")
  y <- x; y$preserved$gene_provenance <- NULL
  expect_error(sx_import_test(sx_seal_test(y)), "preserved")
  y <- x; y$preserved$result$state_ledger$gene_id[1] <- "wrong"
  y$preserved$state_ledger <- y$preserved$result$state_ledger
  expect_error(sx_import_test(sx_seal_test(y)), "state ledger gene-major keys")
  y <- x; y$preserved$result$numeric_matrix[1, 1] <- 100
  expect_error(sx_import_test(sx_seal_test(y)), "primitive numeric payload")
  y <- x; y$preserved$result$gene_provenance$retained_taxon_count[1] <- 4L
  y$preserved$gene_provenance <- y$preserved$result$gene_provenance
  expect_error(sx_import_test(sx_seal_test(y)), "retained-taxon")
})

test_that("scoped B aliases cannot substitute for coordinates or migration history", {
  x <- sx_local_test()
  y <- x; y$labels$active_ledger$alias <- rev(y$labels$active_ledger$alias)
  y <- sx_local_test(y)
  expect_error(sx_import_test(y), "active migration ledger")
  y <- x; y$labels$active_scheme <- list(id = "splitaligner.branch-order.legacy", version = "1")
  y$labels$active_ledger$scheme_id <- y$labels$active_scheme$id; y <- sx_local_test(y)
  expect_error(sx_import_test(y), "active migration ledger")
  y <- x; y$labels$native_scheme$version <- "unknown"
  expect_error(sx_import_test(sx_seal_test(y)), "native label scheme")
  y <- x; y$labels$active_scheme$version <- "unknown"
  expect_error(sx_import_test(sx_seal_test(y)), "Unknown label scheme")
  y <- x; y$coordinates$coordinate_id[1] <- "B1"
  expect_error(sx_import_test(sx_seal_test(y)), "coordinate")
})

test_that("tree instance, split-side and root boundaries cannot be resealed away", {
  x <- sx_local_test()
  y <- x; y$annotation$edges$node[1] <- y$annotation$root_node
  expect_error(sx_import_test(sx_seal_test(y)), "annotation edge table")
  y <- x; y$annotation$root_has_ordinary_incoming_edge <- TRUE
  expect_error(sx_import_test(sx_seal_test(y)), "root incoming-edge")
  y <- x; y$annotation$edges$tree_branch_length[1] <- 123
  expect_error(sx_import_test(sx_seal_test(y)), "annotation edge table")
  y <- x; y$tree$source$bytes[1] <- as.raw(65)
  expect_error(sx_import_test(sx_seal_test(y)), "source bytes/hash")
  y <- x; y$coordinates$side_A_taxa[[1]] <- "NOT_A_TAXON"; y <- sx_local_test(y)
  expect_error(sx_import_test(y), "partition the exact tree domain")
  y <- x; y$coordinates$branch_type[1] <- "internal"; y <- sx_local_test(y)
  expect_error(sx_import_test(y), "Terminal identity")
  y <- x; y$tree$phylo$tip.label[2] <- y$tree$phylo$tip.label[1]; y <- sx_local_test(y)
  expect_error(sx_import_test(y), "unique taxa")
  y <- x; y$tree$phylo$tip.label[1] <- "A;ambiguous"; y <- sx_local_test(y)
  expect_error(sx_import_test(y), "Semicolon-bearing")
  y <- x; y$tree$source <- list(kind = "file", id = "relocated/source", path = "/not/reopened/species.nwk",
    bytes = x$tree$source$bytes, sha256 = x$tree$source$sha256)
  y <- sx_seal_test(y)
  expect_identical(sx_import_test(y)$tree, x$tree$phylo)
})

test_that("new schema, mode, payload and undeclared provenance fail closed", {
  x <- sx_local_test()
  y <- x; y$schema$version <- "0.1.1-development"
  expect_error(sx_import_test(sx_seal_test(y)), class = "lt_error_exchange_version")
  for (mode in c("free", "paired")) {
    y <- x; y$contract$mode <- mode
    expect_error(sx_import_test(sx_seal_test(y)), class = "lt_error_exchange_mode")
  }
  for (field in c("type", "scale")) {
    y <- x; y$contract$payload[[field]] <- "unsupported"
    expect_error(sx_import_test(sx_seal_test(y)), class = "lt_error_exchange_payload")
  }
  y <- x; y$contract$state_schema <- "unknown"
  expect_error(sx_import_test(sx_seal_test(y)), "state_schema")
  y <- x; y$producer$source$provenance_status <- "historical_certified"
  expect_error(sx_import_test(sx_seal_test(y)), "provenance")
  expect_error(lt_read_splitaligner_exchange(tempfile()), class = "lt_error_exchange_file")
})

test_that("synthetic migrations require complete authority chain and explicit opt-in", {
  x <- sx_migration_test()
  expect_error(sx_import_test(x), "opt-in", class = "lt_error_exchange_migration")
  a <- sx_import_test(x, allow_synthetic_migration = TRUE)
  expect_identical(a$matrix$values, x$data$values)
  expect_identical(a$matrix$branch_ids, x$axes$coordinate_ids)
  expect_identical(a$matrix$branch_display_labels, rev(x$labels$native_ledger$alias))
  expect_identical(a$exchange$labels$migrations, x$labels$migrations)
  y <- x; y$labels$migrations <- c(y$labels$migrations, y$labels$migrations)
  expect_error(sx_import_test(sx_seal_test(y), allow_synthetic_migration = TRUE), "Duplicate")
  y <- x; m <- y$labels$migrations[[1]]; m$authority <- list(kind = "source_bound_profile", id = "NOT_ADMITTED_MARINE")
  m$migration_id <- paste0("label-migration2:", sx_h(m[setdiff(names(m), "migration_id")]))
  y$labels$migrations[[1]] <- m
  expect_error(sx_import_test(sx_seal_test(y), allow_synthetic_migration = TRUE), "source-bound migration authority")
  y <- x; m <- y$labels$migrations[[1]]; m$reference_tree_sha256 <- paste(rep("0", 64), collapse = "")
  m$migration_id <- paste0("label-migration2:", sx_h(m[setdiff(names(m), "migration_id")]))
  y$labels$migrations[[1]] <- m
  expect_error(sx_import_test(sx_seal_test(y), allow_synthetic_migration = TRUE), "migration tree identity")
  y <- x; y$producer$source$provenance_status <- "development_unreleased"
  expect_error(sx_import_test(sx_seal_test(y), allow_synthetic_migration = TRUE), "synthetic producer")
})

test_that("composite evidence stays complete without becoming an ordinary edge", {
  x <- sx_composite_test(); a <- sx_import_test(x)
  expect_identical(dim(a$matrix$values), c(2L, 7L))
  expect_identical(dim(a$exchange$preserved$composite_values), c(2L, 1L))
  expect_identical(a$matrix$coord_state[1, 6:7], x$data$coord_state[1, 6:7])
  expect_identical(a$exchange$preserved, x$preserved)
  expect_false(any(grepl("^F\\[", a$matrix$branch_ids)))
  expect_identical(nrow(a$exchange$annotation$edges), 7L)
  y <- x; y$preserved$composite_values[1, 1] <- 10
  expect_error(sx_import_test(sx_seal_test(y)), "composite values")
  y <- x; y$preserved$result$state_ledger$composite_id[6] <- "unknown"
  y$preserved$state_ledger <- y$preserved$result$state_ledger
  expect_error(sx_import_test(sx_seal_test(y)), "Fusion provenance")
  y <- x; y$preserved$result$coordinate_table$coordinate_type[8] <- "primitive"
  expect_error(sx_import_test(sx_seal_test(y)), "coordinate table type")
})

test_that("gene reordering is source-keyed and square transpose is never inferred", {
  x <- sx_local_test(); y <- sx_gene_order_test(x, 2:1)
  expect_identical(sx_import_test(y)$matrix$values, x$data$values[2:1, , drop = FALSE])
  y <- sx_gene_order_test(x, rep(1:2, length.out = 7), paste0("test_gene", 1:7))
  a <- sx_import_test(y)
  z <- y; z$data <- lapply(z$data, t); z$axes$branch_margin <- 1L; z$axes$orientation <- "branch_by_gene"
  b <- sx_import_test(sx_seal_test(z))
  expect_identical(dim(a$matrix$values), c(7L, 7L))
  expect_identical(a$matrix$values, b$matrix$values)
  expect_identical(a$matrix$coord_state, b$matrix$coord_state)
  z$axes <- y$axes
  expect_error(sx_import_test(sx_seal_test(z)), class = "lt_error_exchange_axis")
  z <- y; z$data$values <- t(z$data$values)
  expect_error(sx_import_test(sx_seal_test(z)), class = "lt_error_exchange_axis")
})

test_that("public annotation handles reordered display trees without replacing source reference", {
  x <- sx_local_test(); order <- 7:1; y <- x
  y$tree$phylo$edge <- y$tree$phylo$edge[order, , drop = FALSE]
  y$tree$phylo$edge.length <- y$tree$phylo$edge.length[order]
  attr(y$tree$phylo, "user_extra") <- list(note = "must survive import", flag = TRUE)
  y$tree$source <- list(kind = "provided_phylo", id = "explicit_test_reordered_tree")
  y <- sx_local_test(y)
  expect_error(sx_import_test(y), "annotation edge table")
  y$annotation$edges <- y$annotation$edges[order, , drop = FALSE]
  rownames(y$annotation$edges) <- NULL
  y$annotation$edges$edge_row <- seq_len(7L)
  y <- sx_seal_test(y)
  # Replacing the authoritative reference would also need correct native
  # source/result provenance. The producer owns that check, not LT B numbering.
  expect_error(sx_import_test(y), "native ordered label ledger", class = "lt_error_exchange_upstream")
  a <- sx_import_test(x)
  annotation <- SplitAlignerR::annotate_splitaligner_exchange(a$exchange, "gene01", tree = y$tree$phylo)
  expect_identical(annotation$gene_branch_value, as.numeric(a$matrix$values["gene01", annotation$coordinate_id]))
  expect_identical(annotation$child_node, as.integer(y$tree$phylo$edge[, 2]))
  expect_identical(a$tree, x$tree$phylo)
  y <- x; y$tree$source <- list(kind = "provided_phylo", id = "invalid_node_graph")
  y$tree$phylo$edge[1, 2] <- y$tree$phylo$edge[2, 2]; y <- sx_local_test(y)
  expect_error(sx_import_test(y), "node/edge structure")
})
