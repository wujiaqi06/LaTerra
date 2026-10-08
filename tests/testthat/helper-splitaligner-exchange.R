# The stored small fixture is genuine producer output with its original file SHA.
# Mutation helpers below are TEST-ONLY derived envelopes, never claimed as new
# exporter output, portable identities, or historical migration authority.
sx_fixture_original <- function() readRDS(test_path("fixtures", "splitaligner_exchange_v020.rds"))
# Use the pinned producer's encoder only in mutation fixtures, not a second
# consumer encoder or a public validation bypass.
sx_h <- function(x) getFromNamespace(".sx_hash", "SplitAlignerR")(x)
sx_seal_test <- function(x) {
  components <- c("schema", "producer", "contract", "axes", "identity", "tree", "coordinates", "labels", "data", "preserved", "annotation")
  h <- vapply(x[components], sx_h, character(1))
  x$manifest <- list(hash_schema = "SAR-C14N-1-SHA256", component_sha256 = h, content_sha256 = sx_h(h))
  x
}
sx_local_test <- function(x = sx_fixture_original()) {
  testthat::skip_if_not_installed("SplitAlignerR")
  stopifnot(!length(x$labels$migrations))
  domain <- sx_h(sort(enc2utf8(x$tree$phylo$tip.label), method = "radix"))
  old <- x$coordinates$coordinate_id
  ids <- unname(vapply(x$coordinates$canonical_split_key, function(key) paste0("SAR2:", sx_h(list(
    schema = x$contract$split_key_schema, taxon_domain_sha256 = domain, canonical_key = enc2utf8(key)))), character(1)))
  x$coordinates$coordinate_id <- ids; x$axes$coordinate_ids <- ids
  for (n in names(x$data)) dimnames(x$data[[n]])[[x$axes$branch_margin]] <- ids
  for (n in c("native_ledger", "active_ledger")) x$labels[[n]]$coordinate_id <- ids
  x$annotation$edges$coordinate_id <- ids[match(x$annotation$edges$coordinate_id, old)]
  x$identity <- list(coordinate_system_id = paste0("SAR-CS2:", sx_h(list(taxon_domain_sha256 = domain,
    coordinate_ids = sort(ids, method = "radix")))), taxon_domain_sha256 = domain,
    reference_tree_sha256 = sx_h(x$tree$phylo), ordered_branch_ledger_sha256 = sx_h(x$coordinates),
    ordered_taxon_ledger_sha256 = sx_h(x$tree$ordered_taxon_ledger), ordered_label_ledger_sha256 = sx_h(x$labels$active_ledger))
  instance <- paste0("tree-instance2:", x$identity$reference_tree_sha256)
  x$annotation$tree_instance_id <- instance; x$annotation$edges$tree_instance_id <- rep(instance, nrow(x$annotation$edges))
  sx_seal_test(x)
}
sx_import_test <- function(x, ...) {
  f <- tempfile(fileext = ".rds"); on.exit(unlink(f)); saveRDS(x, f, version = 3L)
  lt_read_splitaligner_exchange(f, ...)
}
sx_cell_test <- function(x, gene, branch, state, value, status, reason) {
  stopifnot(x$axes$branch_margin == 2L)
  x$data$values[gene, branch] <- value
  x$data$source_state[gene, branch] <- state
  x$data$coord_state[gene, branch] <- if (state == "mapped") "observed" else state
  x$data$source_reason[gene, branch] <- reason
  x$data$numeric_status[gene, branch] <- status
  x$data$value_reason[gene, branch] <- if (state == "mapped" && is.na(value)) status else NA_character_
  r <- x$preserved$result
  r$numeric_matrix[gene, branch] <- value; r$state_matrix[gene, branch] <- state
  i <- (gene - 1L) * ncol(r$state_matrix) + branch
  r$state_ledger$state[i] <- state; r$state_ledger$numeric_value[i] <- value
  r$state_ledger$numeric_available[i] <- is.finite(value)
  r$state_ledger$reason_code[i] <- reason; r$state_ledger$numeric_status[i] <- status
  x$preserved$result <- r; x$preserved$state_ledger <- r$state_ledger
  sx_seal_test(x)
}
sx_migration_test <- function(x = sx_local_test()) {
  scope <- list(id = "splitaligner.branch-order.legacy", version = "1")
  target <- x$labels$active_ledger; target$alias <- rev(target$alias); target$scheme_id <- scope$id
  m <- list(schema_version = x$schema$version, from_scheme = x$labels$active_scheme, to_scheme = scope,
    reference_tree_sha256 = x$identity$reference_tree_sha256,
    ordered_branch_ledger_sha256 = x$identity$ordered_branch_ledger_sha256,
    from_ledger = x$labels$active_ledger, to_ledger = target,
    authority = list(kind = "synthetic_fixture", id = "test-only-alias-migration", definition = "Synthetic reversal; not a Marine migration."))
  m$migration_id <- paste0("label-migration2:", sx_h(m))
  x$labels$active_scheme <- scope; x$labels$active_ledger <- target; x$labels$migrations <- list(m)
  x$identity$ordered_label_ledger_sha256 <- sx_h(target)
  x$annotation$edges$display_alias <- target$alias[match(x$annotation$edges$coordinate_id, target$coordinate_id)]
  sx_seal_test(x)
}
sx_gene_order_test <- function(x, index, gene_ids = x$axes$gene_ids[index]) {
  # Complete explicitly keyed derivation, including original result provenance.
  r <- x$preserved$result; old <- x$axes$gene_ids
  for (n in names(x$data)) {
    x$data[[n]] <- x$data[[n]][index, , drop = FALSE]; rownames(x$data[[n]]) <- gene_ids
  }
  for (n in c("state_matrix", "numeric_matrix")) {
    r[[n]] <- r[[n]][index, , drop = FALSE]; rownames(r[[n]]) <- gene_ids
  }
  for (n in c("state_ledger", "composite_ledger", "gene_provenance", "diagnostics")) {
    r[[n]] <- do.call(rbind, lapply(seq_along(index), function(i) {
      a <- r[[n]][r[[n]]$gene_id == old[index[i]], , drop = FALSE]
      a$gene_id <- rep(gene_ids[i], nrow(a)); a
    }))
    rownames(r[[n]]) <- NULL
  }
  r$input$gene_ids <- gene_ids; r$metadata$gene_count <- length(gene_ids)
  x$axes$gene_ids <- gene_ids; x$preserved$result <- r
  for (n in c("state_ledger", "composite_ledger", "gene_provenance")) x$preserved[[n]] <- r[[n]]
  x$preserved$composite_values <- r$numeric_matrix[, r$composite_coordinates$coordinate_id, drop = FALSE]
  sx_seal_test(x)
}
sx_composite_test <- function() {
  x <- sx_local_test()
  for (j in 6:7) x <- sx_cell_test(x, 1, j, "NA_fuse", NA_real_, "composite_coordinate_only", "TEST_ONLY_FUSION")
  r <- x$preserved$result; cid <- "F[B6|B7]"; members <- c("B6", "B7")
  r$state_ledger$composite_id[6:7] <- cid
  r$composite_coordinates <- data.frame(coordinate_id = cid, coordinate_type = "composite",
    member_count = 2L, member_text = "B6|B7", primitive_members = I(list(members)))
  r$composite_ledger <- data.frame(gene_id = x$axes$gene_ids[1], composite_id = cid,
    projected_split = "opaque_test_composite_projection", recovery_status = "mapped",
    numeric_available = TRUE, numeric_value = 9, numeric_status = "finite_numeric")
  r$numeric_matrix <- cbind(r$numeric_matrix, c(9, NA_real_)); colnames(r$numeric_matrix)[8] <- cid
  r$coordinate_table <- rbind(r$coordinate_table, data.frame(coordinate_id = cid, coordinate_type = "composite",
    branch_type = "composite", canonical_reference_split = NA_character_, member_count = 2L,
    member_text = "B6|B7", representation_alias_text = NA_character_, primitive_members = I(list(members))))
  r$metadata$composite_coordinate_count <- 1L
  x$preserved$result <- r
  for (n in c("state_ledger", "composite_coordinates", "composite_ledger")) x$preserved[[n]] <- r[[n]]
  x$preserved$composite_values <- r$numeric_matrix[, cid, drop = FALSE]
  sx_seal_test(x)
}
