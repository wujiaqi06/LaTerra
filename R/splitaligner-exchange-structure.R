# Relational consumer validation is independent of manifest integrity. These
# helpers are private; public import also mandates the versioned upstream validator.
.lt_sx_labels <- function(x, allow_synthetic) {
  l <- x$labels
  .lt_sx_fields(l, c("registry_version", "native_scheme", "native_ledger", "active_scheme", "active_ledger", "migrations"), "labels")
  native <- list(id = "splitalignerr.terminal-first", version = "1")
  legacy <- list(id = "splitaligner.branch-order.legacy", version = "1")
  scheme <- function(s) .lt_sx_need(identical(s, native) || identical(s, legacy), "Unknown label scheme/version; aliases cannot identify scope.", "labels")
  ledger <- function(d, s) {
    .lt_sx_columns(d, c("coordinate_id", "alias", "scheme_id", "scheme_version"), "label ledger")
    .lt_sx_need(identical(names(d), c("coordinate_id", "alias", "scheme_id", "scheme_version")) &&
      identical(d$coordinate_id, x$axes$coordinate_ids) && .lt_sx_keys(d$alias) &&
      all(grepl("^B[1-9][0-9]*$", d$alias)) && all(d$scheme_id == s$id) &&
      all(d$scheme_version == s$version), "Invalid scoped ordered label ledger.", "labels")
  }
  .lt_sx_equal(l$registry_version, x$schema$version, "label registry")
  .lt_sx_equal(l$native_scheme, native, "native label scheme")
  scheme(l$active_scheme); ledger(l$native_ledger, native); ledger(l$active_ledger, l$active_scheme)
  .lt_sx_equal(l$native_ledger$alias, x$coordinates$source_alias, "native source aliases")
  .lt_sx_need(is.list(l$migrations), "Complete ordered migration history is required.", "labels")
  previous <- l$native_ledger; previous_scheme <- native; seen <- character()
  for (m in l$migrations) {
    .lt_sx_fields(m, c("schema_version", "from_scheme", "to_scheme", "reference_tree_sha256",
      "ordered_branch_ledger_sha256", "from_ledger", "to_ledger", "authority", "migration_id"), "migration")
    scheme(m$from_scheme); scheme(m$to_scheme)
    .lt_sx_need(.lt_sx_string(m$migration_id) && !m$migration_id %in% seen &&
      identical(m$schema_version, x$schema$version) && !identical(m$from_scheme, m$to_scheme),
      "Duplicate or invalid migration operation.", "labels")
    .lt_sx_equal(m$from_scheme, previous_scheme, "migration source scheme")
    .lt_sx_equal(m$from_ledger, previous, "migration source ledger")
    .lt_sx_equal(m$reference_tree_sha256, x$identity$reference_tree_sha256, "migration tree identity")
    .lt_sx_equal(m$ordered_branch_ledger_sha256, x$identity$ordered_branch_ledger_sha256, "migration ordered branch identity")
    ledger(m$from_ledger, m$from_scheme); ledger(m$to_ledger, m$to_scheme)
    .lt_sx_need(is.list(m$authority), "Missing migration authority.", "migration")
    if (identical(m$authority$kind, "source_bound_profile")) {
      .lt_sx_fields(m$authority, c("kind", "profile"), "source-bound migration authority")
      p <- m$authority$profile
      .lt_sx_need(inherits(p, "splitaligner_legacy_label_profile") &&
        identical(p$profile_id, "marine302-user-bound-20260908-v1") &&
        identical(p$schema, list(id = "SplitAlignerR.legacy-label-profile",
          version = "0.2.0-development", exchange_version = "0.2.0-development",
          hash_schema = "SAR-C14N-1-SHA256")),
        "Unsupported named legacy profile. Re-export with the admitted source-bound profile; do not hand-match B labels.", "migration")
      # The mandatory upstream public validator reconstructs source bytes and
      # checks the profile, exact tree and from/to ledgers. Do not duplicate it.
    } else {
      .lt_sx_need(identical(m$authority$kind, "synthetic_fixture"),
        "Unknown migration authority; missing source tree/profile cannot be guessed from B labels.", "migration")
      .lt_sx_fields(m$authority, c("kind", "id", "definition"), "synthetic migration authority")
      .lt_sx_need(isTRUE(allow_synthetic) && identical(x$producer$source$provenance_status, "synthetic_fixture") &&
        .lt_sx_string(m$authority$id) && .lt_sx_string(m$authority$definition),
        "Synthetic migration requires explicit opt-in and synthetic producer provenance.", "migration")
    }
    previous <- m$to_ledger; previous_scheme <- m$to_scheme; seen <- c(seen, m$migration_id)
  }
  .lt_sx_equal(l$active_ledger, previous, "active migration ledger")
  .lt_sx_equal(l$active_scheme, previous_scheme, "active migration scheme")
  invisible(x)
}

.lt_sx_coordinates <- function(x) {
  c <- x$coordinates
  .lt_sx_columns(c, c("coordinate_id", "canonical_split_key", "key_schema", "source_alias", "coordinate_type",
    "branch_type", "terminal_taxon", "side_A_taxa", "side_B_taxa", "representation_aliases"), "coordinates")
  .lt_sx_equal(c$coordinate_id, x$axes$coordinate_ids, "ordered coordinate axis")
  .lt_sx_need(.lt_sx_keys(c$canonical_split_key) && .lt_sx_keys(c$source_alias) &&
    all(c$key_schema == x$contract$split_key_schema) && all(c$coordinate_type == "primitive") &&
    all(c$branch_type %in% c("terminal", "internal")), "Invalid scientific coordinate/type/key ledger.", "coordinates")
  for (field in c("side_A_taxa", "side_B_taxa", "representation_aliases")) {
    .lt_sx_need(is.list(c[[field]]) && inherits(c[[field]], "AsIs") &&
      all(vapply(c[[field]], .lt_sx_keys, logical(1))), paste(field, "must retain explicit unique character arrays."), "coordinates")
  }
  .lt_sx_fields(x$tree, c("phylo", "source", "ordered_taxon_ledger"), "tree")
  tree <- x$tree$phylo
  .lt_sx_need(inherits(tree, "phylo") && !inherits(tree, "multiPhylo") &&
    .lt_sx_keys(tree$tip.label) && length(tree$tip.label) >= 3L, "Supply one complete explicit phylo with unique taxa.", "tree")
  .lt_sx_need(!any(grepl(";", tree$tip.label, fixed = TRUE)),
    "Semicolon-bearing taxa cannot be represented losslessly by the current S2 side-string API; no renaming.", "coordinates")
  for (i in seq_len(nrow(c))) {
    a <- c$side_A_taxa[[i]]; b <- c$side_B_taxa[[i]]
    .lt_sx_need(!length(intersect(a, b)) && setequal(c(a, b), tree$tip.label),
      "Split sides must partition the exact tree domain.", "coordinates")
    terminal <- c$branch_type[i] == "terminal"
    expected <- if (terminal) { if (length(a) == 1L) a else b } else NA_character_
    .lt_sx_need(length(expected) == 1L && identical(c$terminal_taxon[i], expected),
      "Terminal identity disagrees with explicit split sides.", "coordinates")
    .lt_sx_need(c$source_alias[i] %in% c$representation_aliases[[i]], "Source alias is absent from representation aliases.", "labels")
  }
  # This is an explicit lossless API conversion of arrays, not split-key parsing.
  data.frame(branch_id = c$coordinate_id, canonical_split_key = c$canonical_split_key,
    side_A_taxa = vapply(c$side_A_taxa, paste, character(1), collapse = ";"),
    side_B_taxa = vapply(c$side_B_taxa, paste, character(1), collapse = ";"),
    branch_type = c$branch_type, terminal_taxon = c$terminal_taxon, stringsAsFactors = FALSE)
}

.lt_sx_preserved <- function(x, d) {
  p <- x$preserved; r <- p$result; genes <- x$axes$gene_ids; aliases <- x$coordinates$source_alias
  .lt_sx_fields(p, c("result", "composite_values", "composite_coordinates", "composite_ledger", "state_ledger", "gene_provenance"), "preserved")
  .lt_sx_need(inherits(r, "splitaligner_result"), "Complete original result is required.")
  required <- c("state_matrix", "numeric_matrix", "state_ledger", "primitive_coordinates", "composite_coordinates",
    "composite_ledger", "gene_provenance", "diagnostics", "species_tree", "conventions", "metadata", "coordinate_table", "input")
  .lt_sx_fields(r, required, "original result")
  for (n in c("composite_coordinates", "composite_ledger", "state_ledger", "gene_provenance"))
    .lt_sx_equal(p[[n]], r[[n]], paste("preserved", n))
  .lt_sx_equal(r$metadata$mode, "fixed", "result mode")
  .lt_sx_equal(r$metadata$core_version, x$producer$core_version, "result core version")
  .lt_sx_equal(r$metadata$state_schema, x$contract$source_state_schema, "result state schema")
  .lt_sx_equal(r$metadata$split_key_schema, x$contract$split_key_schema, "result split schema")
  .lt_sx_need(is.matrix(r$state_matrix) && is.character(r$state_matrix) &&
    identical(dimnames(r$state_matrix), list(genes, aliases)), "Original state matrix has inconsistent ordered axes.", "axis")
  .lt_sx_columns(r$primitive_coordinates, c("coordinate_id", "canonical_split", "branch_type", "side_a_taxa", "side_b_taxa", "primitive_aliases"), "original primitive coordinates")
  fields <- c(coordinate_id = "source_alias", canonical_split = "canonical_split_key", branch_type = "branch_type",
    side_a_taxa = "side_A_taxa", side_b_taxa = "side_B_taxa", primitive_aliases = "representation_aliases")
  for (n in names(fields)) .lt_sx_equal(r$primitive_coordinates[[n]], x$coordinates[[fields[[n]]]], paste("original coordinate", n))
  .lt_sx_columns(r$composite_coordinates, c("coordinate_id", "coordinate_type", "member_count", "member_text", "primitive_members"), "composite coordinates")
  composites <- r$composite_coordinates$coordinate_id
  .lt_sx_need(.lt_sx_keys(composites, TRUE) && !any(composites %in% aliases) &&
    is.matrix(r$numeric_matrix) && is.double(r$numeric_matrix) &&
    identical(dimnames(r$numeric_matrix), list(genes, c(aliases, composites))) &&
    !any(is.nan(r$numeric_matrix) | is.infinite(r$numeric_matrix)), "Invalid primitive/composite numeric axes or payload.", "axis")
  for (i in seq_along(composites)) {
    members <- r$composite_coordinates$primitive_members[[i]]
    .lt_sx_need(.lt_sx_keys(members) && length(members) >= 2L && all(members %in% aliases) &&
      identical(r$composite_coordinates$coordinate_type[i], "composite") &&
      r$composite_coordinates$member_count[i] == length(members), "Invalid composite member definition.")
  }
  v <- r$numeric_matrix[, aliases, drop = FALSE]; dimnames(v) <- dimnames(d$values)
  s <- r$state_matrix; dimnames(s) <- dimnames(d$source_state)
  .lt_sx_equal(d$values, v, "primitive numeric payload")
  .lt_sx_equal(d$source_state, s, "original source states")
  .lt_sx_equal(p$composite_values, r$numeric_matrix[, composites, drop = FALSE], "preserved composite values")
  l <- r$state_ledger
  .lt_sx_columns(l, c("gene_id", "coordinate_id", "projected_split", "projected_side_a_size", "projected_side_b_size",
    "state", "reason_code", "composite_id", "numeric_available", "numeric_value", "numeric_status"), "state ledger")
  .lt_sx_equal(l$gene_id, rep(genes, each = length(aliases)), "state ledger gene-major keys")
  .lt_sx_equal(l$coordinate_id, rep(aliases, times = length(genes)), "state ledger branch keys")
  flatten <- function(v) as.vector(t(v))
  .lt_sx_equal(l$state, flatten(d$source_state), "state ledger states")
  .lt_sx_equal(l$numeric_value, flatten(d$values), "state ledger values")
  .lt_sx_equal(l$numeric_available, is.finite(l$numeric_value), "state ledger numeric availability")
  .lt_sx_equal(l$reason_code, flatten(d$source_reason), "source reason ledger")
  .lt_sx_equal(l$numeric_status, flatten(d$numeric_status), "numeric status ledger")
  fused <- l$state == "NA_fuse"
  .lt_sx_need(all(l$composite_id[fused] %in% composites) && all(is.na(l$composite_id[!fused])), "Fusion provenance is missing or assigned to non-fused cells.")
  cl <- r$composite_ledger
  .lt_sx_columns(cl, c("gene_id", "composite_id", "projected_split", "recovery_status", "numeric_available", "numeric_value", "numeric_status"), "composite ledger")
  keys <- cl[c("gene_id", "composite_id")]
  .lt_sx_need(!anyDuplicated(keys) && all(cl$gene_id %in% genes) && all(cl$composite_id %in% composites), "Invalid composite ledger keys.")
  expected_keys <- unique(data.frame(gene_id = l$gene_id[fused], composite_id = l$composite_id[fused]))
  .lt_sx_need(nrow(keys) == nrow(expected_keys) && all(vapply(seq_len(nrow(keys)), function(i)
    any(expected_keys$gene_id == keys$gene_id[i] & expected_keys$composite_id == keys$composite_id[i]), logical(1))),
    "Composite ledger must exactly cover gene-specific fusion coordinates.")
  .lt_sx_equal(cl$numeric_available, is.finite(cl$numeric_value), "composite availability")
  expected_values <- p$composite_values; expected_values[] <- NA_real_
  if (nrow(cl)) expected_values[cbind(match(cl$gene_id, genes), match(cl$composite_id, composites))] <- cl$numeric_value
  .lt_sx_equal(p$composite_values, expected_values, "composite ledger payload")
  .lt_sx_columns(r$gene_provenance, c("gene_id", "retained_taxon_count", "retained_taxa_key", "retained_taxa"), "gene provenance")
  .lt_sx_equal(r$gene_provenance$gene_id, genes, "gene provenance axis")
  for (i in seq_along(genes)) {
    taxa <- r$gene_provenance$retained_taxa[[i]]
    .lt_sx_need(.lt_sx_keys(taxa) && length(taxa) >= 2L && all(taxa %in% x$tree$phylo$tip.label) &&
      r$gene_provenance$retained_taxon_count[i] == length(taxa), "Invalid retained-taxon provenance.")
  }
  .lt_sx_equal(r$input$gene_ids, genes, "original input gene axis")
  .lt_sx_equal(r$species_tree$coordinates, r$primitive_coordinates, "original species coordinate provenance")
  .lt_sx_need(setequal(r$species_tree$tip_labels, x$tree$phylo$tip.label), "Original species taxon domain differs.")
  .lt_sx_columns(r$coordinate_table, c("coordinate_id", "coordinate_type", "branch_type", "canonical_reference_split",
    "member_count", "member_text", "representation_alias_text", "primitive_members"), "original coordinate table")
  tab <- r$coordinate_table
  .lt_sx_equal(tab$coordinate_id, c(aliases, composites), "coordinate table axis")
  .lt_sx_equal(tab$coordinate_type, c(rep("primitive", length(aliases)), rep("composite", length(composites))), "coordinate table type")
  .lt_sx_equal(tab$branch_type, c(x$coordinates$branch_type, rep("composite", length(composites))), "coordinate table branch type")
  .lt_sx_equal(tab$canonical_reference_split, c(x$coordinates$canonical_split_key, rep(NA_character_, length(composites))), "coordinate table split key")
  .lt_sx_equal(tab$member_count, c(rep(1L, length(aliases)), r$composite_coordinates$member_count), "coordinate table member counts")
  expected_members <- c(lapply(aliases, function(id) id), r$composite_coordinates$primitive_members)
  .lt_sx_equal(unclass(tab$primitive_members), unclass(expected_members), "coordinate table members")
  .lt_sx_columns(r$diagnostics, c("gene_id", "code", "severity", "count", "message"), "diagnostics")
  .lt_sx_need(all(r$diagnostics$gene_id %in% genes) && is.numeric(r$diagnostics$count) &&
    !anyNA(r$diagnostics$count) && all(is.finite(r$diagnostics$count) & r$diagnostics$count >= 0), "Invalid diagnostic provenance.")
  for (field in c("gene_count", "primitive_coordinate_count", "composite_coordinate_count")) {
    expected <- switch(field, gene_count = length(genes), primitive_coordinate_count = length(aliases), composite_coordinate_count = length(composites))
    .lt_sx_equal(r$metadata[[field]], expected, paste("original", field))
  }
  invisible(x)
}

.lt_sx_tree <- function(x, matrix, coordinates) {
  tree <- x$tree$phylo
  .lt_sx_equal(x$tree$ordered_taxon_ledger, data.frame(position = seq_along(tree$tip.label), taxon = tree$tip.label), "ordered taxon ledger")
  source <- x$tree$source
  .lt_sx_need(is.list(source) && .lt_sx_string(source$kind) && .lt_sx_string(source$id), "Explicit tree provenance is required.", "tree")
  if (identical(source$kind, "provided_phylo")) {
    .lt_sx_fields(source, c("kind", "id"), "phylo source")
  } else {
    .lt_sx_need(source$kind %in% c("file", "newick_text"), "Unknown tree-source kind.", "tree")
    .lt_sx_fields(source, c("kind", "id", if (source$kind == "file") "path", "bytes", "sha256"), "tree source")
    if (source$kind == "file") .lt_sx_need(.lt_sx_string(source$path), "Original source path must be recorded, not reopened.", "tree")
    .lt_sx_need(is.raw(source$bytes) && length(source$bytes) > 0L &&
      identical(digest::digest(source$bytes, algo = "sha256", serialize = FALSE), source$sha256), "Tree source bytes/hash disagree.", "tree")
    .lt_sx_need(requireNamespace("ape", quietly = TRUE), "Install ape to validate explicit tree source bytes.", "dependency")
    source_tree <- ape::read.tree(text = rawToChar(source$bytes))
    for (f in c("edge", "tip.label", "Nnode", "edge.length", "node.label", "root.edge"))
      .lt_sx_equal(tree[[f]], source_tree[[f]], paste("tree source geometry", f))
  }
  for (f in c("edge.length", "root.edge")) if (!is.null(tree[[f]]))
    .lt_sx_need(is.numeric(tree[[f]]) && !any(is.nan(tree[[f]]) | is.infinite(tree[[f]])) &&
      length(tree[[f]]) == if (f == "root.edge") 1L else nrow(tree$edge), "Invalid tree geometry lengths.", "tree")
  # Existing LT explicit-taxon-side adapter: no upstream graph/state algorithm.
  join <- .lt_s2_user_coordinates(matrix, tree, coordinates)
  .lt_sx_need(nrow(join) == nrow(tree$edge) && !anyDuplicated(join$tree_edge_index) &&
    setequal(join$tree_edge_index, seq_len(nrow(tree$edge))), "Primitive-to-edge mapping must be exhaustive and bijective.", "tree")
  a <- x$annotation
  .lt_sx_fields(a, c("tree_instance_id", "root_node", "root_has_ordinary_incoming_edge", "root_edge", "edges"), "annotation")
  instance <- paste0("tree-instance2:", x$identity$reference_tree_sha256)
  .lt_sx_equal(a$tree_instance_id, instance, "annotation tree instance")
  .lt_sx_equal(a$root_node, as.integer(setdiff(unique(tree$edge[, 1]), tree$edge[, 2])), "annotation root")
  .lt_sx_equal(a$root_has_ordinary_incoming_edge, FALSE, "root incoming-edge policy")
  .lt_sx_equal(a$root_edge, tree$root.edge, "separate root geometry")
  j <- match(seq_len(nrow(tree$edge)), join$tree_edge_index)
  c <- x$coordinates; ids <- join$branch_id[j]; ix <- match(ids, c$coordinate_id)
  expected <- data.frame(tree_instance_id = instance, coordinate_id = ids,
    canonical_split_key = c$canonical_split_key[ix], source_alias = c$source_alias[ix],
    display_alias = x$labels$active_ledger$alias[match(ids, x$labels$active_ledger$coordinate_id)],
    edge_row = seq_len(nrow(tree$edge)), parent_node = as.integer(tree$edge[, 1]),
    child_node = as.integer(tree$edge[, 2]), node = as.integer(tree$edge[, 2]),
    branch_type = c$branch_type[ix], terminal_taxon = c$terminal_taxon[ix],
    tree_branch_length = if (is.null(tree$edge.length)) rep(NA_real_, nrow(tree$edge)) else tree$edge.length)
  .lt_sx_equal(a$edges, expected, "tree-instance annotation edge table")
  invisible(join)
}

.lt_sx_adapt <- function(x, data, allow_synthetic_migration = FALSE) {
  coordinates <- .lt_sx_coordinates(x)
  .lt_sx_labels(x, allow_synthetic_migration)
  .lt_sx_preserved(x, data)
  payload <- x$contract$payload; payload$origin <- "imported"
  provenance <- list(schema = x$schema, producer = x$producer, identity = x$identity,
    upstream_payload = x$contract$payload, manifest = x$manifest,
    original_orientation = x$axes$orientation, normalization = if (x$axes$branch_margin == 1L) "transpose_all_six_arrays" else "none",
    original_evidence_location = "returned_bundle$exchange (save the complete returned bundle)",
    acceptance = "VALIDATED_TYPED_IMPORT_NOT_RELEASE_CERTIFICATION", scientific_analysis_run = FALSE)
  cp <- list(coordinate_system = x$identity$coordinate_system_id,
    reference_tree_sha256 = x$identity$reference_tree_sha256,
    ordered_branch_ledger_sha256 = x$identity$ordered_branch_ledger_sha256,
    ordered_taxon_ledger_sha256 = x$identity$ordered_taxon_ledger_sha256,
    branch_label_contract = paste(x$labels$active_scheme$id, x$labels$active_scheme$version, sep = "/"),
    source_method = "SplitAlignerR.exchange/import_only", source_software = x$producer$package,
    source_version = x$producer$version, source_commit = x$producer$source$source_commit)
  matrix <- lt_matrix(data$values, data$coord_state, payload, value_reason = data$value_reason,
    gene_ids = x$axes$gene_ids, branch_ids = x$axes$coordinate_ids, coordinate_provenance = cp,
    branch_display_labels = x$labels$active_ledger$alias, metadata = list(splitaligner_import = provenance))
  .lt_sx_tree(x, matrix, coordinates)
  list(matrix = matrix, tree = x$tree$phylo, coordinates = coordinates, exchange = x, provenance = provenance)
}
