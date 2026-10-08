fixture_payload <- function(origin = "imported", type = "branch_length") {
  list(
    type = type,
    name = if (type == "branch_length") "test_branch_length" else "test_rate",
    units = if (type == "branch_length") "substitutions_per_site" else "unitless",
    scale = "linear",
    origin = origin,
    spec_id = if (origin == "derived") "fixture_rate_spec/1.0.0" else NA_character_
  )
}

fixture_coordinate_provenance <- function() {
  list(
    coordinate_system = "fixture projected-split coordinates",
    reference_tree_sha256 = strrep("a", 64L),
    ordered_branch_ledger_sha256 = strrep("b", 64L),
    ordered_taxon_ledger_sha256 = strrep("c", 64L),
    branch_label_contract = "fixture_B1_B2",
    source_method = "fixture imported coordinate mapping",
    source_software = "fixture-mapper",
    source_version = "1.0.0",
    source_commit = strrep("d", 40L)
  )
}

small_lt_matrix <- function() {
  values <- matrix(
    c(0, NA_real_, 2.5, NA_real_),
    nrow = 2,
    dimnames = list(c("gene_a", "gene_b"), c("B1", "B2"))
  )
  coord_state <- matrix(
    c("observed", "NA_fuse", "observed", "NA_topo"),
    nrow = 2,
    dimnames = dimnames(values)
  )
  lt_matrix(
    values = values,
    coord_state = coord_state,
    payload = fixture_payload(),
    coordinate_provenance = fixture_coordinate_provenance()
  )
}
