#' La Terra manifest schema
#'
#' Returns the required future audit fields without inspecting files or
#' computing hashes.
#'
#' @return A named character vector mapping manifest fields to container types.
#' @export
lt_manifest_schema <- function() {
  c(
    input_sha256 = "data.frame(role, path, sha256)",
    configuration_sha256 = "character(1)",
    software_version = "named list",
    environment = "named list",
    ordered_gene_ledger = "data.frame(order, gene_id)",
    ordered_branch_ledger = "data.frame(order, branch_id)",
    taxon_domain_ledger = "data.frame(order, taxon_id, eligible, exclusion_reason)",
    eligibility_masks = "named list",
    specification_ledger = "data.frame(specification_id, kind, version, definition, sha256)",
    filter_ledger = "data.frame(order, stage, decision_id, target_type, target_id, decision, reason)",
    seed_ledger = "data.frame(scope, seed, rng_kind, normal_kind, sample_kind, r_version)",
    fold_ledger = "data.frame(fold_id, taxon_id, role, group)",
    commands = "character vector",
    output_hashes = "data.frame(role, path, sha256)"
  )
}
