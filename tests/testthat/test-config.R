test_that("reference recipes parse and remain non-executing", {
  recipes <- c("marine_2026_legacy.yaml", "aquatic_2026_legacy.yaml")

  for (recipe in recipes) {
    path <- system.file("recipes", recipe, package = "LaTerra")
    config <- lt_read_config(path)
    expect_false(config$execution$enabled)
    expect_false(config$execution$scientific_algorithms_implemented)
    expect_setequal(config$matrix$coord_state$allowed_states, lt_allowed_states())
    expect_identical(config$matrix$payload$type, "branch_length")
    expect_identical(config$matrix$payload$origin, "imported")
    expect_identical(config$matrix$coord_state$mode, "separate_input")
    expect_identical(config$matrix$value_reason$mode, "none")
  }
})

test_that("legacy recipes retain method-specific scientific choices", {
  marine <- lt_read_config(system.file(
    "recipes", "marine_2026_legacy.yaml", package = "LaTerra"
  ))
  aquatic <- lt_read_config(system.file(
    "recipes", "aquatic_2026_legacy.yaml", package = "LaTerra"
  ))

  expect_identical(marine$rate$method, "legacy_gbi_ratio")
  expect_equal(marine$rate$parameters$within_gene_upper_tail_quantile, 0.975)
  expect_identical(marine$model$method, "lasso_binomial")
  expect_identical(marine$model$parameters$engine, "glmnet")
  expect_identical(
    vapply(marine$null_model, `[[`, character(1), "method"),
    c("positive_count_permutation", "predictor_turnover")
  )
  expect_equal(aquatic$trait$coding$excluded_intermediate, 0.5)
  expect_false(aquatic$null_model[[1]]$enabled)
})

test_that("marine authority roles remain distinct", {
  marine <- lt_read_config(system.file(
    "recipes", "marine_2026_legacy.yaml", package = "LaTerra"
  ))
  roles <- vapply(marine$authority$records, `[[`, character(1), "role")
  expect_identical(
    roles,
    c(
      "scientific_specification", "executable_reference", "data_result",
      "source_code_identity"
    )
  )
  expect_identical(
    marine$authority$records[[4]]$peeled_commit,
    "d93e70ee9fdf3ebaaa3ac4ac787856549bfa2a24"
  )
})

test_that("generic plain-matrix configuration is accepted", {
  config <- lt_read_config(system.file(
    "recipes", "marine_2026_legacy.yaml", package = "LaTerra"
  ))
  config$recipe_id <- "generic_plain_tsv"
  config$authority <- NULL
  config$matrix$payload <- list(
    type = "rate_representation", name = "external_rate", units = "unitless",
    scale = "linear", origin = "imported", spec_id = "external/spec/1.0.0"
  )
  config$matrix$coord_state <- list(
    mode = "all_observed",
    allowed_states = lt_allowed_states(),
    numeric_zero_is_observed = TRUE
  )
  config$matrix$value_reason <- list(mode = "none")
  config$model <- list(
    enabled = FALSE,
    method = "not_requested",
    contract_version = "generic/1.0.0",
    parameters = list()
  )
  config$null_model <- NULL

  for (format in c("tsv", "csv", "rds")) {
    config$matrix$input <- list(
      path = paste0("my_matrix.", format), sha256 = NULL, format = format,
      gene_axis = "rows", branch_axis = "columns", gene_id_column = "gene"
    )
    expect_invisible(validate_lt_config(config))
  }
})

test_that("missing coordinates require an explicit coordinate-state source", {
  config <- lt_read_config(system.file(
    "recipes", "marine_2026_legacy.yaml", package = "LaTerra"
  ))
  config$matrix$coord_state <- list(
    mode = "separate_input",
    allowed_states = lt_allowed_states(),
    numeric_zero_is_observed = TRUE
  )
  expect_error(
    validate_lt_config(config),
    "requires an explicit coordinate-state source"
  )
})

test_that("generic modules do not impose a marine model architecture", {
  config <- lt_read_config(system.file(
    "recipes", "marine_2026_legacy.yaml", package = "LaTerra"
  ))
  config$model <- list(
    enabled = TRUE,
    method = "elastic_net",
    contract_version = "example_elastic_net/1.0.0",
    parameters = list(mixing = 0.5)
  )
  config$null_model <- list()
  expect_invisible(validate_lt_config(config))

  schema_path <- system.file("schema", "lt_recipe.schema.yaml", package = "LaTerra")
  schema_text <- paste(readLines(schema_path, warn = FALSE), collapse = "\n")
  expect_false(grepl("glmnet|LASSO|positive_count_permutation|predictor_turnover", schema_text))
})

test_that("YAML configuration round-trips", {
  source <- system.file("recipes", "marine_2026_legacy.yaml", package = "LaTerra")
  config <- lt_read_config(source)
  destination <- tempfile(fileext = ".yaml")

  lt_write_config(config, destination)
  restored <- lt_read_config(destination)
  expect_equal(restored, config)
})

test_that("machine-readable schemas are valid YAML containers", {
  recipe_schema <- yaml::read_yaml(system.file(
    "schema", "lt_recipe.schema.yaml", package = "LaTerra"
  ))
  object_schema <- yaml::read_yaml(system.file(
    "schema", "lt_objects.schema.yaml", package = "LaTerra"
  ))

  expect_identical(recipe_schema$type, "object")
  expect_identical(recipe_schema$properties$schema_version$const, "1.1.0")
  expect_identical(object_schema$schema_version, "1.6.0")
})
