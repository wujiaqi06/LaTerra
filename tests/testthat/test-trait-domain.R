test_that("trait and domain constructors validate aligned objects", {
  trait <- lt_trait(
    "marine_binary",
    c("Taxon_a", "Taxon_b"),
    c(1, 0),
    coding = list(positive = 1, background = 0)
  )
  domain <- lt_domain(
    "all_taxa",
    c("Taxon_b", "Taxon_a"),
    group = c("Taxon", "Taxon")
  )

  expect_s3_class(trait, "lt_trait")
  expect_s3_class(domain, "lt_domain")
  expect_true(validate_lt_trait_domain(trait, domain))
})

test_that("trait/taxon mismatch is detected", {
  trait <- lt_trait("marine_binary", c("Taxon_a", "Taxon_b"), c(1, 0))
  domain <- lt_domain("partial", c("Taxon_a", "Taxon_c"))

  expect_error(
    validate_lt_trait_domain(trait, domain),
    "Trait/taxon mismatch detected"
  )
})

test_that("ineligible taxa require explicit exclusion reasons", {
  domain <- lt_domain(
    "filtered",
    c("Taxon_a", "Taxon_b"),
    eligible = c(TRUE, FALSE),
    exclusion_reason = c(NA_character_, "intermediate trait state")
  )
  expect_invisible(validate_lt_domain(domain))

  expect_error(
    lt_domain(
      "unnamed_exclusion",
      c("Taxon_a", "Taxon_b"),
      eligible = c(TRUE, FALSE)
    ),
    "requires a non-empty"
  )
})

test_that("eligible taxa cannot carry exclusion reasons", {
  expect_error(
    lt_domain(
      "contradictory",
      c("Taxon_a", "Taxon_b"),
      eligible = c(TRUE, TRUE),
      exclusion_reason = c(NA_character_, "should not be present")
    ),
    "Eligible taxa must have an NA"
  )
})
