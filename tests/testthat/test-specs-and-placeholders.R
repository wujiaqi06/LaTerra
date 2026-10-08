test_that("rate and analysis specifications are declarative containers", {
  rate <- lt_rate_spec(
    "legacy_declared_only",
    "not_implemented",
    parameters = list(trim_quantile = 0.975)
  )
  spec <- lt_analysis_spec(
    "skeleton",
    rate,
    branch_state = list(method = "declared_only"),
    screen = list(method = "declared_only"),
    null_model = list(list(method = "declared_only"))
  )

  expect_s3_class(rate, "lt_rate_spec")
  expect_s3_class(spec, "lt_analysis_spec")
  expect_length(spec$null_model, 1L)
  expect_output(print(spec), "scientific execution: not implemented")
})

test_that("scientific module placeholders stop explicitly", {
  expect_error(LaTerra:::lt_screen(), "not implemented")
  expect_error(LaTerra:::lt_model(), "not implemented")
  expect_error(plot(small_lt_matrix()), "not implemented")
})
