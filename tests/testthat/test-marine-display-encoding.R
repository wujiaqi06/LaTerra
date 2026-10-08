test_that("Fig5A interval separator retains UTF-8 bytes in the C locale", {
  old <- Sys.getlocale("LC_CTYPE")
  on.exit(Sys.setlocale("LC_CTYPE", old), add = TRUE)
  suppressWarnings(Sys.setlocale("LC_CTYPE", "C"))
  fold <- list(identity=list(test_species=c("a","b"),y_test=c(0,1)),
    prediction=c(.1,.9),beta=c(g1=1),status="fit")
  ids <- LaTerra:::.lt_marine_sensitivity_ids()
  results <- stats::setNames(rep(list(list(folds=list(fold,fold))),length(ids)),ids)
  table <- LaTerra:::.lt_marine_sensitivity_tables(results)$table
  text <- table$median_selected_predictors_per_fold_IQR[1L]
  expect_identical(charToRaw(text),as.raw(c(0x31,0x20,0x28,0x31,0xe2,0x80,0x93,0x31,0x29)))
  p <- tempfile(); on.exit(unlink(p),add=TRUE)
  LaTerra:::.lt_marine_write_tsv(data.frame(interval=text),p)
  expect_identical(charToRaw(readLines(p)[2L]),charToRaw(text))
  expect_false(grepl("<U+",readLines(p)[2L],fixed=TRUE))
})
