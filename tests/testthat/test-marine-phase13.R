phase13_fixture <- local({
  dest <- NULL
  function() {
    if (is.null(dest)) {
      dest <<- tempfile(); dir.create(dest)
      aliases <- c(phase13_manifest = "manifest.tsv", phase13_main = "main.nwk",
        phase13_support = "support.nwk", phase13_branches = "branches.tsv",
        drop_marine = "marine.tsv", drop_aquatic = "aquatic.tsv")
      spec <- LaTerra:::.lt_marine_m2_authority()$inputs[names(aliases)]
      for (id in names(aliases)) {
        to <- file.path(dest, spec[[id]]$path)
        dir.create(dirname(to), recursive = TRUE, showWarnings = FALSE)
        stopifnot(file.copy(test_path("fixtures", "p13", aliases[[id]]), to))
      }
    }
    dest
  }
})

test_that("005 canonical materialization preserves all six full vectors and ordinals", {
  skip_if_not_installed("ape"); skip_if_not_installed("castor")
  contract <- LaTerra:::.lt_marine_phase13_contract()
  expect_identical(contract$ordinal, c(2L, 3L, 4L, 6L, 7L, 8L))
  for (id in contract$run_id) {
    x <- LaTerra:::.lt_marine_phase13_materialize(phase13_fixture(), id)
    expect_identical(names(x$state), paste0("B", 1:601))
    expect_identical(LaTerra:::.lt_marine_phase13_state_hash(x$state),
      contract$literal_sha256[contract$run_id == id])
    expect_identical(nrow(x$terminal), 302L)
    expect_identical(sum(x$state[x$terminal$branch] != x$terminal$response),
      contract$terminal_disagreements[contract$run_id == id])
    expect_no_error(LaTerra:::.lt_marine_phase13_validate(x, phase13_fixture()))
    expect_identical(x$qualifier, "historical_phase13_screening_mapping_mismatch_preserved")
  }
  for (id in c("fix_marine_binary", "fix_aquatic_v2", "fix_drop_polar_bear", "anything"))
    expect_error(LaTerra:::.lt_marine_phase13_materialize(phase13_fixture(), id), "six admitted")
})

test_that("005 coordinated vector and metadata edits cannot self-authorize", {
  skip_if_not_installed("ape"); skip_if_not_installed("castor")
  x <- LaTerra:::.lt_marine_phase13_materialize(phase13_fixture(), "fix_drop_whale")
  bad <- x; bad$state[1] <- 1
  bad$state_sha256 <- LaTerra:::.lt_marine_phase13_state_hash(bad$state)
  expect_error(LaTerra:::.lt_marine_phase13_validate(bad, phase13_fixture()), "canonical source")
  bad <- x; bad$terminal$response[1] <- 1
  bad$terminal_sha256 <- digest::digest(bad$terminal, algo = "sha256")
  expect_error(LaTerra:::.lt_marine_phase13_validate(bad, phase13_fixture()), "canonical source")
  bad <- x; bad$ordinal <- 1L
  expect_error(LaTerra:::.lt_marine_phase13_validate(bad, phase13_fixture()), "canonical source")
  bad <- x; bad$state <- bad$state[-601L]
  expect_error(LaTerra:::.lt_marine_phase13_validate(bad, phase13_fixture()), "canonical source")
  bad <- x; bad$parser_sha256[] <- "user supplied"
  expect_error(LaTerra:::.lt_marine_phase13_validate(bad, phase13_fixture()), "canonical source")
  expect_false("profile" %in% names(formals(LaTerra:::.lt_marine_phase13_materialize)))
  expect_false("context" %in% names(formals(LaTerra:::.lt_marine_phase13_gloocv)))
})

test_that("005 modified files plus user checksum are rejected before ASR or fitting", {
  skip_if_not_installed("ape"); skip_if_not_installed("castor")
  source <- phase13_fixture()
  dest <- tempfile(); dir.create(dest)
  files <- list.files(source, recursive = TRUE)
  for (f in files) {
    dir.create(dirname(file.path(dest, f)), recursive = TRUE, showWarnings = FALSE)
    file.copy(file.path(source, f), file.path(dest, f))
  }
  a <- LaTerra:::.lt_marine_m2_authority()
  file <- file.path(dest, a$inputs$phase13_main$path)
  writeLines(c(readLines(file), " "), file)
  a$inputs$phase13_main$sha256 <- digest::digest(file = file, algo = "sha256")
  expect_error(LaTerra:::.lt_marine_phase13_materialize(dest, "fix_drop_whale"),
    class = "lt_error_marine_authority_mismatch")
  expect_false(identical(LaTerra:::.lt_marine_m2_authority(), a))
})

test_that("005 literal mismatch never enables ordinary nested_phase13 calls", {
  skip_if_not_installed("ape"); skip_if_not_installed("castor")
  x <- LaTerra:::.lt_marine_phase13_materialize(phase13_fixture(), "fix_drop_whale")
  gbi <- matrix(1, 2L, 601L, dimnames = list(c("user-gene-1", "user-gene-2"), names(x$state)))
  folds <- LaTerra:::.lt_marine_model_folds(x$terminal, "nested_phase13", 2L)
  expect_error(lt_s2_gloocv(gbi, x$state, x$terminal, folds, recipe = "nested_phase13"), "disagree")
  expect_false(any(grepl("skip|historical|override|authority", names(formals(lt_s2_gloocv)))))
})
