test_that("Marine profile is explicit, frozen and excludes output oracles", {
  p <- lt_marine_profile()
  expect_identical(p$profile_id, "marine_S2_frozen")
  expect_identical(p$inputs$branch_archive$sha256,
    "1dcc4bdbec568c3e02a17077a528ecc40ecfe1598c98de23195d215b4e046552")
  expect_false(p$execution_reads_expected_outputs)
  expect_identical(names(p$inputs), c("branch_archive", "traits", "complete_state"))
  expect_identical(p$authority_addendum$id, "TARGET001-AUTHORITY-ADDENDUM001")
  expect_false(p$authority_addendum$changes_numerical_science)
  expect_false(any(grepl("gbi|screen|expected", unlist(p$inputs))))
  expect_identical(p$traits$aquatic_v2, "aquatic_v2_status")
})

test_that("raw token parser keeps zeros distinct from every upstream absence", {
  f <- tempfile(fileext = ".tsv")
  on.exit(unlink(f))
  x <- data.frame(gene = paste0("g", 1:11), B1 = c("0", "1.25", "NA", "NA_fuse",
    "NA_struct", "NA_topo", "residual_NA", "NaN", "Inf", "-Inf", ""))
  utils::write.table(x, f, sep = "\t", row.names = FALSE, quote = FALSE)
  a <- LaTerra:::.lt_marine_parse_matrix(f)
  expect_identical(a$values[1:2, 1L], c(g1 = 0, g2 = 1.25))
  expect_true(all(is.na(a$values[3:11, 1L])))
  expect_identical(a$token_vocabulary[as.integer(a$token_code) + 1L],
    c("finite_numeric", "finite_numeric", "NA", "NA_fuse", "NA_struct", "NA_topo",
      "residual_NA", "NaN", "Inf", "-Inf", "empty"))
  x$B1[2] <- "surprise"
  utils::write.table(x, f, sep = "\t", row.names = FALSE, quote = FALSE)
  expect_error(LaTerra:::.lt_marine_parse_matrix(f), class = "lt_error_marine_matrix_token")
  x$B1[2] <- "1"; x$gene[2] <- x$gene[1]
  utils::write.table(x, f, sep = "\t", row.names = FALSE, quote = FALSE)
  expect_error(LaTerra:::.lt_marine_parse_matrix(f), "duplicat|unique")
})

test_that("approved complete state and historical partial provenance stay separate", {
  f <- tempfile(); p <- tempfile(); s <- tempfile()
  on.exit(unlink(c(f, p, s)))
  writeLines(c("gene\tB1\tB2", "g1\t0\tNA", "g2\tNA\t2"), f)
  writeLines(c("gene\tB1\tB2", "g1\t0\tNA_struct"), p)
  writeLines(c("gene\tB1\tB2", "g1\tobserved\tNA_struct", "g2\tNA_fuse\tobserved"), s)
  raw <- LaTerra:::.lt_marine_parse_matrix(f)
  hist <- LaTerra:::.lt_marine_parse_matrix(p)
  full <- LaTerra:::.lt_marine_parse_complete_state(s)
  before <- serialize(raw, NULL, version = 3)
  out <- LaTerra:::.lt_marine_state_layers(raw, full, hist)
  expect_identical(serialize(raw, NULL, version = 3), before)
  expect_identical(out$historical_partial$gene_ids, "g1")
  expect_identical(out$historical_partial$token_code, hist$token_code)
  expect_identical(out$complete_state$state_code, full$state_code)
  expect_identical(dim(out$complete_state$state_code), c(2L, 2L))
  expect_identical(out$complete_state$vocabulary[as.integer(full$state_code) + 1L],
    c("observed", "NA_fuse", "NA_struct", "observed"))
  wrong <- full; wrong$state_code <- wrong$state_code[2:1, , drop = FALSE]
  expect_error(LaTerra:::.lt_marine_state_layers(raw, wrong, hist), class = "lt_error_marine_axis")
  wrong <- full; wrong$state_code[1, 2] <- as.raw(0)
  expect_error(LaTerra:::.lt_marine_state_layers(raw, wrong, hist), class = "lt_error_marine_state_payload")
  wrong <- hist; wrong$gene_ids <- "missing"
  expect_error(LaTerra:::.lt_marine_state_layers(raw, full, wrong), class = "lt_error_marine_axis")
  writeLines(c("gene\tB1", "g1\tguessed_state"), s)
  expect_error(LaTerra:::.lt_marine_parse_complete_state(s), class = "lt_error_marine_state_token")
})

test_that("S2 baseline reproduces inclusive trimming and exact source operation order", {
  y <- rbind(g1 = c(1:17, 40, 40, 40), g2 = rep(0, 20),
    g3 = c(rep(NA_real_, 19), 3), g4 = rep(NA_real_, 20), g5 = 2:21)
  colnames(y) <- paste0("B", 1:20)
  raw <- list(values = y, gene_ids = rownames(y), branch_ids = colnames(y))
  b <- LaTerra:::.lt_marine_baseline(raw)
  manual <- y
  for (i in seq_len(nrow(y))) {
    if (!all(is.na(y[i, ]))) {
      q <- quantile(y[i, ], .975, na.rm = TRUE, names = FALSE, type = 7)
      manual[i, which(y[i, ] >= q)] <- NA_real_
    }
  }
  avg <- function(x) if (all(is.na(x))) NA_real_ else mean(x[!is.na(x)])
  ge <- apply(manual, 1L, avg); be <- apply(manual, 2L, avg)
  expected <- sweep(sweep(manual, 2L, be, "/"), 1L, ge, "/")
  expected[!is.finite(expected)] <- NA_real_
  expect_identical(b$trimmed_values, manual)
  expect_identical(b$gbi, expected)
  expect_equal(sum(b$trim_mask["g1", ]), 3L)
  expect_equal(sum(b$trim_mask["g2", ]), 20L)
  expect_true(all(is.na(b$gbi["g2", ])))
  expect_true(is.na(b$gene_effect$GE_mean[4L]))
})

test_that("historical strict-prefix FDR is not silently replaced by standard BH", {
  f <- LaTerra:::.lt_marine_historical_fdr
  expect_identical(f(c(.006, .006)), integer())
  expect_true(all(p.adjust(c(.006, .006), "BH") < .01))
  expect_identical(f(c(.005, .009)), integer()) # equality at first threshold
  expect_identical(f(c(.001, .02, .021)), 1L)
  expect_identical(f(c(.0001, .001)), 1:2)
  expect_identical(f(numeric()), integer())
})

test_that("Welch screen uses gene-wise finite values and explicit untested reasons", {
  x <- rbind(g1 = c(2, 3, 4, 5, 1, 1.2, 1.8, 2.1),
             g2 = c(NA, NA, NA, 4, 1, 2, 3, 4),
             g3 = c(1, 2, 3, 4, 4, 6, 7, 8))
  colnames(x) <- paste0("B", 1:8)
  s <- data.frame(branch_id = colnames(x), state = rep(c(0, 1), each = 4),
    screened = TRUE, run_id = "fix_marine_binary", trait = "marine_binary")
  a <- LaTerra:::.lt_marine_screen(x, s)
  expect_equal(nrow(a$tested), 2L)
  expect_identical(a$gene_ledger$reason[2], "insufficient_finite_group_size")
  direct <- t.test(x[1, 1:4], x[1, 5:8])
  j <- match("g1", a$tested$gene)
  expect_identical(a$tested$pvalue[j], direct$p.value)
  expect_equal(a$tested$tvalue[j], unname(direct$statistic))
  expect_identical(a$gene_ledger$direction, c("slow", NA_character_, "fast"))
  expect_identical(a$gene_ledger$n_focal, c(4L, 4L, 4L))
  expect_identical(a$gene_ledger$n_reference, c(4L, 1L, 4L))
  s$branch_id <- rev(s$branch_id)
  expect_error(LaTerra:::.lt_marine_screen(x, s), class = "lt_error_marine_axis")
})

test_that("tree annotation joins supplied split identities rather than branch position", {
  skip_if_not_installed("ape")
  skip_if_not_installed("castor")
  tree <- ape::read.tree(text = "((a,b),c,d);")
  keys <- data.frame(branch_label = paste0("old", 1:5),
    side_A_taxa = c("a", "b", "c", "d", "a;b"),
    side_B_taxa = c("b;c;d", "a;c;d", "a;b;d", "a;b;c", "c;d"),
    canonical_split_key = paste0("key", 1:5),
    branch_type = c(rep("terminal", 4), "internal"),
    terminal_taxon_if_terminal = c("a", "b", "c", "d", NA))
  join <- LaTerra:::.lt_marine_branch_join(tree, keys, keys$branch_label)
  expect_identical(tree$tip.label[join$offspring[1:4]], c("a", "b", "c", "d"))
  traits <- data.frame(species = c("d", "c", "b", "a"), marine_binary = c(0, 0, 1, 1))
  before <- if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) get(".Random.seed", .GlobalEnv) else NULL
  result <- LaTerra:::.lt_marine_annotate(tree, traits, join, "marine_binary", "marine_binary")
  after <- if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) get(".Random.seed", .GlobalEnv) else NULL
  expect_identical(before, after)
  # Explicitly exercise both cold, absent RNG and pre-existing RNG streams.
  had_seed <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
  saved_seed <- if (had_seed) get(".Random.seed", .GlobalEnv) else NULL
  on.exit({
    if (had_seed) assign(".Random.seed", saved_seed, .GlobalEnv) else
      if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  if (had_seed) rm(".Random.seed", envir = .GlobalEnv)
  cold <- LaTerra:::.lt_marine_annotate(tree, traits, join, "marine_binary", "marine_binary")
  expect_false(exists(".Random.seed", .GlobalEnv, inherits = FALSE))
  for (seed in c(1L, 20260520L, 99L)) {
    set.seed(seed); before_seeded <- get(".Random.seed", .GlobalEnv)
    seeded <- LaTerra:::.lt_marine_annotate(tree, traits, join, "marine_binary", "marine_binary")
    expect_identical(get(".Random.seed", .GlobalEnv), before_seeded)
    expect_identical(seeded, cold)
  }
  expect_equal(result$state[1:4], c(1, 1, 0, 0))
  expect_identical(result$branch_id, keys$branch_label)
  traits$species[1] <- "missing"
  expect_error(LaTerra:::.lt_marine_annotate(tree, traits, join, "marine_binary", "marine_binary"),
    class = "lt_error_marine_trait")
})

test_that("public entry rejects overwrite and missing data without hidden fallback", {
  expect_error(lt_run_marine(tempfile(), tempfile()), class = "lt_error_marine_missing_input")
  root <- tempfile(); dir.create(root); on.exit(unlink(root, recursive = TRUE))
  expect_error(lt_run_marine(root, root), class = "lt_error_marine_output_exists")
  skip_if_not_installed("ape"); skip_if_not_installed("castor")
  p <- lt_marine_profile()
  expect_error(LaTerra:::.lt_marine_inputs(root, p, root), class = "lt_error_marine_missing_input")
  f <- file.path(root, "wrong.tsv"); writeLines("wrong input", f)
  expect_error(LaTerra:::.lt_marine_verify(f, paste(rep("0", 64), collapse = ""), "traits"),
    class = "lt_error_marine_authority_mismatch")
})
