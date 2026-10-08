.lt_abort <- function(message) {
  stop(message, call. = FALSE)
}

.lt_abort_class <- function(message, subclass, fields = list()) {
  .lt_assert_scalar_character(subclass, "subclass")
  .lt_assert_named_list(fields, "fields")
  condition <- structure(
    c(list(message = message, call = NULL), fields),
    class = c(subclass, "lt_error", "error", "condition")
  )
  stop(condition)
}

.lt_assert_scalar_character <- function(x, name, allow_na = FALSE) {
  ok <- is.character(x) && length(x) == 1L
  if (!allow_na) {
    ok <- ok && !is.na(x) && nzchar(x)
  }
  if (!ok) {
    .lt_abort(sprintf("`%s` must be a non-empty character scalar.", name))
  }
  invisible(x)
}

.lt_assert_nullable_scalar_character <- function(x, name) {
  ok <- is.character(x) && length(x) == 1L &&
    (is.na(x) || nzchar(trimws(x)))
  if (!ok) {
    .lt_abort(sprintf("`%s` must be NA or a non-empty character scalar.", name))
  }
  invisible(x)
}

.lt_assert_sha256 <- function(x, name, allow_na = FALSE) {
  ok <- is.character(x) && length(x) == 1L
  if (allow_na && ok && is.na(x)) {
    return(invisible(x))
  }
  ok <- ok && !is.na(x) && grepl("^[0-9a-f]{64}$", x)
  if (!ok) {
    .lt_abort(sprintf("`%s` must be a lowercase 64-character SHA-256.", name))
  }
  invisible(x)
}

.lt_assert_named_list <- function(x, name) {
  if (!is.list(x)) {
    .lt_abort(sprintf("`%s` must be a list.", name))
  }
  if (length(x) > 0L && (is.null(names(x)) || any(!nzchar(names(x))))) {
    .lt_abort(sprintf("`%s` must be a named list.", name))
  }
  invisible(x)
}

.lt_assert_unique_ids <- function(x, name) {
  if (!is.character(x) || anyNA(x) || any(!nzchar(x))) {
    .lt_abort(sprintf("`%s` must contain non-empty character identifiers.", name))
  }
  duplicates <- unique(x[duplicated(x)])
  if (length(duplicates) > 0L) {
    .lt_abort(sprintf(
      "`%s` contains duplicate identifiers: %s.",
      name,
      paste(utils::head(duplicates, 5L), collapse = ", ")
    ))
  }
  invisible(x)
}

.lt_validate_component_list <- function(x, name) {
  .lt_assert_named_list(x, name)
  invisible(x)
}

.lt_assert_data_frame_columns <- function(x, columns, name) {
  if (!is.data.frame(x) || !identical(names(x), columns)) {
    .lt_abort(sprintf(
      "`%s` must be a data frame with columns: %s.",
      name,
      paste(columns, collapse = ", ")
    ))
  }
  invisible(x)
}

.lt_not_implemented <- function(component) {
  .lt_abort(sprintf("%s: not implemented", component))
}

.lt_count_label <- function(n, singular, plural = paste0(singular, "s")) {
  paste(n, if (identical(n, 1L)) singular else plural)
}

.lt_hash <- function(x) {
  digest::digest(x, algo = "sha256", serialize = TRUE)
}

# Boundary-object identities freeze the R serialization version explicitly.
# This is intentionally separate from .lt_hash(): changing the historical
# package-wide hash contract would invalidate pre-existing scientific objects.
.lt_hash_serialization_v3 <- function(x) {
  bytes <- serialize(x, NULL, version = 3)
  digest::digest(bytes, algo = "sha256", serialize = FALSE)
}

.lt_hex <- function(x) {
  ifelse(is.na(x), "NA", sprintf("%a", x))
}

# Deterministic compensated pairwise summation. Each merge uses an error-free
# TwoSum transform and carries the residual through a frozen adjacent-pair tree.
.lt_pairwise_sum <- function(x) {
  x <- as.double(x)
  n <- length(x)
  if (n == 0L) return(0)
  if (any(!is.finite(x))) return(sum(x))
  s <- x
  c <- numeric(n)
  while (length(s) > 1L) {
    n <- length(s)
    m <- n %/% 2L
    odd <- n %% 2L == 1L
    ia <- seq.int(1L, 2L * m, by = 2L)
    ib <- ia + 1L
    a <- s[ia]
    b <- s[ib]
    z <- a + b
    bp <- z - a
    e <- (a - (z - bp)) + (b - bp)
    carry <- c[ia] + c[ib] + e
    z2 <- z + carry
    bp2 <- z2 - z
    e2 <- (z - (z2 - bp2)) + (carry - bp2)
    if (odd) {
      s <- c(z2, s[n])
      c <- c(e2, c[n])
    } else {
      s <- z2
      c <- e2
    }
  }
  s[[1L]] + c[[1L]]
}
