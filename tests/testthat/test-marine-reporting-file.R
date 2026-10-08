test_that("workbook finalization preserves every logical value and resolves package parts", {
  skip_if_not_installed("openxlsx"); skip_if_not_installed("xml2"); skip_if_not_installed("zip")
  tmp <- tempfile("lt_xlsx_test_"); dir.create(tmp); on.exit(unlink(tmp, recursive = TRUE))
  before <- file.path(tmp, "before.xlsx"); after <- file.path(tmp, "after.xlsx")
  wb <- openxlsx::createWorkbook(); openxlsx::addWorksheet(wb, "typed")
  d <- data.frame(key = c("g1", "g2", "g3"), value = c(0, NA_real_, -.035),
    flag = c(TRUE, FALSE, NA), label = c("normal", "", "mixed signs"))
  openxlsx::writeData(wb, "typed", d); openxlsx::saveWorkbook(wb, before)
  prov <- LaTerra:::.lt_marine_finalize_table_file(before, after)
  expect_true(prov$non_relationship_parts_byte_identical)
  expect_false(prov$cell_values_modified)
  expect_identical(openxlsx::read.xlsx(before, skipEmptyCols = FALSE, skipEmptyRows = FALSE),
    openxlsx::read.xlsx(after, skipEmptyCols = FALSE, skipEmptyRows = FALSE))
  after2 <- file.path(tmp, "after2.xlsx")
  again <- LaTerra:::.lt_marine_finalize_table_file(after, after2)
  expect_length(again$removed_unused_relationships, 0L)
  expect_error(LaTerra:::.lt_marine_finalize_table_file(before, after), class = "lt_error_marine_output_exists")
  # An actually referenced missing drawing must fail, not silently disappear.
  parts <- file.path(tmp, "parts"); dir.create(parts); utils::unzip(before, exdir = parts)
  rel <- file.path(parts, "xl/worksheets/_rels/sheet1.xml.rels")
  doc <- xml2::read_xml(rel)
  link <- xml2::xml_find_first(doc, "//*[local-name()='Relationship' and contains(@Type, '/drawing')]")
  if (!inherits(link, "xml_missing")) {
    sheet <- file.path(parts, "xl/worksheets/sheet1.xml")
    text <- paste(readLines(sheet, warn = FALSE), collapse = "")
    text <- sub("</worksheet>", paste0('<drawing r:id="', xml2::xml_attr(link, "Id"), '"/></worksheet>'), text, fixed = TRUE)
    writeLines(text, sheet)
    broken <- file.path(tmp, "referenced.xlsx")
    zip::zipr(broken, list.files(parts, recursive = TRUE, all.files = TRUE), root = parts, mode = "mirror")
    expect_error(LaTerra:::.lt_marine_finalize_table_file(broken, file.path(tmp, "forbidden.xlsx")),
      class = "lt_error_marine_reporting_file")
    expect_false(file.exists(file.path(tmp, "forbidden.xlsx")))
  }
})

test_that("workbook part resolution fails closed on escaping or external local paths", {
  f <- LaTerra:::.lt_marine_table_part_path
  expect_identical(f("xl/worksheets", "../drawings/drawing1.xml"), "xl/drawings/drawing1.xml")
  expect_identical(f("xl/worksheets", "/xl/styles.xml"), "xl/styles.xml")
  for (p in c("../../../outside", "https://example.com/x", "../x\\y", "x#fragment", ""))
    expect_error(f("xl/worksheets", p), class = "lt_error_marine_reporting_file")
})

test_that("long annotation layout expands rows without altering logical contents", {
  v <- matrix(list("short", paste(rep("long annotation", 20), collapse = " "), NULL,
    "first\nsecond", 0, TRUE), nrow = 3)
  before <- serialize(v, NULL, version = 3)
  h <- LaTerra:::.lt_marine_table_row_heights(v, c(25, 20))
  expect_length(h, 3L)
  expect_true(h[2] > h[1])
  expect_true(all(h >= 18 & h <= 409.5))
  expect_identical(serialize(v, NULL, version = 3), before)
})
