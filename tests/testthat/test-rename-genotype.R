test_that("the single-cell genotype column is found by name, not position", {
  # The same sample column at two different positions, as in the 82-column and
  # 73-column SCAN2 schemas.
  wide   <- data.frame(a = 1, b = 2, `936ECA` = "0/1", c = 3, check.names = FALSE)
  narrow <- data.frame(`936ECA` = "0/1", a = 1, check.names = FALSE)

  expect_equal(colnames(rename_single_cell_genotype(wide, "936ECA"))[3], "sc_gp")
  expect_equal(colnames(rename_single_cell_genotype(narrow, "936ECA"))[1], "sc_gp")
})

test_that("a missing single-cell column is an error, not a silent rename", {
  expect_error(
    rename_single_cell_genotype(data.frame(a = 1, b = 2), "936ECA"),
    "no column named '936ECA'")
})
