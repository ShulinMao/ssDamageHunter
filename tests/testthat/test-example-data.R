# The shipped example backs the README quickstart and the roxygen examples.
# These tests keep the documented behaviour from silently drifting.

example_candidates <- function() {
  path <- system.file("extdata", "example_candidates.csv", package = "ssDamageHunter")
  expect_true(nzchar(path))
  read.csv(path)
}

test_that("the example ships with the columns the pipeline needs", {
  d <- example_candidates()
  expect_true(all(c("cell_id", "chr", "pos", "refnt", "altnt", "muttype",
                    "mutsig", "dp", "scalt", "af", "ab", "pass") %in% colnames(d)))
  expect_gt(nrow(d), 0)
  # af must be consistent with the read counts it was derived from.
  expect_equal(d$af, round(d$scalt / d$dp, 6))
})

test_that("the example flows through the documented pipeline", {
  d <- example_candidates()
  d <- cbind(d, ss_damage_likelihoods(d$af, d$ab, d$dp, d$scalt))
  kept <- single_strand_damage_filter(d)

  expect_gt(nrow(kept), 0)
  expect_lt(nrow(kept), nrow(d))
  # Every cell and both mutation types are represented, so the example
  # exercises the per-cell and per-muttype error budgets.
  expect_setequal(as.character(unique(kept$cell_id)), c("cell_A", "cell_B", "cell_C"))
  expect_setequal(as.character(unique(kept$muttype)), c("snv", "indel"))
})

test_that("the example covers all three hypotheses the model discriminates", {
  d <- example_candidates()
  d <- cbind(d, ss_damage_likelihoods(d$af, d$ab, d$dp, d$scalt))
  best <- c("ds", "ss", "err")[apply(
    d[, c("double_strand_likelihood", "single_strand_likelihood", "error_likelihood")],
    1, which.max)]
  expect_setequal(unique(best), c("ds", "ss", "err"))
})

test_that("sites SCAN2 already called are excluded despite looking single-stranded", {
  d <- example_candidates()
  d <- cbind(d, ss_damage_likelihoods(d$af, d$ab, d$dp, d$scalt))

  scan2_called <- d[d$pass, ]
  expect_gt(nrow(scan2_called), 0)
  # They would otherwise qualify: the single-stranded model fits them best.
  expect_true(all(scan2_called$single_strand_likelihood >
                    scan2_called$double_strand_likelihood))
  expect_equal(nrow(single_strand_damage_filter(d[d$pass, ])), 0)
})

test_that("the example sample sheet documents the expected input columns", {
  path <- system.file("extdata", "example_sample_sheet.csv", package = "ssDamageHunter")
  expect_true(nzchar(path))
  expect_equal(colnames(read.csv(path)), c("cell_id", "case_id"))
})
