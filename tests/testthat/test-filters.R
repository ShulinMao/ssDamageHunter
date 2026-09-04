# A candidate table with the columns single_strand_damage_filter() requires.
make_candidates <- function(...) {
  defaults <- list(
    cell_id = "cell1", muttype = "snv", pass = FALSE,
    af_ab_diff = 0, chr = "chr1", pos = 1L,
    single_strand_likelihood = 0.5,
    double_strand_likelihood = 0.01,
    error_likelihood = 0.001,
    mutsig = "A[C>T]G"
  )
  given <- list(...)
  n <- max(vapply(given, length, integer(1)), 1L)
  as.data.frame(c(
    lapply(given, function(v) rep_len(v, n)),
    lapply(defaults[setdiff(names(defaults), names(given))], rep_len, n)
  ), stringsAsFactors = FALSE)
}

test_that("each likelihood criterion can reject on its own", {
  expect_equal(nrow(single_strand_damage_filter(make_candidates())), 1)

  # allele frequency too far from the single-stranded expectation
  expect_equal(nrow(single_strand_damage_filter(make_candidates(af_ab_diff = 0.2))), 0)
  # single-stranded likelihood below the absolute floor
  expect_equal(nrow(single_strand_damage_filter(
    make_candidates(single_strand_likelihood = 0.04, error_likelihood = 1e-4))), 0)
  # a double-stranded mutation explains it better
  expect_equal(nrow(single_strand_damage_filter(
    make_candidates(double_strand_likelihood = 0.9))), 0)
  # an error explains it better
  expect_equal(nrow(single_strand_damage_filter(
    make_candidates(error_likelihood = 0.9))), 0)
  # already called as a somatic mutation by SCAN2
  expect_equal(nrow(single_strand_damage_filter(make_candidates(pass = TRUE))), 0)
})

test_that("the cumulative error budget truncates the ranked list", {
  # Four sites at 0.002 each. Cumulative sums are 0.002, 0.004, 0.006, 0.008;
  # the default budget of 0.005 admits the first two.
  candidates <- make_candidates(error_likelihood = rep(0.002, 4), pos = 1:4)
  expect_equal(nrow(single_strand_damage_filter(candidates)), 2)

  # A larger budget admits more.
  expect_equal(
    nrow(single_strand_damage_filter(candidates, cum_error_likelihood = 0.007)), 3)
})

test_that("sites are ranked by increasing error likelihood, best first", {
  candidates <- make_candidates(error_likelihood = c(0.004, 0.001, 0.0005),
                                pos = 1:3)
  kept <- single_strand_damage_filter(candidates)
  expect_equal(sort(kept$pos), c(2L, 3L))
  expect_false(1L %in% kept$pos)
})

test_that("the budget is per cell, not shared across cells", {
  candidates <- make_candidates(
    cell_id = rep(c("cell1", "cell2"), each = 4),
    error_likelihood = rep(0.002, 8),
    pos = 1:8)
  kept <- single_strand_damage_filter(candidates)
  expect_equal(as.vector(table(kept$cell_id)), c(2L, 2L))
})

test_that("SNVs and indels get independent budgets", {
  candidates <- make_candidates(
    muttype = rep(c("snv", "indel"), each = 4),
    error_likelihood = rep(0.002, 8),
    pos = 1:8)
  kept <- single_strand_damage_filter(candidates)
  expect_equal(as.vector(table(kept$muttype)), c(2L, 2L))
})

test_that("recurrent_site_filter drops sites seen in more than one cell", {
  candidates <- make_candidates(cell_id = c("cell1", "cell2"), pos = c(10L, 10L))
  kept <- recurrent_site_filter(candidates, candidates)
  expect_equal(nrow(kept), 0)

  unique_sites <- make_candidates(cell_id = c("cell1", "cell2"), pos = c(10L, 20L))
  expect_equal(nrow(recurrent_site_filter(unique_sites, unique_sites)), 2)
})

test_that("recurrent_site_filter flags cells with excess recurrent T>A calls", {
  # cell1 carries 6 recurrent T>A sites; cell2 carries one non-recurrent call.
  recurrent <- make_candidates(cell_id = "cell1", pos = 1:6, mutsig = "A[T>A]G")
  other     <- make_candidates(cell_id = "cell2", pos = 100L)
  damage    <- rbind(recurrent, other)
  # Recurrence is defined over the full candidate set: duplicate cell1's sites.
  somatic   <- rbind(damage, make_candidates(cell_id = "cellX", pos = 1:6))

  kept <- recurrent_site_filter(somatic, damage)
  expect_false(any(kept$cell_qc_filter[kept$cell_id == "cell1"]))
  expect_true(all(kept$cell_qc_filter[kept$cell_id == "cell2"]))
})

test_that("cells with recurrent sites that are not T>A are not flagged", {
  recurrent <- make_candidates(cell_id = "cell1", pos = 1:6, mutsig = "A[C>T]G")
  somatic   <- rbind(recurrent, make_candidates(cell_id = "cellX", pos = 1:6))

  kept <- recurrent_site_filter(somatic, recurrent)
  expect_true(all(kept$cell_qc_filter))
})
