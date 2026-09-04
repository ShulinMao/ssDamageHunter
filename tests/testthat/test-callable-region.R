test_that("effective callable region weights each depth by its detection probability", {
  callable <- c(100, 200, 300)
  detect   <- c(0.0, 0.5, 1.0)
  expect_equal(effective_callable_region_size(callable, detect), 0 + 100 + 300)
})

test_that("mismatched lengths are an error, not silent recycling", {
  expect_error(effective_callable_region_size(c(1, 2, 3), c(0.5, 0.5)),
               "same length")
})

test_that("is_single_strand_call accepts the single-stranded expectation", {
  # 15 alt of 60 reads at allele balance 0.5 is exactly 0.25.
  expect_true(is_single_strand_call(alt_dp = 15, dp = 60, ab = 0.5))
})

test_that("is_single_strand_call rejects the heterozygous expectation", {
  expect_false(is_single_strand_call(alt_dp = 30, dp = 60, ab = 0.5))
})

test_that("is_single_strand_call is vectorised over allele balance", {
  got <- is_single_strand_call(alt_dp = 15, dp = 60, ab = c(0.5, 0.9))
  expect_length(got, 2)
  expect_true(got[1])
})
