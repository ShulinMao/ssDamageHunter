test_that("the haplotype is chosen independently at each scale", {
  # ab = 0.6, af = 0.30, so the two candidate expectations are:
  #   double-stranded  0.60 vs 0.40  -> 0.40 is closer to af
  #   single-stranded  0.30 vs 0.20  -> 0.30 is closer to af
  #   error            0.15 vs 0.10  -> 0.15 is closer to af
  # The double-stranded hypothesis therefore picks the other haplotype from the
  # other two. Choosing one haplotype for all three would be wrong.
  got <- ss_damage_likelihoods(af = 0.30, ab = 0.6, dp = 40, scalt = 12)

  expect_equal(got$double_strand_likelihood, dbinom(12, 40, 1.00 * 0.4))
  expect_equal(got$single_strand_likelihood, dbinom(12, 40, 0.50 * 0.6))
  expect_equal(got$error_likelihood,         dbinom(12, 40, 0.25 * 0.6))
  expect_equal(got$af_ab_diff, 0.30 - 0.30)
})

test_that("the haplotype choice flips to 1 - ab when af is closer to it", {
  # ab = 0.8, af = 0.10. Single-stranded expectations are 0.40 and 0.10, so the
  # 1 - ab haplotype wins.
  got <- ss_damage_likelihoods(af = 0.10, ab = 0.8, dp = 50, scalt = 5)

  expect_equal(got$single_strand_likelihood, dbinom(5, 50, 0.5 * 0.2))
  expect_equal(got$error_likelihood,         dbinom(5, 50, 0.25 * 0.2))
  expect_equal(got$af_ab_diff, 0)

  # The choice is made independently at each scale: for the double-stranded
  # hypothesis the expectations are 0.8 and 0.2, and 0.2 is still closer.
  expect_equal(got$double_strand_likelihood, dbinom(5, 50, 0.2))
})

test_that("af_ab_diff is signed", {
  balanced <- ss_damage_likelihoods(af = 0.25, ab = 0.5, dp = 40, scalt = 10)
  expect_equal(balanced$af_ab_diff, 0)

  above <- ss_damage_likelihoods(af = 0.30, ab = 0.5, dp = 40, scalt = 12)
  expect_gt(above$af_ab_diff, 0)

  below <- ss_damage_likelihoods(af = 0.20, ab = 0.5, dp = 40, scalt = 8)
  expect_lt(below$af_ab_diff, 0)
})

test_that("ss_damage_likelihoods is vectorised", {
  got <- ss_damage_likelihoods(af    = c(0.30, 0.10),
                               ab    = c(0.60, 0.80),
                               dp    = c(40, 50),
                               scalt = c(12, 5))

  expect_equal(nrow(got), 2)
  expect_equal(got$single_strand_likelihood,
               c(dbinom(12, 40, 0.30), dbinom(5, 50, 0.10)))
})

test_that("a single-stranded event is most likely under the single-stranded model", {
  # Depth 60, allele balance 0.5. A single-stranded damage event is expected at
  # af = 0.25, i.e. 15 alt reads.
  got <- ss_damage_likelihoods(af = 15 / 60, ab = 0.5, dp = 60, scalt = 15)

  expect_gt(got$single_strand_likelihood, got$double_strand_likelihood)
  expect_gt(got$single_strand_likelihood, got$error_likelihood)
})

test_that("a heterozygous germline-like site favours the double-stranded model", {
  # 30 of 60 reads at allele balance 0.5 is the double-stranded expectation.
  got <- ss_damage_likelihoods(af = 30 / 60, ab = 0.5, dp = 60, scalt = 30)

  expect_gt(got$double_strand_likelihood, got$single_strand_likelihood)
})
