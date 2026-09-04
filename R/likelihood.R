#' Likelihoods of the three competing hypotheses for a somatic candidate
#'
#' For each candidate site, three likelihoods are computed to distinguish the
#' origin of the alternate allele: a true double-stranded somatic mutation, a
#' single-stranded DNA damage event, or an amplification/sequencing error.
#'
#' Let \code{af} be the observed alternate allele frequency in the single cell,
#' \code{ab} the allele balance of the locus estimated by SCAN2, \code{dp} the
#' total read depth and \code{scalt} the number of alternate-allele supporting
#' reads.
#'
#' Under the double-stranded hypothesis (a mutation present on both strands of
#' a DNA molecule) the expected alternate allele frequency is \code{ab} or
#' \code{1 - ab}, and
#' \deqn{L_{ds} = Binom(scalt; dp, f_{ds})}
#' where \eqn{f_{ds} = ab} if \eqn{|af - ab| < |af - (1 - ab)|} and
#' \eqn{f_{ds} = 1 - ab} otherwise.
#'
#' Under the single-stranded hypothesis (damage present on only one strand,
#' which upon PTA amplification is represented at approximately half the
#' expected double-stranded frequency) the expected frequency is
#' \eqn{0.5 \times ab} or \eqn{0.5 \times (1 - ab)}, and \eqn{L_{ss}} is
#' computed the same way. The error hypothesis uses a factor of 0.25.
#'
#' @param af Numeric vector. Observed alternate allele frequency in the single
#'   cell (\code{scalt / dp}).
#' @param ab Numeric vector. Allele balance of the locus estimated by SCAN2.
#' @param dp Integer vector. Total read depth in the single cell.
#' @param scalt Integer vector. Alternate-allele supporting reads in the single
#'   cell.
#'
#' @return A \code{data.frame} with one row per input site and the columns
#'   \code{double_strand_likelihood}, \code{single_strand_likelihood},
#'   \code{error_likelihood} and \code{af_ab_diff}. The last is the signed
#'   difference between the observed \code{af} and the single-stranded
#'   expectation, used by \code{\link{single_strand_damage_filter}}.
#'
#' @examples
#' # A single site sitting exactly at the single-stranded expectation.
#' ss_damage_likelihoods(af = 0.25, ab = 0.5, dp = 40, scalt = 10)
#'
#' # The shipped synthetic example, covering all three hypotheses.
#' candidates <- read.csv(system.file("extdata", "example_candidates.csv",
#'                                    package = "ssDamageHunter"))
#' head(ss_damage_likelihoods(candidates$af, candidates$ab,
#'                            candidates$dp, candidates$scalt))
#'
#' @export
ss_damage_likelihoods <- function(af, ab, dp, scalt) {
  # Of the two haplotypes (frequency `ab` and `1 - ab`), attribute the
  # candidate to whichever one the observed `af` sits closer to.
  expected_af <- function(scale) {
    ifelse(abs(af - scale * ab) < abs(af - scale * (1 - ab)),
           scale * ab,
           scale * (1 - ab))
  }

  data.frame(
    double_strand_likelihood = stats::dbinom(scalt, size = dp, prob = expected_af(1.00)),
    single_strand_likelihood = stats::dbinom(scalt, size = dp, prob = expected_af(0.50)),
    error_likelihood         = stats::dbinom(scalt, size = dp, prob = expected_af(0.25)),
    af_ab_diff               = af - expected_af(0.50)
  )
}
