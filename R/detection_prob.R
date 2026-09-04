#' Would a site with this depth and allele balance be called single-stranded?
#'
#' The decision rule of \code{\link{single_strand_damage_filter}}, evaluated on
#' hypothetical rather than observed data. Used to compute detection
#' sensitivity by sweeping over depths and empirical allele balances.
#'
#' Two approximations distinguish this from the production filter, both
#' appropriate to a sensitivity simulation:
#' \itemize{
#'   \item the expected frequencies use \code{ab} directly rather than choosing
#'     between \code{ab} and \code{1 - ab}, since the simulation draws \code{ab}
#'     from the empirical distribution and the two branches are symmetric;
#'   \item the per-cell cumulative error budget is replaced by the per-site
#'     proxy \code{error_likelihood < max_error_likelihood}.
#' }
#'
#' @param alt_dp Number of alternate-allele supporting reads.
#' @param dp Total read depth.
#' @param ab Allele balance (vector; the rule is evaluated for each element).
#' @param max_af_ab_diff,min_single_strand_likelihood,max_error_likelihood
#'   Thresholds, matching the defaults of
#'   \code{\link{single_strand_damage_filter}}.
#'
#' @return Logical vector, the same length as \code{ab}.
#'
#' @export
is_single_strand_call <- function(alt_dp, dp, ab,
                                  max_af_ab_diff = 0.1,
                                  min_single_strand_likelihood = 0.05,
                                  max_error_likelihood = 0.005) {
  double_strand_likelihood <- stats::dbinom(alt_dp, dp, ab)
  single_strand_likelihood <- stats::dbinom(alt_dp, dp, 0.5 * ab)
  error_likelihood         <- stats::dbinom(alt_dp, dp, 0.25 * ab)

  af_ab_diff <- alt_dp / dp - 0.5 * ab

  abs(af_ab_diff) < max_af_ab_diff &
    single_strand_likelihood > min_single_strand_likelihood &
    single_strand_likelihood > double_strand_likelihood &
    single_strand_likelihood > error_likelihood &
    error_likelihood < max_error_likelihood
}

#' Probability of detecting single-stranded damage, as a function of depth
#'
#' Marginalises the decision rule over the alternate-read count and the
#' empirical allele-balance distribution of the cell:
#' \deqn{P(detect | dp) = \sum_{alt} P(alt | dp) \, E_{ab}[ 1(call | alt, dp, ab) ]}
#' where \eqn{P(alt | dp) = Binom(alt; dp, 0.25)} — the expected alternate
#' allele frequency of a single-stranded damage event after PTA amplification.
#'
#' The allele-balance distribution is taken from high-confidence heterozygous
#' germline variants in the same cell (SCAN2 resampled training sites with a
#' \code{0/1} bulk genotype passing the cell's own depth and alt-read minima).
#'
#' @param results A SCAN2 object.
#' @param muttype \code{"snv"} or \code{"indel"}.
#' @param max_het_germline_count Number of heterozygous germline variants to
#'   sample for the allele-balance distribution.
#' @param ... Threshold overrides passed to \code{\link{is_single_strand_call}}.
#'
#' @return Numeric vector of detection probabilities, one per single-cell depth
#'   from the cell's \code{min.sc.dp} up to the maximum depth in its depth
#'   profile. Aligns element-wise with
#'   \code{\link{callable_region_profile}(results, muttype)}.
#'
#' @details Sampling is random; set \code{set.seed()} for reproducibility.
#'
#' @export
prob_detect_ss_damage <- function(results,
                                  muttype,
                                  max_het_germline_count = 1000,
                                  ...) {
  max_dp <- dim(results@depth.profile$dptab)[2] - 1

  min_sc_dp     <- results@static.filter.params[[muttype]]$min.sc.dp
  min_sc_alt_dp <- results@static.filter.params[[muttype]]$min.sc.alt
  min_bulk_dp   <- results@static.filter.params[[muttype]]$min.bulk.dp

  dp_distribution <- rowSums(
    results@depth.profile$dptab[(min_sc_dp + 1):(max_dp + 1), (min_bulk_dp + 1):(max_dp + 1)],
    na.rm = TRUE)

  # allelic imbalance distribution
  # high-confidence heterozygous germline variants
  het_germline_filter <-
    results@gatk$muttype == muttype &
    results@gatk$resampled.training.site &
    results@gatk$bulk.gt == "0/1" &
    results@gatk$scalt > min_sc_alt_dp &
    results@gatk$dp > min_sc_dp
  het_germline <- results@gatk[which(het_germline_filter), ]
  het_germline <- het_germline[sample(seq_len(nrow(het_germline)),
                                      max_het_germline_count, replace = TRUE), ]

  # probability of detecting single-stranded damage at each depth
  prob_detect_ssDNA_damage_all_dp <- numeric(max_dp - min_sc_dp + 1)

  idx_dp <- 1
  for (sc_dp in seq(min_sc_dp, max_dp)) {
    alt_range <- seq(min_sc_alt_dp, sc_dp)
    # Scale the number of allele-balance draws with how much of the genome sits
    # at this depth, with a floor of 50, to keep rare depths cheap.
    ab <- het_germline$ab[1:ceiling(max(
      min(max_het_germline_count, dp_distribution[as.character(sc_dp)] * 0.05), 50))]

    # p(detect ssDNA damage | alt_dp, dp)
    conditional_prob_alt_dp_dp <- vapply(
      alt_range,
      function(sc_alt_dp) mean(is_single_strand_call(sc_alt_dp, sc_dp, ab, ...)),
      numeric(1)
    )

    # p(alt_dp | dp)
    prob_alt_dp_dp <- stats::dbinom(alt_range, sc_dp, 0.25)
    # p(detect ssDNA damage | dp)
    prob_detect_ssDNA_damage_all_dp[idx_dp] <- sum(prob_alt_dp_dp * conditional_prob_alt_dp_dp)
    idx_dp <- idx_dp + 1
  }
  prob_detect_ssDNA_damage_all_dp
}
