#' Select single-stranded DNA damage candidates by likelihood
#'
#' Applies the single-stranded decision rule to a table of somatic candidates
#' annotated by \code{\link{call_ss_damage}}. A candidate is retained when the
#' observed allele frequency is close to the single-stranded expectation, the
#' single-stranded likelihood exceeds an absolute floor, it beats both
#' competing hypotheses, and the site was not already called as a
#' double-stranded somatic mutation by SCAN2 (\code{pass}).
#'
#' Surviving candidates are then ranked within each cell by increasing error
#' likelihood, and truncated at the point where the cumulative error likelihood
#' reaches \code{cum_error_likelihood}. This bounds the expected number of
#' false calls per cell rather than per site. SNVs and indels are ranked
#' independently.
#'
#' @param single_strand_candidates Data frame of somatic candidates with the
#'   columns added by \code{\link{call_ss_damage}}: \code{af_ab_diff},
#'   \code{single_strand_likelihood}, \code{double_strand_likelihood},
#'   \code{error_likelihood}, plus \code{pass}, \code{muttype} and
#'   \code{cell_id}.
#' @param max_af_ab_diff Maximum absolute deviation of the observed allele
#'   frequency from the single-stranded expectation.
#' @param min_single_strand_likelihood Absolute floor on the single-stranded
#'   likelihood.
#' @param cum_error_likelihood Per-cell budget for the cumulative error
#'   likelihood; the expected number of error calls per cell.
#'
#' @return The subset of \code{single_strand_candidates} passing the filter,
#'   with an added \code{cum_error} column.
#'
#' @examples
#' candidates <- read.csv(system.file("extdata", "example_candidates.csv",
#'                                    package = "ssDamageHunter"))
#' candidates <- cbind(candidates,
#'                     ss_damage_likelihoods(candidates$af, candidates$ab,
#'                                           candidates$dp, candidates$scalt))
#' single_strand_damage_filter(candidates)
#'
#' @export
single_strand_damage_filter <- function(single_strand_candidates,
                                        max_af_ab_diff = 0.1,
                                        min_single_strand_likelihood = 0.05,
                                        cum_error_likelihood = 0.005) {
  filter <-
    abs(single_strand_candidates$af_ab_diff) < max_af_ab_diff &
    single_strand_candidates$single_strand_likelihood > min_single_strand_likelihood &
    single_strand_candidates$single_strand_likelihood > single_strand_candidates$double_strand_likelihood &
    single_strand_candidates$single_strand_likelihood > single_strand_candidates$error_likelihood &
    !single_strand_candidates$pass

  single_strand_candidates$single_strand_likelihood_filter <- filter

  # Rank by error likelihood within each cell and cut off once the cumulative
  # error likelihood exceeds the per-cell budget. SNVs and indels are ranked
  # separately so that one mutation type cannot consume the other's budget.
  cumulative_error_cutoff <- function(muttype) {
    candidates <- single_strand_candidates[
      single_strand_candidates$single_strand_likelihood_filter &
        single_strand_candidates$muttype == muttype, ]
    candidates %>%
      dplyr::group_by(.data$cell_id) %>%
      dplyr::arrange(.data$error_likelihood, .by_group = TRUE) %>%
      dplyr::mutate(cum_error = cumsum(.data$error_likelihood)) %>%
      dplyr::filter(.data$cum_error < cum_error_likelihood)
  }

  rbind(cumulative_error_cutoff("snv"),
        cumulative_error_cutoff("indel"))
}

#' Remove recurrent sites and flag low-quality cells
#'
#' Sites recurring across cells are unlikely to be independent damage events
#' and are more consistent with systematic artifacts. This function (1) flags
#' cells carrying an excess of recurrent T>A calls, the signature of the
#' artifact this filter targets, and (2) drops single-stranded calls at any
#' site seen in more than one cell.
#'
#' A cell is flagged when it has at least \code{count_recurrent_threshold}
#' recurrent sites (or a fraction above \code{frac_recurrent_thresold}) *and*
#' the fraction of its calls that are both recurrent and T>A exceeds
#' \code{frac_TA_recurrent_threshold}. Flagged cells are labelled via
#' \code{cell_qc_filter} rather than removed, so that downstream steps can
#' decide how to treat them.
#'
#' @param somatic_candidates Data frame of all somatic candidates (before the
#'   single-stranded filter). Recurrence is defined over this larger set.
#' @param single_strand_candidates Data frame of called single-stranded damage.
#' @param count_recurrent_threshold Minimum number of recurrent sites in a cell.
#' @param frac_recurrent_thresold Minimum fraction of recurrent sites in a cell.
#' @param frac_TA_recurrent_threshold Minimum fraction of calls that are both
#'   recurrent and T>A.
#'
#' @return \code{single_strand_candidates} with recurrent sites removed and a
#'   logical \code{cell_qc_filter} column added (\code{FALSE} = flagged cell).
#'
#' @export
recurrent_site_filter <- function(somatic_candidates,
                                  single_strand_candidates,
                                  count_recurrent_threshold = 5,
                                  frac_recurrent_thresold = 0.2,
                                  frac_TA_recurrent_threshold = 0.15) {
  sites <- stringr::str_c(somatic_candidates$chr, somatic_candidates$pos, sep = "-")
  # get recurrent sites
  recurrent_sites <- names(table(sites))[
    table(sites) > 1
  ]

  # label cells with too many recurrent sites
  single_strand_candidates$site <-
    stringr::str_c(single_strand_candidates$chr, single_strand_candidates$pos, sep = "-")
  cells_to_remove <- single_strand_candidates %>%
    dplyr::group_by(.data$cell_id) %>%
    dplyr::summarise(
      count_recurrent = sum(.data$site %in% recurrent_sites),
      frac_recurrent = mean(.data$site %in% recurrent_sites),
      frac_TA_recurrent = mean((.data$site %in% recurrent_sites) & stringr::str_detect(.data$mutsig, "T>A"))
    ) %>%
    dplyr::filter((.data$count_recurrent >= count_recurrent_threshold | .data$frac_recurrent > frac_recurrent_thresold),
                  .data$frac_TA_recurrent > frac_TA_recurrent_threshold) %>%
    dplyr::pull(.data$cell_id)
  single_strand_candidates$cell_qc_filter <- !(single_strand_candidates$cell_id %in% cells_to_remove)

  # remove recurrent single strand damage
  sites <- stringr::str_c(single_strand_candidates$chr, single_strand_candidates$pos, sep = "-")
  # get recurrent single strand damage
  recurrent_sites <- names(table(sites))[
    table(sites) > 1
  ]
  single_strand_candidates <- single_strand_candidates[!(single_strand_candidates$site %in% recurrent_sites), ]
  single_strand_candidates <- single_strand_candidates[, colnames(single_strand_candidates) != "site"]

  single_strand_candidates
}
