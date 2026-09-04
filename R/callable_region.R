#' Callable genome size as a function of single-cell depth
#'
#' Sums the SCAN2 depth profile over the region that passes the cell's own
#' static depth filters (single-cell depth >= \code{min.sc.dp} and bulk depth
#' >= \code{min.bulk.dp}), keeping the single-cell depth axis intact.
#'
#' @param results A SCAN2 object.
#' @param muttype \code{"snv"} or \code{"indel"}; the two use different depth
#'   minima.
#'
#' @return Named numeric vector: callable base pairs at each single-cell depth,
#'   from \code{min.sc.dp} upward. Aligns element-wise with
#'   \code{\link{prob_detect_ss_damage}(results, muttype)}.
#'
#' @export
callable_region_profile <- function(results, muttype) {
  dptab <- results@depth.profile$dptab
  max_col <- dim(dptab)[2]

  min_sc_dp   <- results@static.filter.params[[muttype]]$min.sc.dp
  min_bulk_dp <- results@static.filter.params[[muttype]]$min.bulk.dp

  rowSums(dptab[(min_sc_dp + 1):max_col, (min_bulk_dp + 1):max_col], na.rm = TRUE)
}

#' Effective callable region size
#'
#' The callable genome weighted by the probability of actually detecting a
#' single-stranded damage event there:
#' \deqn{\sum_{dp} \mathrm{callable}(dp) \times P(\mathrm{detect} | dp)}
#'
#' This is the denominator for single-stranded damage burden, and corrects for
#' the fact that low-depth regions are callable but insensitive.
#'
#' @param callable_profile Numeric vector from
#'   \code{\link{callable_region_profile}}.
#' @param detection_probability Numeric vector from
#'   \code{\link{prob_detect_ss_damage}}, the same length and depth alignment.
#'
#' @return A single number: effective callable base pairs.
#'
#' @export
effective_callable_region_size <- function(callable_profile, detection_probability) {
  if (length(callable_profile) != length(detection_probability)) {
    stop("callable_profile (", length(callable_profile), ") and ",
         "detection_probability (", length(detection_probability), ") must be ",
         "the same length; they are indexed by the same single-cell depths")
  }
  sum(callable_profile * detection_probability)
}
