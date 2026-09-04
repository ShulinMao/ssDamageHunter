#' Standard SCAN2 quality tests applied to somatic candidates
#'
#' The candidate set is filtered by CIGAR-based hard/soft clipping tests, a low
#' mapping quality filter, a sequencing depth filter, a minimum number of
#' supporting alternate reads in the single cell, a maximum number of alternate
#' reads in the matched bulk, and a dbSNP overlap filter.
#'
#' @export
SCAN2_QUALITY_TESTS <- c(
  "cigar.hs.test",      # CIGAR-based hard clipping test
  "cigar.id.test",      # CIGAR-based indel test
  "lowmq.test",         # low mapping quality filter
  "dp.test",            # sequencing depth filter
  "min.sc.alt.test",    # minimum supporting alternate reads in the single cell
  "max.bulk.alt.test",  # maximum alternate reads in the matched bulk
  "dbsnp.test",         # dbSNP overlap filter
  "csf.test",
  "static.filter"
)

#' Quality tests applied to heterozygous germline training sites
#'
#' Germline sites are used as a negative control: they are true double-stranded
#' variants, so any of them classified as single-stranded is a misclassification.
#' Only the tests that make sense for a germline site are applied -- the bulk
#' and dbSNP filters would remove them by construction.
#'
#' @export
SCAN2_GERMLINE_QUALITY_TESTS <- c(
  "cigar.hs.test",
  "cigar.id.test",
  "dp.test",
  "min.sc.alt.test"
)

#' Call single-stranded DNA damage candidates from a SCAN2 object
#'
#' Extracts the somatic candidate set from a SCAN2 object, applies the standard
#' quality filters, and annotates every surviving site with the three competing
#' likelihoods (see \code{\link{ss_damage_likelihoods}}). The returned table is
#' the input to \code{\link{single_strand_damage_filter}}.
#'
#' @param results A SCAN2 object, as loaded from a \code{scan2_object.rda}
#'   produced by \code{scan2 call_mutations}.
#' @param cell_id Identifier to record in the \code{cell_id} column. Defaults to
#'   the object's own \code{single.cell} slot.
#' @param quality_tests Character vector of logical columns in
#'   \code{results@gatk} that must all be \code{TRUE}. Defaults to
#'   \code{\link{SCAN2_QUALITY_TESTS}}.
#' @param site_filter Name of the logical column selecting which sites to
#'   consider. \code{"somatic.candidate"} is the somatic candidate set;
#'   \code{"resampled.training.site"} selects heterozygous germline sites,
#'   used as a negative control (see
#'   \code{\link{SCAN2_GERMLINE_QUALITY_TESTS}}).
#'
#' @return A table of somatic candidates passing the quality filters, with the
#'   added columns \code{double_strand_likelihood},
#'   \code{single_strand_likelihood}, \code{error_likelihood},
#'   \code{af_ab_diff}, \code{sc_gp} and \code{cell_id}.
#'
#' @details
#' The single-cell genotype column in \code{results@gatk} is named after the
#' sample itself, so its position varies between SCAN2 versions (index 29 in
#' 82-column objects, 26 in 73-column ones). It is renamed to \code{sc_gp} by
#' matching \code{results@single.cell}, never by position.
#'
#' Sites where any quality test is \code{NA} are dropped. This matches the
#' behaviour of the original analysis scripts, which relied on
#' \code{data.table} treating \code{NA} in a logical subscript as \code{FALSE};
#' here it is made explicit with \code{which()} so the result does not depend on
#' whether \code{results@gatk} is a \code{data.table} or a \code{data.frame}.
#'
#' @export
call_ss_damage <- function(results,
                           cell_id = results@single.cell,
                           quality_tests = SCAN2_QUALITY_TESTS,
                           site_filter = "somatic.candidate") {
  gatk <- results@gatk

  missing_tests <- setdiff(c(quality_tests, site_filter), colnames(gatk))
  if (length(missing_tests) > 0) {
    stop("SCAN2 object is missing required column(s): ",
         paste(missing_tests, collapse = ", "))
  }

  somatic_mut <- gatk[which(gatk[[site_filter]]), ]

  # All quality tests must pass. NA counts as a failure.
  keep <- Reduce(`&`, lapply(quality_tests, function(test) somatic_mut[[test]]))
  somatic_mut <- somatic_mut[which(keep), ]

  likelihoods <- ss_damage_likelihoods(
    af    = somatic_mut$af,
    ab    = somatic_mut$ab,
    dp    = somatic_mut$dp,
    scalt = somatic_mut$scalt
  )
  for (column in colnames(likelihoods)) {
    somatic_mut[[column]] <- likelihoods[[column]]
  }

  somatic_mut <- rename_single_cell_genotype(somatic_mut, results@single.cell)
  somatic_mut$cell_id <- cell_id
  somatic_mut
}

#' Rename the sample-named single-cell genotype column to \code{sc_gp}
#'
#' @param x A table extracted from \code{results@gatk}.
#' @param single_cell The value of \code{results@single.cell}.
#'
#' @return \code{x} with the column named \code{single_cell} renamed to
#'   \code{sc_gp}.
#'
#' @keywords internal
#' @export
rename_single_cell_genotype <- function(x, single_cell) {
  position <- match(single_cell, colnames(x))
  if (is.na(position)) {
    stop("no column named '", single_cell, "' (results@single.cell) in the ",
         "SCAN2 table; cannot identify the single-cell genotype column")
  }
  colnames(x)[position] <- "sc_gp"
  x
}
