#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# ssDamageHunter -- call single-stranded DNA damage from SCAN2 objects.
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(scan2)
  library(ssDamageHunter)
})

USAGE <- "
Usage: ssdamagehunter.R --out FILE (--scan2-object FILE... | --scan2-dir DIR)

Call single-stranded DNA damage from the somatic candidate set of one or more
SCAN2 objects.

Input (one of):
  --scan2-object FILE   A scan2_object.rda file. Repeat for several cells, or
                        pass a glob your shell expands.
  --scan2-dir DIR       Search DIR recursively for scan2_object.rda files.

Output:
  --out FILE            Where to write the calls. A .rds path writes the full
                        table; any other extension writes a tab-separated file.
  --candidates FILE     Also write the complete annotated candidate table,
                        before the likelihood filter.

Filter parameters (defaults match the manuscript):
  --max-af-ab-diff X       Max |af - single-stranded expectation|.  [0.1]
  --min-ss-likelihood X    Floor on the single-stranded likelihood. [0.05]
  --cum-error X            Per-cell cumulative error budget.        [0.005]
  --no-recurrent-filter    Skip the cross-cell recurrent-site filter. It needs
                           several cells to be meaningful, and is skipped
                           automatically for a single cell.

Other:
  -h, --help            Show this message.

Example:
  ssdamagehunter.R --scan2-dir /data/scan2/call_mutations --out damage.rds
"

# ------------------------------------------------------------------------------
# Argument parsing
# ------------------------------------------------------------------------------

parse_args <- function(argv) {
  opts <- list(
    scan2_object = character(0),
    scan2_dir = NULL,
    out = NULL,
    candidates = NULL,
    max_af_ab_diff = 0.1,
    min_ss_likelihood = 0.05,
    cum_error = 0.005,
    recurrent_filter = TRUE
  )

  # Flags that take a value, mapped to their name in `opts`.
  valued <- c("--scan2-object"      = "scan2_object",
              "--scan2-dir"         = "scan2_dir",
              "--out"               = "out",
              "--candidates"        = "candidates",
              "--max-af-ab-diff"    = "max_af_ab_diff",
              "--min-ss-likelihood" = "min_ss_likelihood",
              "--cum-error"         = "cum_error")
  numeric_opts <- c("max_af_ab_diff", "min_ss_likelihood", "cum_error")

  i <- 1
  while (i <= length(argv)) {
    arg <- argv[i]
    if (arg %in% c("-h", "--help")) {
      cat(USAGE)
      quit(status = 0)
    } else if (arg == "--no-recurrent-filter") {
      opts$recurrent_filter <- FALSE
      i <- i + 1
    } else if (arg %in% names(valued)) {
      if (i == length(argv)) stop(arg, " needs a value", call. = FALSE)
      name <- valued[[arg]]
      value <- argv[i + 1]
      if (name %in% numeric_opts) {
        value <- suppressWarnings(as.numeric(value))
        if (is.na(value)) stop(arg, " needs a number, got '", argv[i + 1], "'", call. = FALSE)
      }
      # --scan2-object accumulates; everything else is last-wins.
      opts[[name]] <- if (name == "scan2_object") c(opts[[name]], value) else value
      i <- i + 2
    } else {
      stop("unknown argument: ", arg, "\n", USAGE, call. = FALSE)
    }
  }
  opts
}

opts <- parse_args(commandArgs(trailingOnly = TRUE))

if (is.null(opts$out)) {
  cat(USAGE)
  stop("--out is required", call. = FALSE)
}

# ------------------------------------------------------------------------------
# Collect input objects
# ------------------------------------------------------------------------------

object_paths <- opts$scan2_object
if (!is.null(opts$scan2_dir)) {
  object_paths <- c(object_paths,
                    list.files(opts$scan2_dir, pattern = "^scan2_object\\.rda$",
                               recursive = TRUE, full.names = TRUE))
}
object_paths <- unique(object_paths)

if (length(object_paths) == 0) {
  stop("no SCAN2 objects found; pass --scan2-object or --scan2-dir", call. = FALSE)
}
missing <- object_paths[!file.exists(object_paths)]
if (length(missing) > 0) {
  stop("no such file: ", paste(missing, collapse = ", "), call. = FALSE)
}
message("Found ", length(object_paths), " SCAN2 object(s)")

# ------------------------------------------------------------------------------
# Call
# ------------------------------------------------------------------------------

candidates <- list()
for (path in object_paths) {
  load(path)  # defines `results`
  message("  ", results@single.cell)
  candidates[[results@single.cell]] <- call_ss_damage(results)
  rm(results)
  gc()
}
candidates <- do.call(rbind, candidates)
message(nrow(candidates), " quality-passing candidates from ",
        length(object_paths), " cell(s)")

damage <- single_strand_damage_filter(
  candidates,
  max_af_ab_diff = opts$max_af_ab_diff,
  min_single_strand_likelihood = opts$min_ss_likelihood,
  cum_error_likelihood = opts$cum_error)
message(nrow(damage), " single-stranded damage call(s)")

# The recurrent-site filter compares sites across cells, so it does nothing
# useful on a single cell.
if (opts$recurrent_filter && length(object_paths) > 1) {
  damage <- recurrent_site_filter(candidates, damage)
  message(nrow(damage), " call(s) after the recurrent-site filter")
} else if (opts$recurrent_filter) {
  message("Skipping the recurrent-site filter (needs more than one cell)")
}

# ------------------------------------------------------------------------------
# Write
# ------------------------------------------------------------------------------

write_table <- function(x, path) {
  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    saveRDS(x, path)
  } else {
    write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE)
  }
  message("Wrote ", path)
}

write_table(damage, opts$out)
if (!is.null(opts$candidates)) {
  write_table(candidates, opts$candidates)
}
