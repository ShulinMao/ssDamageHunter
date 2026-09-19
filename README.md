# ssDamageHunter

Calls single-stranded DNA damage in single-cell WGS data.

## Layout

```
R/               package source
inst/scripts/    command-line entry point
inst/extdata/    synthetic example data
data-raw/        regenerates the example
tests/           unit tests
```

## Install

```r
remotes::install_github("ShulinMao/ssDamageHunter", dependencies = TRUE)
```

## Usage

```r
library(ssDamageHunter)

load("candidate.rda")
candidates <- call_ss_damage(results)
damage <- single_strand_damage_filter(candidates)
```

`recurrent_site_filter(candidates, damage)` adds the cross-cell filter once you
have more than one cell.

To estimate a burden, divide by the effective callable region size: the
callable genome weighted by the chance of detecting a single-stranded event at
each depth.

```r
callable <- callable_region_profile(results, "snv")
detect   <- prob_detect_ss_damage(results, "snv")

effective_callable_region_size(callable, detect)
```

## Method

Each candidate is scored under three hypotheses, each a binomial on the observed
alternate-read count at its expected allele frequency:

* double-stranded somatic mutation: the locus allele balance
* single-stranded damage: half that
* amplification or sequencing error: a quarter

A site is called single-stranded when that model fits best, the observed
frequency is near its expectation, and SCAN2 has not already called the site
somatic. Calls are then sorted within each cell by error likelihood and cut off
at a cumulative error budget, which bounds the expected false calls per cell
rather than per site. Sites seen in more than one cell are dropped, and cells
with an excess of recurrent calls are flagged.


## Demo

The package ships a small synthetic table to
check the install and see what the model does. It has 21 sites across three
cells, built so that every case the classifier has to separate is present: sites
at the single-stranded frequency, heterozygous germline-like sites at twice
that, low-frequency errors at half, and two sites SCAN2 has already called
somatic mutations.

```r
library(ssDamageHunter)

candidates <- read.csv(system.file("extdata", "example_candidates.csv",
                                   package = "ssDamageHunter"))
candidates <- cbind(candidates,
                    ss_damage_likelihoods(candidates$af, candidates$ab,
                                          candidates$dp, candidates$scalt))
single_strand_damage_filter(candidates)
```

Eight of the 21 sites survive, spread over all three cells and both mutation
types:

```
  cell_id chr       pos muttype    dp scalt    af
  cell_A  chr1  1050000 snv        80    20 0.25
  cell_A  chr1  1000000 snv        60    15 0.25
  cell_B  chr1  1250000 snv        64    21 0.328
  cell_B  chr1  1125000 snv        72    18 0.25
  cell_C  chr1  1200000 snv       100    25 0.25
  cell_A  chr1  1100000 indel      60    16 0.267
  cell_B  chr1  1175000 indel      56    14 0.25
  cell_C  chr1  1225000 indel      96    24 0.25
```

## License

MIT
