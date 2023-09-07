## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(maldipickr)

## ----examples-process_spectra-------------------------------------------------
# Get an example directory of six Bruker MALDI Biotyper spectra
directory_biotyper_spectra <- system.file(
  "toy-species-spectra",
  package = "maldipickr"
)
# Import the six spectra
spectra_list <- import_biotyper_spectra(directory_biotyper_spectra)
# Transform the spectra signals according to Strejcek et al. (2018)
processed <- process_spectra(spectra_list)
# Overview of the list architecture that is returned
#  with the list of processed spectra, peaks identified and the
#  metadata table
str(processed, max.level = 2)
# A detailed view of the metadata with the median signal-to-noise
#  ratio (SNR) and the number of peaks
processed$metadata

## ----examples-merge_processed_spectra-----------------------------------------
# Get an example directory of six Bruker MALDI Biotyper spectra
directory_biotyper_spectra <- system.file(
  "toy-species-spectra",
  package = "maldipickr"
)
# Import the six spectra
spectra_list <- import_biotyper_spectra(directory_biotyper_spectra)
# Transform the spectra signals according to Strejcek et al. (2018)
processed <- process_spectra(spectra_list)
# Merge the spectra to produce the feature matrix
fm <- merge_processed_spectra(list(processed))
# The feature matrix has 6 spectra as rows and
#  35 peaks as columns
dim(fm)
# Notice the difference when the interpolation is turned off
fm_no_interpolation <- merge_processed_spectra(
  list(processed),
  interpolate_missing = FALSE
)
sum(fm == 0) # 0
sum(fm_no_interpolation == 0) # 68

# Multiple runs can be aggregated using list()
# Merge the spectra to produce the feature matrix
fm_all <- merge_processed_spectra(list(processed, processed, processed))
# The feature matrix has 3×6=18 spectra as rows and
#  35 peaks as columns
dim(fm_all)

## ----similarity, eval = FALSE-------------------------------------------------
#  # A. Compute the similarity matrix on the transposed feature matrix
#  #   using Pearson correlation coefficient
#  sim_matrix <- stats::cor(t(fm), method = "pearson")
#  
#  # B.1 Install the coop package
#  # install.packages("coop")
#  
#  # B.2 Compute the similarity matrix on the rows of the feature matrix
#  sim_matrix <- coop::tcosine(fm)

## ----examples-delineate_with_similarity---------------------------------------
# Toy similarity matrix between the six example spectra of
#  three species. The cosine metric is used and a value of
#  zero indicates dissimilar spectra and a value of one
#  indicates identical spectra.
cosine_similarity <- matrix(
  c(
    1, 0.79, 0.77, 0.99, 0.98, 0.98,
    0.79, 1, 0.98, 0.79, 0.8, 0.8,
    0.77, 0.98, 1, 0.77, 0.77, 0.77,
    0.99, 0.79, 0.77, 1, 1, 0.99,
    0.98, 0.8, 0.77, 1, 1, 1,
    0.98, 0.8, 0.77, 0.99, 1, 1
  ),
  nrow = 6,
  dimnames = list(
    c(
      "species1_G2", "species2_E11", "species2_E12",
      "species3_F7", "species3_F8", "species3_F9"
    ),
    c(
      "species1_G2", "species2_E11", "species2_E12",
      "species3_F7", "species3_F8", "species3_F9"
    )
  )
)
# Delineate clusters based on a 0.92 threshold applied
#  to the similarity matrix
delineate_with_similarity(cosine_similarity, threshold = 0.92)

## ----examples-set_reference_spectra-------------------------------------------
# Get an example directory of six Bruker MALDI Biotyper spectra
# Import the six spectra and
# Transform the spectra signals according to Strejcek et al. (2018)
processed <- system.file(
  "toy-species-spectra",
  package = "maldipickr"
) %>%
  import_biotyper_spectra() %>%
  process_spectra()

# Toy similarity matrix between the six example spectra of
#  three species. The cosine metric is used and a value of
#  zero indicates dissimilar spectra and a value of one
#  indicates identical spectra.
cosine_similarity <- matrix(
  c(
    1, 0.79, 0.77, 0.99, 0.98, 0.98,
    0.79, 1, 0.98, 0.79, 0.8, 0.8,
    0.77, 0.98, 1, 0.77, 0.77, 0.77,
    0.99, 0.79, 0.77, 1, 1, 0.99,
    0.98, 0.8, 0.77, 1, 1, 1,
    0.98, 0.8, 0.77, 0.99, 1, 1
  ),
  nrow = 6,
  dimnames = list(
    c(
      "species1_G2", "species2_E11", "species2_E12",
      "species3_F7", "species3_F8", "species3_F9"
    ),
    c(
      "species1_G2", "species2_E11", "species2_E12",
      "species3_F7", "species3_F8", "species3_F9"
    )
  )
)
# Delineate clusters based on a 0.92 threshold applied
#  to the similarity matrix
clusters <- delineate_with_similarity(
  cosine_similarity,
  threshold = 0.92
)

# Set reference spectra with the toy example
set_reference_spectra(clusters, processed$metadata)

## ----example-delineate_with_identification------------------------------------
report_unknown <- read_biotyper_report(
  system.file("biotyper_unknown.csv", package = "maldipickr")
)
delineate_with_identification(report_unknown)

## ----examples-import_spede_clusters-------------------------------------------
# Reformat the output from SPeDE table
# https://github.com/LM-UGent/SPeDE
import_spede_clusters(
  system.file("spede.csv", package = "maldipickr")
)

## ----examples-pick_spectra----------------------------------------------------
# 0. Load a toy example of a tibble of clusters created by
#   the `delineate_with_similarity` function.
clusters <- readRDS(
  system.file("clusters_tibble.RDS",
    package = "maldipickr"
  )
)
# 1. By default and if no other metadata are provided,
#   the function picks reference spectra for each clusters.
#
# N.B: The spectra `name` and `to_pick` columns are moved to the left
# only for clarity using the `relocate()` function.
#
pick_spectra(clusters) %>%
  dplyr::relocate(name, to_pick) # only for clarity

# 2.1 Simulate OD600 values with uniform distribution
#  for each of the colonies we measured with
#  the Bruker MALDI Biotyper
set.seed(104)
metadata <- dplyr::transmute(
  clusters,
  name = name, OD600 = runif(n = nrow(clusters))
)
metadata

# 2.2 Pick the spectra based on the highest
#   OD600 value per cluster
pick_spectra(clusters, metadata, "OD600") %>%
  dplyr::relocate(name, to_pick) # only for clarity

# 3.1 Say that the wells on the right side of the plate are
#   used for negative controls and should not be picked.
metadata <- metadata %>% dplyr::mutate(
  well = gsub(".*[A-Z]([0-9]{1,2}$)", "\\1", name) %>%
    strtoi(),
  is_edge = is_well_on_edge(
    well_number = well, plate_layout = 96, edges = "right"
  )
)

# 3.2 Pick the spectra after discarding (or soft masking)
#   the spectra indicated by the `is_edge` column.
pick_spectra(clusters, metadata, "OD600",
  soft_mask_column = "is_edge"
) %>%
  dplyr::relocate(name, to_pick) # only for clarity

# 4.1 Say that some spectra were picked before
#   (e.g., in the column F) in a previous experiment.
# We do not want to pick clusters with those spectra
#   included to limit redundancy.
metadata <- metadata %>% dplyr::mutate(
  picked_before = grepl("_F", name)
)
# 4.2 Pick the spectra from clusters without spectra
#   labeled as `picked_before` (hard masking).
pick_spectra(clusters, metadata, "OD600",
  hard_mask_column = "picked_before"
) %>%
  dplyr::relocate(name, to_pick) # only for clarity

