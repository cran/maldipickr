## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(maldipickr)

## ----examples-read_biotyper_report--------------------------------------------
# Get a example Bruker report
biotyper <- system.file("biotyper.csv", package = "maldipickr")
# Import the report as a tibble
report_tibble <- read_biotyper_report(biotyper)
# Display the tibble
report_tibble

## ----examples-read_many_biotyper_reports--------------------------------------
# List of Bruker MALDI Biotyper reports
reports_paths <- system.file(
  c("biotyper.csv", "biotyper.csv", "biotyper.csv"),
  package = "maldipickr"
)
# Read the list of reports and combine them in a single tibble
read_many_biotyper_reports(
  reports_paths,
  report_ids = c("first", "second", "third"),
  # Additional metadata below are passed to dplyr::mutate
  growth_temperature = 37.0
)

## ----examples-import_biotyper_spectra-----------------------------------------
# Get an example directory of six Bruker MALDI Biotyper spectra
directory_biotyper_spectra <- system.file(
  "toy-species-spectra",
  package = "maldipickr"
)
# Import the six spectra
spectra_list <- import_biotyper_spectra(directory_biotyper_spectra)
# Display the list of spectra
spectra_list

## ----examples-check_spectra---------------------------------------------------
# Get an example directory of six Bruker MALDI Biotyper spectra
directory_biotyper_spectra <- system.file(
  "toy-species-spectra",
  package = "maldipickr"
)
# Import the six spectra
spectra_list <- import_biotyper_spectra(directory_biotyper_spectra)
# Display the list of checks, with FALSE where no anomaly is detected
check_spectra(spectra_list)
# The overall sanity can be checked with Reduce
Reduce(any, check_spectra(spectra_list)) # Should be FALSE

