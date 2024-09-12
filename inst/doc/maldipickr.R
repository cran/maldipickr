## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(maldipickr)

## ----quickstart_report_data, eval = TRUE--------------------------------------
report_tbl <- read_biotyper_report(
  system.file("biotyper_unknown.csv", package = "maldipickr")
)
report_tbl %>%
  dplyr::select(name, bruker_species, bruker_log) %>% knitr::kable()

## ----quickstart_report_filter, eval = TRUE------------------------------------
report_tbl <- report_tbl %>%
  dplyr::mutate(
      bruker_species = dplyr::if_else(bruker_log >= 2, bruker_species,
                                      "not reliable identification")
  )
knitr::kable(report_tbl)

## ----quickstart_report_delineate, eval = TRUE---------------------------------
report_tbl %>%
  delineate_with_identification() %>%
  pick_spectra(report_tbl, criteria_column = "bruker_log") %>%
  dplyr::relocate(name, to_pick, bruker_species) %>% 
  knitr::kable()

## ----quickstart_spectra_data, eval = TRUE-------------------------------------
spectra_dir <- system.file("toy-species-spectra", package = "maldipickr")

processed <- spectra_dir %>%
  import_biotyper_spectra() %>%
  process_spectra()

## ----quickstart_spectra_delineate, eval = TRUE--------------------------------
processed %>%
  list() %>%
  merge_processed_spectra() %>%
  coop::tcosine() %>%
  delineate_with_similarity(threshold = 0.92) %>%
  set_reference_spectra(processed$metadata) %>%
  pick_spectra() %>%
  dplyr::relocate(name, to_pick) %>% 
  knitr::kable()

## ----session, eval = TRUE-----------------------------------------------------
sessionInfo()

