% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pick_spectra.R
\name{pick_spectra}
\alias{pick_spectra}
\title{Cherry-pick Bruker MALDI Biotyper spectra}
\usage{
pick_spectra(
  cluster_df,
  metadata_df = NULL,
  criteria_column = NULL,
  hard_mask_column = NULL,
  soft_mask_column = NULL,
  is_descending_order = TRUE,
  is_sorted = FALSE
)
}
\arguments{
\item{cluster_df}{A tibble with clusters information
from the \link{delineate_with_similarity} or the \link{import_spede_clusters} function.}

\item{metadata_df}{Optional tibble with relevant metadata to guide the
picking process (e.g., OD600).}

\item{criteria_column}{Optional character indicating the column in \code{metadata_df}
to be used as a criteria.}

\item{hard_mask_column}{Column name in the \code{cluster_df} or \code{metadata_df} tibble indicating whether the spectra, \strong{and the clusters to which they belong} should be discarded (\code{TRUE}) or not (\code{FALSE}) before the picking decision.}

\item{soft_mask_column}{Column name in the \code{cluster_df} or \code{metadata_df} tibble indicating whether the spectra should be discarded (\code{TRUE}) or not (\code{FALSE}) before the picking decision.}

\item{is_descending_order}{Optional logical indicating whether to sort the \code{criteria_column} from the highest-to-lowest value (\code{TRUE}) or lowest-to-highest (\code{FALSE}).}

\item{is_sorted}{Optional logical to indicate that the \code{cluster_df} is
already sorted by cluster based on (usually multiple) internal criteria to
pick the first of each cluster. This flag is \strong{overridden} if a \code{metadata_df}
is provided.}
}
\value{
A tibble with as many rows as \code{cluster_df} with an additional logical
column named \code{to_pick} to indicate whether the colony associated to the spectra
should be picked. If \code{metadata_df} is provided, then additional columns from
this tibble are added to the returned tibble.
}
\description{
Using the clusters information, and potential additional metadata as external
criteria, spectra are labeled as to be picked for each cluster. Note that some
spectra and therefore clusters can be explicitly removed (\emph{masked})
from the picking decision if they have been previously picked
or should be discarded, using logical columns in the metadata table.
If no metadata are provided, the reference spectra of each cluster will be picked.
}
\examples{
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
pick_spectra(clusters) \%>\%
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
pick_spectra(clusters, metadata, "OD600") \%>\%
  dplyr::relocate(name, to_pick) # only for clarity

# 3.1 Say that the wells on the right side of the plate are
#   used for negative controls and should not be picked.
metadata <- metadata \%>\% dplyr::mutate(
  well = gsub(".*[A-Z]([0-9]{1,2}$)", "\\\\1", name) \%>\%
    strtoi(),
  is_edge = is_well_on_edge(
    well_number = well, plate_layout = 96, edges = "right"
  )
)

# 3.2 Pick the spectra after discarding (or soft masking)
#   the spectra indicated by the `is_edge` column.
pick_spectra(clusters, metadata, "OD600",
  soft_mask_column = "is_edge"
) \%>\%
  dplyr::relocate(name, to_pick) # only for clarity

# 4.1 Say that some spectra were picked before
#   (e.g., in the column F) in a previous experiment.
# We do not want to pick clusters with those spectra
#   included to limit redundancy.
metadata <- metadata \%>\% dplyr::mutate(
  picked_before = grepl("_F", name)
)
# 4.2 Pick the spectra from clusters without spectra
#   labeled as `picked_before` (hard masking).
pick_spectra(clusters, metadata, "OD600",
  hard_mask_column = "picked_before"
) \%>\%
  dplyr::relocate(name, to_pick) # only for clarity
}
\seealso{
\link{delineate_with_similarity}, \link{set_reference_spectra}. For a useful utility function to soft-mask specific spectra: \link{is_well_on_edge}.
}
