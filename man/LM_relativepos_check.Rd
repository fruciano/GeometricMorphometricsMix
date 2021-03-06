% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LM_relativepos_check.R
\name{LM_relativepos_check}
\alias{LM_relativepos_check}
\title{Check the relative positions for a set of landmarks, compared to a reference specimen}
\usage{
LM_relativepos_check(
  Dataset,
  LM_to_check,
  reference_specimen = 1,
  only_by_landmark_order = FALSE,
  use_specimen_names = TRUE
)
}
\arguments{
\item{Dataset}{k x p x n array
of k landmarks and p dimensions for n observations (specimens).}

\item{LM_to_check}{Vector with the indices of which landmarks of Dataset should be tested}

\item{reference_specimen}{The index of the specimen to use as reference (by default, the first).}

\item{only_by_landmark_order}{If TRUE, each landmark in LM_to_check will be tested relative to the next in LM_to_check;
otherwise all the pairwise relative positions between landmarks will be run}

\item{use_specimen_names}{whether the names of specimens in Dataset should be used for the output}
}
\value{
The function outputs a list with the following elements:
 \describe{
  \item{potentially_problematic_idx}{Indices of potentially problematic specimens}
  \item{potentially_problematic_LMs}{List of potentially problematic landmarks and coordinates
  for the specimens in potentially_problematic_idx}
}
}
\description{
The function compares the relative position of a set of landmarks to the one observed in a reference specimen.
It also outputs which coordinates do not much the pattern observed in the reference specimen.
}
\details{
Compare relative positions of landmarks to the one observed in a reference specimen.
This function is useful to identify specimens with switched landmark positions and similar problems
when a dataset has been collected using consistent criteria
(e.g., a set of fish body pictures, with the mouth-tail axis approximately horizontal for all specimens).

For instance, if we want to check that landmarks 1, 2, and 3 are all aligned along the y coordinate in 2D,
we will use a specimen (usually the first), checking that it has been landmarked correctly.
Then, we will run the function on landmarks 1, 2, and 3 (by setting LM_to_check=c(1, 2, 3)).
The function will output a list of which specimens seem to be problematic and which landmarks and coordinates for these specimens seem problematic.
The putatively problematic coordinates will be indicated with "0".
Clearly, if we are only interested in 1, 2, and 3 being along y, it will not matter if in some case
we will find "0" under Coord1, but only if we find a "0" under Coord2.
}
\section{Notice}{

Generally, it is better to use this function with a small number
of carefully selected landmarks, instead of running it on all landmarks at the same time.

The parameter only_by_landmark_order allows to set whether all
combinations of landmarks should be tested (default), or not.

If only_by_landmark_order is set to TRUE, each landmark in LM_to_check
will be only tested with the next one. For instance if
LM_to_check=c(1,2,3) the function will only compare the first landmark with the second,
and the second with the third.

If the function cannot find potentially problematic specimens, a message
to this effect will be presented.
}

