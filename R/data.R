#' 2xLPS mouse dataset from frontal cortex, striatum, and hippocampal subregions
#'
#' @format ## `data_2xLPS`
#' A data frame with 46,104 rows and 35 columns:
#' \describe{
#'   \item{Antibody}{Microglia marker used in immunofluorescent image}
#'   \item{MouseID}{Subject identifier}
#'   \item{Sex}{Sex}
#'   \item{Treatment}{PBS vehicle or 2xLPS i.p. injection at 0.5 mg/kg, spaced 24 hours apart}
#'   \item{BrainRegion}{frontal cortex, striatum, or hippocampus}
#'   \item{Subregion}{frontal cortex: infralimbic (IL), prelimbic (PL), anterior cingulate cortex (ACC); hippocampus: CA1, CA2, CA3, dentate gyrus (DG); stratium: caudate putamen (CP), nucleus accumbens (NA)}
#'   ...
#' }
"data_2xLPS"

#' EAE mouse dataset from frontal cortex, striatum, and hippocampal subregions
#'
#' @format ## `data_EAE`
#' A data frame with 46,104 rows and 35 columns:
#' \describe{
#'   \item{Antibody}{Microglia marker used in immunofluorescent image}
#'   \item{MouseID}{Subject identifier}
#'   \item{BrainRegion}{frontal cortex, striatum, or hippocampus}
#'   \item{Subregion}{frontal cortex: infralimbic (IL), prelimbic (PL), anterior cingulate cortex (ACC); hippocampus: CA1, CA2, CA3, dentate gyrus (DG); stratium: caudate putamen (CP), nucleus accumbens (NA)}
#'   ...
#' }
"data_EAE"

#' Fuzzy k-means dataset from 2xLPS
#'
#' @format ## `data_fuzzykmeans`
#' A data frame with 46,104 rows and 35 columns:
#' \describe{
#'   \item{Antibody}{Microglia marker used in immunofluorescent image}
#'   \item{MouseID}{Subject identifier}
#'   \item{BrainRegion}{frontal cortex, striatum, or hippocampus}
#'   \item{Subregion}{frontal cortex: infralimbic (IL), prelimbic (PL), anterior cingulate cortex (ACC); hippocampus: CA1, CA2, CA3, dentate gyrus (DG); stratium: caudate putamen (CP), nucleus accumbens (NA)}
#'   ...
#' }
"data_fuzzykmeans"

#' Fkm dataset from 2xLPS tided up
#'
#' @format ## `data_2xLPS_fuzzykmeans`
#' A data frame with 46,104 rows and 35 columns:
#' \describe{
#'   \item{Antibody}{Microglia marker used in immunofluorescent image}
#'   \item{MouseID}{Subject identifier}
#'   \item{BrainRegion}{frontal cortex, striatum, or hippocampus}
#'   \item{Subregion}{frontal cortex: infralimbic (IL), prelimbic (PL), anterior cingulate cortex (ACC); hippocampus: CA1, CA2, CA3, dentate gyrus (DG); stratium: caudate putamen (CP), nucleus accumbens (NA)}
#'   ...
#' }
"data_2xLPS_fuzzykmeans"
