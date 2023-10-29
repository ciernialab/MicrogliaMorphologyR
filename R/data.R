#' 2xLPS mouse dataset
#'
#' Mouse microglia cells from frontal cortex, striatum, and hippocampal subregions. Mice were given two daily 0.5 mg/kg LPS intraperitoneal injections or PBS vehicle injections and brains were collected 3 hours after the final injection.
#'
#' @format ## `data_2xLPS_mouse`
#' A data frame with 43,332 rows and 35 columns:
#' \describe{
#'   \item{Antibody}{Cx3cr1, Iba1, P2ry12}
#'   \item{MouseID}{1, 2, 3, 4, 5, 6}
#'   \item{Sex}{F, M}
#'   \item{Treatment}{PBS, 2xLPS}
#'   \item{BrainRegion}{frontal cortex (FC), striatum (STR), or hippocampus (HC)}
#'   \item{Subregion}{frontal cortex: infralimbic (IL), prelimbic (PL), anterior cingulate cortex (ACC); hippocampus: CA1, CA2, CA3, dentate gyrus (DG); stratium: caudate putamen (CP), nucleus accumbens (NA)}
#'   \item{ID}{Individual Cell ID}
#'   \item{UniqueID}{Unique descriptor for each cell in dataset}
#'   \item{27 morphology features: Columns 9-35}{Foreground pixels, Density of foreground pixels in hull area, Span ratio of hull: major/minor axis, Maximum span across hull, Area, Perimeter, Circularity, Width of bounding rectangle, Height of bounding rectangle, Maximum radius from hull's center of mass, Max/min radii from hull's center of mass, Relative variation in radii from hull's center of mass, Mean radius, Diameter of bounding circle, Maximum radius from circle's center of mass, Max/min radii from circle's center of mass, Relative variation in radii from circle's center of mass, Mean radius from circle's center of mass, # of branches, # of junctions, # of end point voxels, # of junction voxels, # of slab voxels, Average branch length # of triple points, # of quadruple points, Maximum branch length}
#'   ...
#' }
"data_2xLPS_mouse"

#' 2xLPS fuzzy k-means soft clustering dataset
#'
#' Mouse microglia cells from frontal cortex, striatum, and hippocampal subregions. Mice were given two daily 0.5 mg/kg LPS intraperitoneal injections or PBS vehicle injections and brains were collected 3 hours after the final injection.
#'
#' @format ## `data_2xLPS_mouse_fuzzykmeans`
#' A data frame with 43,332 rows and 43 columns:
#' \describe{
#'   \item{Antibody}{Cx3cr1, Iba1, P2ry12}
#'   \item{MouseID}{1, 2, 3, 4, 5, 6}
#'   \item{Sex}{F, M}
#'   \item{Treatment}{PBS, 2xLPS}
#'   \item{BrainRegion}{frontal cortex (FC), striatum (STR), or hippocampus (HC)}
#'   \item{Subregion}{frontal cortex: infralimbic (IL), prelimbic (PL), anterior cingulate cortex (ACC); hippocampus: CA1, CA2, CA3, dentate gyrus (DG); stratium: caudate putamen (CP), nucleus accumbens (NA)}
#'   \item{ID}{Individual Cell ID}
#'   \item{UniqueID}{Unique descriptor for each cell in dataset}
#'   \item{27 morphology features: Columns 9-35}{Foreground pixels, Density of foreground pixels in hull area, Span ratio of hull: major/minor axis, Maximum span across hull, Area, Perimeter, Circularity, Width of bounding rectangle, Height of bounding rectangle, Maximum radius from hull's center of mass, Max/min radii from hull's center of mass, Relative variation in radii from hull's center of mass, Mean radius, Diameter of bounding circle, Maximum radius from circle's center of mass, Max/min radii from circle's center of mass, Relative variation in radii from circle's center of mass, Mean radius from circle's center of mass, # of branches, # of junctions, # of end point voxels, # of junction voxels, # of slab voxels, Average branch length # of triple points, # of quadruple points, Maximum branch length}
#'   \item{PC 1-3 loadings: Columns 36-38}{PC1, PC2, PC3}
#'   \item{Fuzzy k-means cluster membership scores: Columns 39-42}{Cluster 1, Cluster 2, Cluster 3, Cluster 4}
#'   \item{Hard clustering assignment: Cluster}{1, 2, 3, 4}
#'   ...
#' }
"data_2xLPS_mouse_fuzzykmeans"

#' 1xLPS mouse dataset
#'
#' Mouse microglia cells from regions of interest encompassing hippocampus, white matter tracks, and overlaying cortex. Mice were given a single 1 mg/kg LPS intraperitoneal injection or PBS vehicle and brains were collected 24 hours later.
#'
#' @format ## `data_1xLPS_mouse`
#' A data frame with 15,293 rows and 34 columns:
#' \describe{
#'   \item{Antibody}{Iba1}
#'   \item{MouseID}{1_F, 1_M, 2_F, 2_M, 3_F, 4_F, 4_M, 5_M}
#'   \item{Sex}{F, M}
#'   \item{Treatment}{PBS, LPS}
#'   \item{ID}{Individual Cell ID}
#'   \item{UniqueID}{Unique descriptor for each cell in dataset}
#'   \item{27 morphology features: Columns 7-33}{Foreground pixels, Density of foreground pixels in hull area, Span ratio of hull: major/minor axis, Maximum span across hull, Area, Perimeter, Circularity, Width of bounding rectangle, Height of bounding rectangle, Maximum radius from hull's center of mass, Max/min radii from hull's center of mass, Relative variation in radii from hull's center of mass, Mean radius, Diameter of bounding circle, Maximum radius from circle's center of mass, Max/min radii from circle's center of mass, Relative variation in radii from circle's center of mass, Mean radius from circle's center of mass, # of branches, # of junctions, # of end point voxels, # of junction voxels, # of slab voxels, Average branch length # of triple points, # of quadruple points, Maximum branch length}
#'   ...
#' }
"data_1xLPS_mouse"

#' 2xLPS mouse dataset from frontal cortex, striatum, and hippocampal subregions
#'
#' @format ## `data_ImageTypeComparison`
#' A data frame with 60 rows and 13 columns:
#' \describe{
#'   \item{CellID}{ameboid1-5, hypertrophic1-5, ramified1-5, rod1-5}
#'   \item{MorphologyClass}{ameboid, hypertrophic, ramified, rod}
#'   \item{ID}{1, 2, 3, 4, 5}
#'   \item{ImageType}{2D, 3D, EDF}
#'   \item{AnalyzeSkeleton measures: Columns 5-13}{# of branches, # of junctions, # of end point voxels, # of junction voxels, # of slab voxels, Average branch length # of triple points, # of quadruple points, Maximum branch length}
#'   ...
#' }
"data_ImageTypeComparison"

