# 2xLPS fuzzy k-means soft clustering dataset

Mouse microglia cells from frontal cortex, striatum, and hippocampal
subregions. Mice were given two daily 0.5 mg/kg LPS intraperitoneal
injections or PBS vehicle injections and brains were collected 3 hours
after the final injection.

## Usage

``` r
data_2xLPS_mouse_fuzzykmeans
```

## Format

### `data_2xLPS_mouse_fuzzykmeans`

A data frame with 43,332 rows and 43 columns:

- Antibody:

  Cx3cr1, Iba1, P2ry12

- MouseID:

  1, 2, 3, 4, 5, 6

- Sex:

  F, M

- Treatment:

  PBS, 2xLPS

- BrainRegion:

  frontal cortex (FC), striatum (STR), or hippocampus (HC)

- Subregion:

  frontal cortex: infralimbic (IL), prelimbic (PL), anterior cingulate
  cortex (ACC); hippocampus: CA1, CA2, CA3, dentate gyrus (DG);
  stratium: caudate putamen (CP), nucleus accumbens (NA)

- ID:

  Individual Cell ID

- UniqueID:

  Unique descriptor for each cell in dataset

- 27 morphology features: Columns 9-35:

  Foreground pixels, Density of foreground pixels in hull area, Span
  ratio of hull: major/minor axis, Maximum span across hull, Area,
  Perimeter, Circularity, Width of bounding rectangle, Height of
  bounding rectangle, Maximum radius from hull's center of mass, Max/min
  radii from hull's center of mass, Relative variation in radii from
  hull's center of mass, Mean radius, Diameter of bounding circle,
  Maximum radius from circle's center of mass, Max/min radii from
  circle's center of mass, Relative variation in radii from circle's
  center of mass, Mean radius from circle's center of mass, \# of
  branches, \# of junctions, \# of end point voxels, \# of junction
  voxels, \# of slab voxels, Average branch length \# of triple points,
  \# of quadruple points, Maximum branch length

- PC 1-3 scores: Columns 36-38:

  PC1, PC2, PC3

- Fuzzy k-means cluster membership scores: Columns 39-42:

  Cluster 1, Cluster 2, Cluster 3, Cluster 4

- Hard clustering assignment: Cluster:

  1, 2, 3, 4
