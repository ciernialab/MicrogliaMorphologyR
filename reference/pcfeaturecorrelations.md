# What morphology features describe the PCs?

'pcfeaturecorrelations' allows you to perform correlations between the
top PCs and morphology features to assess which features are driving the
variability captured by each of the PCs.

## Usage

``` r
pcfeaturecorrelations(
  pca_data,
  pc.start,
  pc.end,
  feature.start,
  feature.end,
  rthresh,
  pthresh,
  title
)
```

## Arguments

- pc.start:

  is first PC you want included in analysis (e.g., PC1)

- pc.end:

  is last PC you want included in analysis (e.g., PC4)

- feature.start:

  is first column number of morphology measures

- feature.end:

  is last column number of morphology measures

- rthresh:

  is cutoff threshold for significant correlation values

- pthresh:

  is cutoff threshold for significant p-values

- title:

  is what you want to name your heatmap

- data:

  is your input data frame
