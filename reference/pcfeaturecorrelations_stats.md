# Underlying stats for correlation heatmap exploring which morphologies describe the PCs

This function is complementary to the 'pcfeaturecorrelations' function.
'pcfeaturecorrelations_stats' outputs a dataframe containing the
correlation and pvalues related to your pcfeaturecorrelations heatmap.

## Usage

``` r
pcfeaturecorrelations_stats(
  pca_data,
  pc.start,
  pc.end,
  feature.start,
  feature.end,
  rthresh,
  pthresh
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

- data:

  is your input data frame
