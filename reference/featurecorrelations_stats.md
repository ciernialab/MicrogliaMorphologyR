# Underlying stats for correlation heatmap across morphology features

This function is complementary to the 'featurecorrelations' function.
'featurecorrelations_stats' outputs a dataframe containing the
correlation and pvalues related to your featurecorrelations heatmap.

## Usage

``` r
featurecorrelations_stats(data, featurestart, featureend, rthresh, pthresh)
```

## Arguments

- data:

  is your input data frame

- featurestart:

  is first column number of morphology measures

- featureend:

  is last column number of morphology measures

- rthresh:

  is cutoff threshold for significant correlation values

- pthresh:

  is cutoff threshold for significant p-values
