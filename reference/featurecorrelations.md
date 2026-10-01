# Correlation heatmap across morphology features

'featurecorrelations' allows you to generate a heatmap depicting
significant correlations across features

## Usage

``` r
featurecorrelations(data, featurestart, featureend, rthresh, pthresh, title)
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

- title:

  is what you want to name your heatmap
