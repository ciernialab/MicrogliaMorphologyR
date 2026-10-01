# Dimensionality reduction using PCA

'pcadata' allows you to perform PCA analysis and update your dataframe
to include PCs of interest

## Usage

``` r
pcadata(data, featurestart, featureend, pc.start, pc.end)
```

## Arguments

- data:

  is your input data frame

- featurestart:

  is first column number of morphology measures

- featureend:

  is last column number of morphology measures

- pc.start:

  is first PC you want included (e.g., PC1)

- pc.end:

  is last PC you want included (e.g., PC4)
