# What are your morphology cluster percentages across variables of interest?

'clusterpercentage' groups your data by variables of interest then
calculates percentage of each morphology cluster within those groups.

## Usage

``` r
clusterpercentage(data, clustercol, ...)
```

## Arguments

- data:

  is your input data frame

- clustercol:

  is the name of your column which contains cluster IDs. Make sure to
  put this in quotes.

- ...:

  list out your variables of interest without combining (e.g., Cohort,
  Sex, Treatment)
