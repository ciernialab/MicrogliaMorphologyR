# Exploratory data analysis: normalization methods

'normalize_logplots' allows you to check distributions of features after
log normalization (e.g., log10(value+0.1))

## Usage

``` r
normalize_logplots(data, x)
```

## Arguments

- data:

  is your final dataframe which contains gathered data (measure, value)
  format

- x:

  is the constant value you want to add (e.g., 0, 0.1, 1)
