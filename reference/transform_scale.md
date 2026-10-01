# Transform data: Scaling

'transform_scale' allows you to center and scale your dataframe within
each morphology measure to calculate z-scores using the `scale` function
in R.

## Usage

``` r
transform_scale(data, start, end)
```

## Arguments

- data:

  is your input dataframe

- start:

  is first column number of data you want to scale

- end:

  is last column number of data you want to scale
