# Merge FracLac and AnalyzeSkeleton results together into one final cell-level data frame

'merge_data' merges your FracLac and AnalyzeSkeleton results by Name and
ID, gets rid of non-numerical data, and brushes up feature names. The
final output is a dataframe, where every row is a cell and every column
is an identifier or one of 27 unique morphology features.

## Usage

``` r
merge_data(fraclac, skeleton)
```

## Arguments

- fraclac:

  is your tidied fraclac output

- skeleton:

  is your tidied skeleton output
