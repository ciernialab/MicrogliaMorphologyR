# Tidy up metadata

'metadata_columns' allows you to extract individual metadata columns
from the UniqueID column. Each piece of metadata should be separated by
a common deliminator, 'sep' parameter

## Usage

``` r
metadata_columns(data, metadata, sep)
```

## Arguments

- data:

  is your tidied dataframe containing all of your numerical morphology
  data

- metadata:

  is your metadata listed out in order that they appear in 'Name'
  column: e.g., c("Cohort","Sex","MouseID")

- sep:

  is the character, "\_" or "-", that separates your metadata
