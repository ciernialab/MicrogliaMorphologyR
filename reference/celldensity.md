# Calculate microglia density

'celldensity' loads in your Areas.csv file from MicrogliaMorphology,
which contains the areas of your images or rois within images, then
calculates the density of microglia cells for each image

## Usage

``` r
celldensity(AreasPath, SamplesizeDF)
```

## Arguments

- AreasPath:

  is the path to your Areas.csv file output from MicrogliaMorphology

- SampleSizeDF:

  is the dataframe of cell numbers that you generated from the
  'samplesize' function within MicrogliaMorphologyR
