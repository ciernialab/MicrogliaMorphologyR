MicrogliaMorphologyR
================

**Created**: 26 June, 2023 by Jenn Kim  
**Last updated**: 26 June, 2023

## Welcome to MicrogliaMorphologyR!

MicrogliaMorphologyR is an R package for microglia morphology analysis,
that is complimentary to ImageJ macro
[MicrogliaMorphology](https://github.com/ciernialab/MicrogliaMorphology).
Using MicrogliaMorphologyR, you can perform exploratory data analysis
and visualization of 27 different morphology features and perform
dimensionality reduction, clustering, and statistical analysis of your
data.

#### If you are using this tool, please cite the following publications:

-   Insert manuscript link

## Instructions on how to use MicrogliaMorphologyR

### install and load package

``` r
BiocManager::install('ciernialab/MicrogliaMorphologyR')
library(MicrogliaMorphologyR)
```

### load in your fraclac and skeleton data, tidy, and merge into final data frame

``` r
fraclac.dir <- "/Volumes/NINC/CierniaLab/Jenn/2022-2023/Paper1_MicrogliaMorphology/Collaborations/Jun_SzeleLab/Microglia_FracLac/20230313121607/"
skeleton.dir <- "/Volumes/NINC/CierniaLab/Jenn/2022-2023/Paper1_MicrogliaMorphology/Collaborations/Jun_SzeleLab/Microglia_SkeletonResults/"

fraclac <- fraclac_tidying(fraclac.dir)
skeleton <- skeleton_tidying(skeleton.dir)

data <- merge_data(fraclac, skeleton)

finaldata <- metadata_columns(data.test, c("Cohort","MouseID","Treatment","Sex","BrainRegion","SubRegion"),"_")
```

### exploratory data visualization

``` r
data2 <- finaldata %>% gather(measure, value, 4:ncol(finaldata))

# check for outliers
outliers_boxplots(data2)
outliers_distributions(data2)

# checking different normalization features
normalize_logplots(data2)
```
