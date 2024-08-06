# MicrogliaMorphologyR 1.0

* First release of MicrogliaMorphologyR with BioRxiv preprint. November 2023.
* [Link to github release](https://github.com/ciernialab/MicrogliaMorphologyR/releases/tag/v1.0)

# MicrogliaMorphologyR 1.1

* Updated package with eNeuro publication. August 2024.
* [Link to github release]()

## Major changes
* Created github pages website and moved tutorial to dedicated section
* `celldensity()`: new function to calculate microglia density
* `featurecorrelations_stats()`: new function to obtain r and p-values underlying `featurecorrelations` output
* `pcfeaturecorrelations_stats()`: new function to obtain r and p-values underlying `pcfeaturecorrelations` output
* `stats_morphologymeasures.animal()`: added new argument 'type' to specify fixed or mixed effects modeling, removed levene test results (won't work when "+ variable" detected in model) and model summary (only summarizes model fit to last morphology measure in loop) in output

## Minor changes
* `transform_log()`: updated deprecated dplyr function funs()
* Updated tutorial to incoporate new changes and new functions
* Updated package documentation: authors, publication citation
