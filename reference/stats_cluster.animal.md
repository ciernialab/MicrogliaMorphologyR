# Stats analysis: cluster-level changes

Linear mixed model to statistically assess how your experimental
variables of interest influence cluster percentages, at animal level.
output\[1\] contains anova results, output\[2\] and output\[3\] contain
posthoc results, output\[4\] contains model fit checks (read more under
Details section), and output\[5\] contains information about the model.

## Usage

``` r
stats_cluster.animal(data, model, posthoc1, posthoc2, adjust)
```

## Arguments

- data:

  is your input data frame

- model:

  is your linear mixed model (e.g., Value ~ Cluster\*Treatment +
  (1\|MouseID))

- posthoc1:

  is your first set of posthoc comparisons (e.g., ~Cluster\|Treatment)

- posthoc2:

  is your second set of posthoc comparisons (e.g., ~Cluster)

- adjust:

  is your method of multiple test correction (from `emmeans` package:
  "tukey","scheffe","sidak","dunnettx","mvt","holm", "hochberg",
  "hommel", "bonferroni", "BH", "BY", "fdr","none"). See "P-value
  adjustments" section under ?emmeans::summary.emmGrid for more
  information.

## Details

The stats_cluster.animal function fits a generalized linear mixed model
on your dataset to a beta distribution, which is suitable for values
like percentages or probabilities that are constrained to a range of
0-1, using the `glmmTMB` package. Part of the output includes a check of
the model fit using the `DHARMa` package, which "uses a simulation-based
approach to create readily interpretable scaled (quantile) residuals for
fitted (generalized) linear mixed models." The function creates two
`DHARMa` plots, contained in output\[4\]. You can read more about how to
interpret model fit using `DHARMa` by reading the package
[vignette](https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html).
