# Stats analysis: individual morphology measures

Linear mixed model to statistically assess how your experimental
variables of interest influence each morphology measure, at the subject
level. output\[1\] contains anova results, output\[2\] and output\[3\]
contain posthoc results, output\[4\] contains qqplots to check normality
assumptions, and output\[5\] contains shapiro test results.

## Usage

``` r
stats_morphologymeasures.animal(data, model, type, posthoc1, posthoc2, adjust)
```

## Arguments

- data:

  is your input data frame

- model:

  is your model (e.g., Value ~ Treatment\*Sex + (1\|MouseID))

- type:

  is the type of model you want to fit. "lm" for fixed effects or "lmer"
  for mixed effects (from the `lme4` package)

- posthoc1:

  is your posthoc comparisons (e.g., when considering sex:
  ~Treatment\|Sex)

- posthoc2:

  is your posthoc comparisons (e.g., when not considering sex:
  ~Treatment)

- adjust:

  is your method of multiple test correction (from `emmeans` package:
  "tukey","scheffe","sidak","dunnettx","mvt","holm", "hochberg",
  "hommel", "bonferroni", "BH", "BY", "fdr","none"). See "P-value
  adjustments" section under ?emmeans::summary.emmGrid for more
  information.

## Details

The stats_morphologymeasures.animal function fits a linear model using
the 'lm' function for fixed effects or 'lmer' function from the `lme4`
package for mixed effects for each morphology measure individually
within your dataset. Posthocs are run for each morphology measure
individually and bound together into the final dataframe that is output
by this function.
