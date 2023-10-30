#' Tidy up your FracLac output
#'
#' 'fraclac_tidying' loads in your FracLac output and cleans up names of cells
#'
#' @param dir is your directory containing the 'Hull and Circles Results.txt' file from FracLac output
#' @export
fraclac_tidying <- function(dir){
  setwd(dir)
  FracLac1 <- read.csv("Hull and Circle Results.txt", sep='\t')
  y <- FracLac1 %>%
    separate(1, into=c("Name","trash"), sep="tif_thresholdedtif_") %>%
    separate(trash, into=c("ID","trash1"), sep="tif") %>%
    select(-trash1)
  y
}

#' Concatenate your individual AnalyzeSkeleton output files into one dataframe
#'
#' 'skeleton_tidying' loads in all of your AnalyzeSkeleton output and clean up IDs to match FracLac data for downstream merging of data at cell level.
#'
#' @param dir is your directory containing all of your individual AnalyzeSkeleton .csv result files output from MicrogliaMorphology
#' @export
skeleton_tidying <- function(dir){
  setwd(dir)
  files <- list.files()
  nf = length(files)
  good = vector("list", nf)
  for(i in 1:nf){
    tmp = try(read.csv(files[i]))
    if(inherits(tmp,"try-error")) next
    tmp$Name=files[i]
    if(nrow(tmp)==1)
      good[[i]]=tmp
  }
  finaldf2 <- do.call(rbind,good)

  y <- finaldf2 %>%
    separate(Name, into=c("Name","trash"), sep=".tif_thresholded.tif_") %>%
    separate(trash, into=c("ID","junk2"), sep=".tif_results") %>% select(-junk2) %>%
    unite("UniqueID",c("Name","ID"),sep="_",remove=FALSE)
  y
}

#' Merge FracLac and AnalyzeSkeleton results together into one final cell-level data frame
#'
#' 'merge_data' merges your FracLac and AnalyzeSkeleton results by Name and ID, gets rid of non-numerical data, and brushes up feature names.
#' The final output is a dataframe, where every row is a cell and every column is an identifier or one of 27 unique morphology features.
#'
#' @param fraclac is your tidied fraclac output
#' @param skeleton is your tidied skeleton output
#' @export
merge_data <- function(fraclac, skeleton){
  data <- inner_join(fraclac, skeleton, by=c("Name","ID"))

  data <- data %>% select(-c(TOTAL.PIXELS,
                             Hull.s.Centre.of.Mass,
                             Method.Used.to.Calculate.Circle,
                             Hull.s.Centre.of.Mass,
                             Circle.s.Centre, X))
  data <- data[, c(1:2,30,3:29)]

  # brush up feature names
  feature_names <- c("Foreground pixels",
                     "Density of foreground pixels in hull area",
                     "Span ratio of hull (major/minor axis)",
                     "Maximum span across hull",
                     "Area",
                     "Perimeter",
                     "Circularity",
                     "Width of bounding rectangle",
                     "Height of bounding rectangle",
                     "Maximum radius from hull's center of mass",
                     "Max/min radii from hull's center of mass",
                     "Relative variation (CV) in radii from hull's center of mass",
                     "Mean radius",
                     "Diameter of bounding circle",
                     "Maximum radius from circle's center of mass",
                     "Max/min radii from circle's center of mass",
                     "Relative variation (CV) in radii from circle's center of mass",
                     "Mean radius from circle's center of mass",
                     "# of branches",
                     "# of junctions",
                     "# of end point voxels",
                     "# of junction voxels",
                     "# of slab voxels",
                     "Average branch length",
                     "# of triple points",
                     "# of quadruple points",
                     "Maximum branch length")

  colnames(data)[4:30] <- feature_names
  data
}


#' Tidy up metadata
#'
#' 'metadata_columns' allows you to extract individual metadata columns from the UniqueID column.
#' Each piece of metadata should be separated by a common deliminator, 'sep' parameter
#'
#' @param data is your tidied dataframe containing all of your numerical morphology data
#' @param metadata is your metadata listed out in order that they appear in 'Name' column: e.g., c("Cohort","Sex","MouseID")
#' @param sep is the character, "_" or "-", that separates your metadata
#' @export
metadata_columns <- function(data, metadata, sep){
  final <- data %>% separate(Name, into=metadata, sep=sep)
  final
}

#' Exploratory data analysis: outlier detection
#'
#' 'outliers_boxplots' generates boxplots of morphology features to visualize those with broader ranges which might dominate the analysis downstream.
#'
#' @param data is your final dataframe which contains gathered data (measure, value) format
#' @export
outliers_boxplots <- function(data){
  data %>%
    ggplot(aes(x=measure, y=value)) +
    geom_boxplot() +
    geom_jitter(width=0.15)+
    theme(axis.text.x = element_text(angle = 90)) +
    ggtitle("Range of values for each morphology measure")
}

#' Exploratory data analysis: outlier detection
#'
#' 'outliers_distributions' generates plots of feature distributions to determine which features have skewed distributions and how to normalize your values for analysis downstream.
#'
#' @param data is your final dataframe which contains gathered data (measure, value) format
#' @export
outliers_distributions <- function(data){
  data %>%
    ggplot(aes(x = value)) +
    geom_histogram(bins = 40) +
    facet_wrap(~measure, scales = "free_x") +
    ggtitle("Distributions of morphology measures")
}

#' Exploratory data analysis: normalization methods
#'
#' 'normalize_logplots' allows you to check distributions of features after log normalization (e.g., log10(value+0.1))
#'
#' @param data is your final dataframe which contains gathered data (measure, value) format
#' @param x is the constant value you want to add (e.g., 0, 0.1, 1)
#' @export
normalize_logplots <- function(data,x){
  data %>%
    ggplot(aes(x = log(value+x))) +
    geom_histogram(bins = 40) +
    facet_wrap(~measure, scales = "free") +
    #theme(strip.text.x = element_text(size=8)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_y_continuous(name = "Count") +
    ggtitle("Distribution of measures after log transforming")
}

#' Exploratory data analysis: scaled plots
#'
#' 'normalize_scaled' allows you to check distributions of features after scaling (e.g., scale(value))
#'
#' @param data is your final dataframe which contains gathered data (measure, value) format
#' @export
normalize_scaled <- function(data){
  data %>%
    ggplot(aes(x = scale(value))) +
    geom_histogram(bins = 40) +
    facet_wrap(~measure, scales = "free") +
    #theme(strip.text.x = element_text(size=8)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_continuous(name = "Morphology measure (scale(value))") +
    scale_y_continuous(name = "Count") +
    ggtitle("Distributions of measures after scaling")
}

#' Exploratory data analysis: min-max scaled plots
#'
#' 'normalize_minmax' allows you to check distributions of features after min-max normalizing to get your values for each measure constrained to 0-1 range.
#'
#' @param data is your final dataframe which contains gathered data (measure, value) format
#' @export
normalize_minmax <- function(data){
  data %>%
    ggplot(aes(x = minmax(value))) +
    geom_histogram(bins = 40) +
    facet_wrap(~measure, scales = "free") +
    #theme(strip.text.x = element_text(size=8)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_continuous(name = "Morphology measure (minmax(value))") +
    scale_y_continuous(name = "Count") +
    ggtitle("Distribution of measures after min-max normalizing")
}

#' Log scale
#'
#' 'logscale' allows you to log scale values
#'
#' @param x is your input values
#' @param y is the constant value you want to add (e.g., 0, 0.1, 1)
#' @export
logscale <- function(x,y){
  data <- log(x+y)
  data}

#' Transform data: Log scale
#'
#' 'transform_log' allows you to log scale your dataframe
#'
#' @param data is your input dataframe
#' @param x is the constant value you want to add (e.g., 0.1, 1)
#' @param start is first column number of morphology measures
#' @param end is last column number of morphology measures
#' @export
transform_log <- function(data,x,start,end){
  data %>% mutate_at(start:end, funs(logscale(.,x)))
}

#' Min-max scaling
#'
#' 'minmax' allows you to min-max normalize values
#'
#' @param x is your input values
#' @export
minmax <- function(x){
  data <- (x- min(x)) /(max(x)-min(x))
  data}

#' Transform data: Min-max scaling
#'
#' 'transform_minmax' allows you to min-max scale your dataframe within each morphology measure
#'
#' @param data is your input dataframe
#' @param start is first column number of morphology measures
#' @param end is last column number of morphology measures
#' @export
transform_minmax <- function(data,start,end){
  data %>% mutate_at(start:end, minmax)
}

#' Transform data: Scaling
#'
#' 'transform_scale' allows you to center and scale your dataframe within each morphology measure
#' to calculate z-scores using the `scale` function in R.
#'
#' @param data is your input dataframe
#' @param start is first column number of data you want to scale
#' @param end is last column number of data you want to scale
#' @export
transform_scale <- function(data,start,end){
  data %>% mutate_at(start:end, scale)
}

#' Sample size
#'
#' 'samplesize' allows you to obtain sample size of your experimental groups based on variables of interest
#'
#' @param data is your input dataframe
#' @param ... list out your variables of interest without combining (e.g., Cohort, Sex, Treatment)
#' @export
samplesize <- function(data,...){
  data %>% group_by(...) %>% summarise(num=n())
}

#' Correlation heatmap across morphology features
#'
#' 'featurecorrelations' allows you to generate a heatmap depicting significant correlations across features
#'
#' @param data is your input data frame
#' @param featurestart is first column number of morphology measures
#' @param featureend is last column number of morphology measures
#' @param rthresh is cutoff threshold for significant correlation values
#' @param pthresh is cutoff threshold for significant p-values
#' @param title is what you want to name your heatmap
#' @export
featurecorrelations <- function(data,featurestart,featureend,rthresh,pthresh,title){
  # correlations with p-values
  hi <- rcorr(as.matrix(data[,featurestart:featureend]), type="spearman")
  mat1 <- hi$r
  filter <- which(abs(mat1)<rthresh)

  mat2 <- hi$P
  mat2 <- round(mat2,3)

  mat2[mat2 < pthresh] <- "*" # significant p-values
  mat2[mat2 > pthresh] <- "" # insignificant p-values
  mat2[is.na(mat2)] <- "" # NAs (should be 27)
  mat2[filter] <- "" # correlation values that are less than 0.5 (weak correlations)

  # make heatmap
  pheatmap(hi$r, display_numbers = mat2, border_color=NA,
           fontsize_number=12,
           #fontsize=10, fontsize_row=10, fontsize_col=10,
           main=title)

  #ComplexHeatmap::pheatmap(hi$r, display_numbers = mat2, fontsize_number=12,
  #                         fontsize=10, fontsize_row=10, fontsize_col=10,
  #                         legend = TRUE, heatmap_legend_param = list(title = "testing", annotation_name_side = "top"),
  #                         main=title)
}

#' Dimensionality reduction using PCA and elbow/scree plot to get variance explained
#'
#' 'pcadata_elbow' allows you to perform PCA analysis and update your dataframe
#'
#' @param data is your input data frame
#' @param featurestart is first column number of morphology measures
#' @param featureend is last column number of morphology measures
#' @export
pcadata_elbow <- function(data, featurestart, featureend){
  prin_comp <- prcomp(data[,featurestart:featureend], scale=TRUE)
  fviz_eig(prin_comp, addlabels=TRUE, main="Variance explained by PCs")
}

#' Dimensionality reduction using PCA
#'
#' 'pcadata' allows you to perform PCA analysis and update your dataframe to include PCs of interest
#'
#' @param data is your input data frame
#' @param featurestart is first column number of morphology measures
#' @param featureend is last column number of morphology measures
#' @param pc.start is first PC you want included (e.g., PC1)
#' @param pc.end is last PC you want included (e.g., PC4)
#' @export
pcadata <- function(data, featurestart, featureend, pc.start, pc.end){
  prin_comp <- prcomp(data[,featurestart:featureend], scale=TRUE)
  components <- prin_comp[["x"]]
  components <- data.frame(components)
  pca_data <- cbind(prin_comp$x[,pc.start:pc.end], data)
  pca_data
}


#' What morphology features describe the PCs?
#'
#' 'pcfeaturecorrelations' allows you to perform correlations between the top PCs and morphology features to assess which features are driving the variability captured by each of the PCs.
#'
#' @param data is your input data frame
#' @param pc.start is first PC you want included in analysis (e.g., PC1)
#' @param pc.end is last PC you want included in analysis (e.g., PC4)
#' @param feature.start is first column number of morphology measures
#' @param feature.end is last column number of morphology measures
#' @param rthresh is cutoff threshold for significant correlation values
#' @param pthresh is cutoff threshold for significant p-values
#' @param title is what you want to name your heatmap
#' @export
pcfeaturecorrelations <- function(pca_data, pc.start, pc.end, feature.start, feature.end, rthresh, pthresh, title){
  pc.start2 <- colnames(pca_data[pc.start])
  pc.end2 <- colnames(pca_data[pc.end])
  morphologyfeatures <- colnames(pca_data[,feature.start:feature.end])
  PCs <- colnames(pca_data[,pc.start:pc.end])
  master_correlations <- NULL
  master_pvalues <- NULL
  for(p in PCs){
    all_correlations <- NULL
    all_pvalues <- NULL
    for(m in morphologyfeatures){
      # correlations
      master <- rcorr(as.matrix(pca_data[,c(p, m)]), type="spearman")
      correlations <- as.data.frame(master$r)
      correlations <- rownames_to_column(correlations, "PC")
      correlations <- correlations %>% dplyr::select("PC", m) %>% .[1,]
      names(correlations)[names(correlations) == m] <- "value"
      correlations$measure <- paste(m)
      all_correlations <- rbind(all_correlations, correlations)

      # p-values
      pvalues <- as.data.frame(master$P)
      pvalues <- rownames_to_column(pvalues, "PC")
      pvalues <- pvalues %>% dplyr::select("PC", m) %>% .[1,]
      names(pvalues)[names(pvalues) == m] <- "value"
      pvalues$measure <- paste(m)
      all_pvalues <- rbind(all_pvalues, pvalues)
    }
    master_correlations <- rbind(master_correlations, all_correlations)
    master_pvalues <- rbind(master_pvalues, all_pvalues)
  }

  # heatmap of correlations PC vs. measures

  # correlations matrix
  master_correlations_PC = list()
  for (p in PCs){
    master_correlations_PC[[p]] =
      master_correlations %>% filter(PC==p) %>% dplyr::rename(!!p :=value) %>% select(-PC)
  }

  testing <- master_correlations_PC %>% reduce(inner_join, by="measure", keep=FALSE)

  heatmap_PC_correlations <- testing[,c(2,1,3:ncol(testing))]
  heatmap_PC_correlations <- heatmap_PC_correlations %>% column_to_rownames(var="measure")

  rownames(heatmap_PC_correlations) <- morphologyfeatures

  # pvalues matrix
  master_pvalues_PC = list()
  for (p in PCs){
    master_pvalues_PC[[p]] =
      master_pvalues %>% filter(PC==p) %>% dplyr::rename(!!p :=value) %>% select(-PC)
  }

  testing <- master_pvalues_PC %>% reduce(inner_join, by="measure", keep=FALSE)

  heatmap_PC_pvalues <- testing[,c(2,1,3:ncol(testing))]
  heatmap_PC_pvalues <- heatmap_PC_pvalues %>% column_to_rownames(var="measure")

  rownames(heatmap_PC_pvalues) <- morphologyfeatures

  # re-formatting for heatmap to only show significant values
  bat1 <- heatmap_PC_correlations %>% select(pc.start2:pc.end2) %>% as.matrix()
  filter <- which(abs(bat1)<rthresh)

  bat2 <- heatmap_PC_pvalues %>% select(pc.start2:pc.end2) %>% as.matrix()
  bat2 <- round(bat2,3)

  bat2[bat2 < pthresh] <- "*" # significant p-values
  bat2[bat2 > pthresh] <- "" # insignificant p-values
  bat2[is.na(bat2)] <- "" # NAs (should be 27)
  bat2[filter] <- ""

  pheatmap(bat1, display_numbers = bat2, border_color=NA,
           fontsize_number=12,
           #fontsize=14, fontsize_row=12, fontsize_col=12,
           angle_col=0, main=title)
}

#' Explore how experimental variables describe data
#'
#' 'plots_expvariable' colors points in pc space by experimental variables
#'
#' @param data is your input data frame
#' @param pc.xaxis is the pc values you want on x-axis (e.g, PC1)
#' @param pc.yaxis is the pc values you want on y-axis (e.g., PC2)
#' @export
plots_expvariable <- function(data, pc.xaxis, pc.yaxis){
  pc.xaxis <- sym(pc.xaxis)
  pc.yaxis <- sym(pc.yaxis)

  variables <- unique(as.character(data$variable))
  plots_variables = list()
  for(v in variables){
    plots_variables[[v]] =
      data %>%
      filter(variable==v) %>%
      ggplot(aes(x=!!pc.xaxis, y=!!pc.yaxis, color=value)) +
      geom_point(alpha=1/5) +
      theme_classic() +
      labs(title=v)
  }

  do.call("grid.arrange", c(plots_variables))
}


#' Plot to show morphology clusters in PC space
#'
#' 'clusterplots' visualizes k-means morphology clusters in PC space
#'
#' @param data is your input data frame
#' @param pc.xaxis is the pc values you want on x-axis (e.g, PC1)
#' @param pc.yaxis is the pce values you want on y-axis (e.g., PC2)
#' @export
clusterplots <- function(data, pc.xaxis, pc.yaxis){
  pc.xaxis <- sym(pc.xaxis)
  pc.yaxis <- sym(pc.yaxis)
  data %>% ggplot(aes(x = !!pc.xaxis, y = !!pc.yaxis, color = as.character(Cluster))) +
    geom_point() +
    stat_ellipse() +
    ggtitle("K-means clusters") +
    labs(color="Cluster") +
    theme_classic()
}

#' Cluster-specific average morphology measures
#'
#' 'clusterfeatures' generates heatmap visualization of average cluster-specific morphology measures relative to other clusters
#'
#' @param data is your input data frame
#' @param featurestart is first column number of morphology measures
#' @param featureend is last column number of morphology measures
#' @export
clusterfeatures <- function(data, featurestart, featureend){
  heatmap <- data %>% group_by(Cluster) %>% summarise(across(featurestart:featureend, ~ mean(.x)))
  #heatmap$`k2$cluster` <- paste0("Cluster ", heatmap$`k2$cluster`)

  heatmap <- column_to_rownames(heatmap, var="Cluster")
  pheatmap(t(heatmap), scale="row", cluster_cols=FALSE, cluster_rows=TRUE,
           border_color=NA, angle_col=45,
           main="Cluster-specific measures")
          #fontsize=12, fontsize_row=12, fontsize_col=12
}

#' What are your morphology cluster percentages across variables of interest?
#'
#' 'clusterpercentage' groups your data by variables of interest then calculates percentage of each morphology cluster within those groups.
#'
#' @param data is your input data frame
#' @param clustercol is the name of your column which contains cluster IDs. Make sure to put this in quotes.
#' @param ... list out your variables of interest without combining (e.g., Cohort, Sex, Treatment)
#' @export
clusterpercentage <- function(data, clustercol,...){
  clustercol <- sym(clustercol)
  data %>%
    group_by(..., !!clustercol) %>%
    count() %>%
    group_by(...) %>%
    mutate(percentage = n/sum(n))
}

#' How do morphology cluster percentages vary across variables of interest?
#'
#' 'clusterpercentage_boxplots' generates boxplots that depict cluster shifts with experimental variables of interest
#'
#' @param data is your input data frame
#' @param ... list out your variables of interest without combining (e.g., Cohort, Sex, Treatment)
#' @export
clusterpercentage_boxplots <- function(data,...){
  y <- unname(sapply(rlang::enexprs(...), as.character))
  z <- paste0(y, collapse="*")
  data %>%
    ggplot(aes(x=Cluster, y=percentage*100, fill=as.character(Cluster))) +
    facet_wrap(as.formula(paste("~",z))) +
    geom_boxplot(position="dodge") +
    geom_point() +
    ggtitle("K-means cluster percentages") +
    labs(fill="Cluster") +
    theme_bw(base_size=14) +
    theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1)) +
    xlab("Cluster") +
    ylab("Percentage")
}

#' Stats analysis: individual morphology measures
#'
#' Linear mixed model to statistically assess how your experimental variables of interest
#' influence each morphology measure, at the animal level
#'
#' The stats_morphologymeasures.animal function fits a linear model using the `lm` function
#' for each morphology measure individually within your dataset. Posthocs are run for each morphology measure
#' individually and bound together into the final dataframe that is output by this function.
#'
#' @param data is your input data frame
#' @param model is your linear mixed model (e.g., Value ~ Treatment*Sex + (1|MouseID))
#' @param posthoc1 is your posthoc comparisons (e.g., when considering sex: ~Treatment|Sex)
#' @param posthoc2 is your posthoc comparisons (e.g., when not considering sex: ~Treatment)
#' @param adjust is your method of multiple test correction (from `emmeans` package: "tukey","scheffe","sidak","dunnettx","mvt","holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr","none"). See "P-value adjustments" section under ?emmeans::summary.emmGrid for more information.
#' @export
### function begins here
stats_morphologymeasures.animal <- function(data,model,posthoc1,posthoc2,adjust){

  y.model <- as.character(model)
  z.model <- as.character(posthoc1)
  a.model <- as.character(posthoc2)

  measure <- unique(as.character(data$Measure))
  log_ggqqplots = list()
  final.output = list()

  # for writing out results to
  anova.out <- NULL
  posthoc.out <- NULL
  posthoc.out2 <- NULL
  levene_out <- NULL
  shapiro_out <- NULL

  for(m in measure){

    tmp <- data %>% filter(Measure == m)
    tmp <- as.data.frame(tmp)

    print(m)

    # linear mixed effects model
    # summary(model) gives you
    options(contrasts=c("contr.sum","contr.poly"))
    model <- lm(as.formula(paste(y.model)), data=tmp)

    ### Test ANOVA assumptions
    # visual check of distribution
    log_ggqqplots[[m]] =
      ggqqplot(residuals(model)) +
      labs(title=m)

    # anova
    #anova = anova(model, test="F", type="III")
    anova = car::Anova(model)
    anova$measure <- paste(m)

    # shapiro test for normality of residuals: passes if p>0.05
    shapiro <- shapiro_test(residuals(model))
    if(shapiro$p.value > 0.05) {
      shapiro$pass <- c("pass")
    } else {
      shapiro$pass <- c("fail")
    }

    shapiro$measure <- paste(m)

    # levene test for homogeneity of variances: passes if p>0.05
    tryCatch(
      {
        levene <- tmp %>% levene_test(as.formula(paste(y.model)))
        if (levene$p > 0.05) {
          levene$pass <- c("pass")
        }
        else {
          levene$pass <- c("fail")
        }
        levene$measure <- paste(m)
        levene_out <- rbind(levene_out, levene)
      },
      error=function(e){
        message('An error occurred with the Levene test: it can only handle interaction terms in model. Function will return an empty dataframe in output[[5]]. All other analyses should run as expected (e.g., Anova, posthocs, Shapiro test, etc.)')
        print(e)
      }
    )

    #posthocs test 1
    refgrid <- ref.grid(model)
    ph <- emmeans(refgrid, as.formula(paste(z.model)))
    ph2 <- contrast(ph, method="pairwise", adjust="none")
    ph3 <- test(ph2, by=NULL, adjust=adjust)
    outsum <- as.data.frame(ph3)
    outsum$measure <- paste(m)

    #posthocs test 2
    refgrid <- ref.grid(model)
    ph <- emmeans(refgrid, as.formula(paste(a.model)))
    ph2 <- contrast(ph, method="pairwise", adjust="none")
    ph3 <- test(ph2, by=NULL, adjust=adjust)
    outsum2 <- as.data.frame(ph3)
    outsum2$measure <- paste(m)

    # save everything before looping to next subregion
    anova.out <- rbind(anova.out,anova)
    posthoc.out <- rbind(posthoc.out,outsum)
    posthoc.out2 <- rbind(posthoc.out2,outsum2)
    shapiro_out <- rbind(shapiro_out, shapiro)
  }

  anova.out$Significant <- ifelse(anova.out$`Pr(>F)` < 0.05, "significant", "ns")
  posthoc.out$Significant <- ifelse(posthoc.out$`p.value` < 0.05, "significant", "ns")
  posthoc.out2$Significant <- ifelse(posthoc.out2$`p.value` < 0.05, "significant", "ns")

  final.output[[1]] = as.data.frame(anova.out)
  final.output[[2]] = posthoc.out
  final.output[[3]] = posthoc.out2
  #final.output[[4]] = do.call("grid.arrange", c(log_ggqqplots, ncol=8))
  final.output[[4]] = log_ggqqplots
  final.output[[5]] = as.data.frame(levene_out)
  final.output[[6]] = as.data.frame(shapiro_out)
  final.output[[7]] = print(model)

  final.output
}

#' Stats analysis: cluster-level changes
#'
#' Linear mixed model to statistically assess how your experimental variables of interest
#' influence cluster percentages, at animal level.
#'
#' The stats_cluster.animal function fits a generalized linear mixed model on your dataset
#' to a beta distribution, which is suitable for values like percentages or probabilities
#' that are constrained to a range of 0-1, using the `glmmTMB` package. Part of the output
#' includes a check of the model fit using the `DHARMa` package, which "uses a simulation-based
#' approach to create readily interpretable scaled (quantile) residuals for fitted (generalized)
#' linear mixed models." The function creates two `DHARMa` plots, contained in output[[4]].
#' You can read more about how to interpret model fit using `DHARMa` by reading the
#' package [vignette](https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html).
#'
#' @param data is your input data frame
#' @param model is your linear mixed model (e.g., Value ~ Cluster*Treatment + (1|MouseID))
#' @param posthoc1 is your first set of posthoc comparisons (e.g., ~Cluster|Treatment)
#' @param posthoc2 is your second set of posthoc comparisons (e.g., ~Cluster)
#' @param adjust is your method of multiple test correction (from `emmeans` package: "tukey","scheffe","sidak","dunnettx","mvt","holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr","none"). See "P-value adjustments" section under ?emmeans::summary.emmGrid for more information.
#' @export
stats_cluster.animal <- function(data, model, posthoc1, posthoc2, adjust){

  y.model <- as.character(model)
  z.model <- as.character(posthoc1)
  a.model <- as.character(posthoc2)

  #modelcheck = NULL
  final.output = list()

  # linear mixed effects model
  options(contrasts=c("contr.sum","contr.poly"))
  #model <- lmerTest::lmer(as.formula(paste(y.model)), data=data)
  #model <- glmmadmb(as.formula(paste(y.model)), data=data, family="beta")
  model <- glmmTMB(as.formula(paste(y.model)), data=data, family=beta_family(link="logit"))

  ### Test ANOVA assumptions
  # visual check of distribution
  # log_ggqqplots = ggqqplot(residuals(model))
  res = simulateResiduals(model)
  plot(res, rank = T)
  modelcheck <- recordPlot()
  plot.new()
  modelcheck

  # anova
  anova = car::Anova(model)
  anova

  # posthocs 1: considering sex
  ph <- emmeans(model, as.formula(paste(z.model)))
  ph2 <- contrast(ph, method="pairwise", adjust="none")
  ph3 <- test(ph2, by=NULL, adjust=adjust)
  posthoc1 <- as.data.frame(ph3)
  posthoc1

  # posthocs 2: not considering sex
  ph <- emmeans(model, as.formula(paste(a.model)))
  ph2 <- contrast(ph, method="pairwise", adjust="none")
  ph3 <- test(ph2, by=NULL, adjust=adjust)
  posthoc2 <- as.data.frame(ph3)
  posthoc2

  # posthocs not looking within cluster
  #test(pairs(ph), by=NULL, adjust="dunnet")

  # annotate as significant for easy filtering
  #anova$Significant <- ifelse(anova$`Pr(>F)` < 0.05, "significant", "ns")
  posthoc1$Significant <- ifelse(posthoc1$`p.value` < 0.05, "significant", "ns")
  posthoc2$Significant <- ifelse(posthoc2$`p.value` < 0.05, "significant", "ns")

  final.output[[1]] = anova
  final.output[[2]] = posthoc1
  final.output[[3]] = posthoc2
  #final.output[[4]] = do.call("grid.arrange", c(log_ggqqplots, ncol=8))
  final.output[[4]] = modelcheck
  final.output[[5]] = print(model)

  final.output
}

#' Data visualization: individual morphology measures
#'
#' Linear mixed model to statistically assess how your experimental variables of interest
#' influence each morphology measure, at the cell level
#'
#' @param data is your input data frame
#' @param group is the variable you want to group your data by
#' @export
cell.level_boxplots <- function(data,group){
  group <- sym(group)
  data %>%
    ggplot(aes(x=Measure, y=Value, fill=!!group)) +
    facet_wrap(~Measure, scales="free") +
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), trim=FALSE) +
    xlab("Morphology Measure") +
    ylab("Value") +
    ggtitle("Individual morphology measure changes at cell level")
}


