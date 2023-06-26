#' Tidying up your FracLac output
#' 
#' This function allows you to load in FracLac output and clean up names of cells
#' 
#' @param x is your directory containing the 'Hull and Circles Results.txt' file
#' @param y is the variable name for your processed output containing your tidied fraclac data 

fraclac_tidying <- function(x){
  setwd(x)
  FracLac1 <- read.csv("Hull and Circle Results.txt", sep='\t')
  y <- FracLac1 %>%
  separate(1, into=c("Name","trash"), sep="_thresholded-roitif_") %>%
  separate(trash, into=c("ID","trash1"), sep="tif") %>%
  select(-trash1)
  y
}

#' Concatenating your individual AnalyzeSkeleton output files into one dataframe
#' 
#' This function allows you to load in all of your AnalyzeSkeleton output and clean up IDs to match FracLac data for downstream merging
#' 
#' @param x is your directory containing all of your individual AnalyzeSkeleton .csv result files
#' @param y is the variable name for your processed output containing your tidied analyzeskeleton data 

skeleton_tidying <- function(x){
  setwd(x)
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
    separate(Name, into=c("Name","trash"), sep="_thresholded-roi.tif_") %>%
    separate(trash, into=c("ID","junk2"), sep=".tif_results") %>% select(-junk2) %>%
    unite("UniqueID",c("Name","ID"),sep="_",remove=FALSE)
  y
}

#' Merge fraclac and skeleton results together
#' 
#' This function allows you to merge your fraclac and skeleton results by Name and ID, get rid of non-numerical data, and brush up feature names
#' 
#' @param x is your tidied fraclac output
#' @param y is your tidied skeleton output

merge_data <- function(x,y){
  data <- inner_join(x, y, by=c("Name","ID"))
  
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


#' Would you like to get some metadata?
#'
#' This function allows you to extract metadata columns from the UniqueID column.
#' 
#' @param x is your tidied dataframe containing all of your numerical morphology data
#' @param y is your metadata listed out in order that they appear in 'Name' column: e.g., c("Cohort","Sex","MouseID")
#' @param z is the character, "_" or "-", that separates your metadata
metadata_columns <- function(x,y,z){
  final <- x %>% separate(Name, into=y, sep=z)
  final
}

#' Exploratory data analysis: outlier detection
#'
#' This function allows you to use boxplots of features to look for those with broader ranges which might dominate the data analysis
#' 
#' @param x is your final dataframe which contains gathered data (measure, value) format
outliers_boxplots <- function(x){
  x %>% 
    ggplot(aes(x=measure, y=value)) +
    geom_boxplot() +
    geom_jitter(width=0.15)+
    theme(axis.text.x = element_text(angle = 90))
}

#' Exploratory data analysis: outlier detection
#'
#' This function allows you to check distributions of features to determine if and how to normalize
#' 
#' @param x is your final dataframe which contains gathered data (measure, value) format
outliers_distributions <- function(x){
  x %>%
    ggplot(aes(x = value)) +
    geom_histogram(bins = 40) +
    facet_wrap(~measure, scales = "free_x", ncol=8)
}

#' Exploratory data analysis: normalization methods
#'
#' This function allows you to check distributions of features after log normalization (log10(x+1))
#' 
#' @param x is your final dataframe which contains gathered data (measure, value) format
normalize_logplots <- function(x){
  x %>% 
    ggplot(aes(x = log(value+1))) +
    geom_histogram(bins = 40) +
    facet_wrap(~measure, scales = "free", ncol=7) +
    theme(strip.text.x = element_text(size=8)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size=8)) +
    scale_x_continuous(name = "Morphology measure (log(value+1))") +
    scale_y_continuous(name = "Count")
}



