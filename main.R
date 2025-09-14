#' [General]
#' TODO:
#' - Figure out a way to support 
#' 
#' - Create an overall "merge_dataset" function with options to 
#'   include proper integration as well as cell cycle regression
#'   functionality
#'   
#' - Implement the initial versions of interactive or
#'   intermediate diagnostics reports that prompt the user
#'   for how to proceed based on the current state.
#'   (i.e. Impact of cell cycle on initial dimension
#'   reduction components, or the elbow plot of PC variability
#'   when deciding how many dims to include in UMAP/tSNE)
#'   
#' - Utilize `samples@misc` slot to store run specifications
#'   (i.e. whether to use pre- or post-batch effect corrections)
#'  
#' [Generating Results]
#' TODO:
#' - Automate generation of cluster markers for all cluster
#'   sources within the sample. Save entire list as well as 
#'   top-n markers as .csv files.
#'   
#' [Things to store in `@misc`]
#' - Species
#' - Results directory
#'   
#'   
#'   





# --== Initialize Pipeline ==--
CONFIG <- list(
  # [File Path Locations]
  "srcDir" = "./src",
  "resultsHome" = "./results",
  
  # Used to access genome references and marker gene file
  "dataHomeLocal" = "/media/michael/Nguyen-Lab/sequencing-data",
  "dataHomeRemote" = "../sequencing-data",
  
  "defaultExclude" = c("meta", "Raw_FASTQ", "all-sample", "process"),
  
  # [Analysis Metrics]
  "useFilter" = TRUE,
  
  "minFeatRNA" = 200,
  "mtExprLimit" = 20,        # The highest % of MT gene reads allowed
  
  "numPCs" = 30,             # Amount of principal components to generate
  
  "clusterResolution" = 1
)

ref <- celldex::HumanPrimaryCellAtlasData()

source(paste(CONFIG$srcDir, "util.R", sep="/"))
.initialize_env(remote=FALSE)

MARKERS$target <- c("FOXK1")


# --== Run Pipeline ==--
SAMPLE_DATA <- list(
  "datasetName" = "GSE175604",
  "sampleDir" = "/nguyen/GE8070/Deep/output/combined_human/",
  "species" = "human", # Options: ("human", "mouse")
  "dataFormat" = "parsebio" # Options: ("10X", "counts", "rdata", "parsebio")
)

SAMPLE_DATA <- list(
  "datasetName" = "GUT_DSS_MOUSE",
  "sampleDir" = "/nguyen/GE8070/Deep/output/combined_mouse/",
  "species" = "mouse", # Options: ("human", "mouse")
  "dataFormat" = "parsebio" # Options: ("10X", "counts", "rdata", "parsebio")
)

mouse_samples <- c("14_GUT_DSS_FK1FL_2", "15_GUT_DSS_FK1FL_3",
                   "16_GUT_DSS_FK1-FOXP3_4", "17_GUT_DSS_FK1-FOXP3_5")

samples <- load_dataset(SAMPLE_DATA, includeList=mouse_samples)

samples <- perform_integration(samples)

samples <- run_dim_reduction(samples, reduction.target="orig.pca")

visualize_batch_effect(samples)








