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
#' [.init_env vs .update_env]
#' - Test
#'

# --== Initialize Pipeline ==--
library(yaml)
CONFIG <- yaml::read_yaml("config.yml")

source(paste(CONFIG$srcDir, "env.R", sep="/"))
.init_env(remote=FALSE)

