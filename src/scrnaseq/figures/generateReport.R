#' Sum of vector elements
#'
#' `sum` returns the sum of all the values present in its arguments.
#'
#' This is a generic function: methods can be defined for it directly
#' or via the [Summary()] group generic. For this to work properly,
#' the arguments `...` should be unnamed, and dispatch is on the
#' first argument.
#' --== Generate Complete Report ==--
generateReport <- function(seuratObject, saveFolder=NULL) {

  .validateSeuratObject(seuratObject)


  if (is.null(saveFolder)) {
    saveFolder = file.path(CONFIG$resultsHome, CONFIG$timestamp)
  }


  # Folder 1: Preprocessing  
  visualize_batch_effect(seuratObject, saveFolder=saveFolder)
  
  visualize_umap_clusters(seuratObject, saveFolder=saveFolder)
  

  # Folder 2: Differential Gene Expression




  # Folder 3: Gene Network Analysis

}