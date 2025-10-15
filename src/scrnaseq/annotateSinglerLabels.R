#' --== Annotate Cells Using SingleR ==--
#' 
#' 
annotateSinglerLabels <- function(seuratObject, species="human", assay="SCT") {

  if (species == "human") {
    ref <- REFERENCES$celldex.human
  } else if (species == "mouse") {
    ref <- REFERENCES$celldex.mouse
  }

  sce <- as.SingleCellExperiment(seuratObject, assay = assay)


  # Main Cell Type Labeling
  singler_results_main <- SingleR(
    test = sce,
    ref = ref,
    labels = ref$label.main
  )

  seuratObject[[paste0("SingleR_main.", tolower(assay))]] <- singler_results_main$labels

  # Fine Cell Type Labeling
  singler_results_fine <- SingleR(
    test = sce,
    ref = ref,
    labels = ref$label.fine
  )
  
  seuratObject[[paste0("SingleR_fine.", tolower(assay))]] <- singler_results_fine$labels
  
  return(seuratObject)
}
