#' --== Manually Annotate Clusters ==--
#' 
#' 
manualAnnotation <- function(seuratObject, markerGenes=NULL, useSingleR=FALSE) {
  # Initialize a new column for cell type annotations
  seuratObject$cell_type <- NA
  
  # Loop through each cluster and annotate cells
  for (cell_type in names(markerGenes)) {

    # Identify cells expressing marker genes for the current cell type
    marker_genes_expr <- markerGenes[[cell_type]]
    cell_indices <-
      WhichCells(seuratObject, expression = marker_genes_expr)
    
    # Assign the cell type label to the corresponding cells
    seuratObject$cell_type[cell_indices] <- cell_type
  }
  
  # Return the annotated Seurat object
  Idents(seuratObject) <- "cell_type"

  return(seuratObject)
}
