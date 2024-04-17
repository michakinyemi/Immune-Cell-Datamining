

# --== convertFeatures() ==--
# Input/Usage
#

convertToSymbol <- function(data) {
  features <- data@Dimnames[[1]]
  geneTable <- mapIds(
    org.Mm.eg.db,
    keys = rownames(data),
    column = "SYMBOL",
    keytype = "ENSEMBL"
  )
  
  inCommon <- rownames(data) %in% names(geneTable)
  rownames(data)[inCommon] <-
    unlist(geneTable[rownames(data)[inCommon]])
  
  return(data)
}


# --== annotateCellTypes() ==--
#

annotateCellTypes <- function(seuratObject, markerGenes) {
  # Initialize a new column for cell type annotations
  seurat_object$cell_type <- NA
  
  # Loop through each cell type and annotate cells
  for (cell_type in names(markerGenes)) {
    # Identify cells expressing marker genes for the current cell type
    marker_genes_expr <- markerGenes[[cell_type]]
    cell_indices <-
      WhichCells(seurat_object, expression = marker_genes_expr)
    
    # Assign the cell type label to the corresponding cells
    seurat_object$cell_type[cell_indices] <- cell_type
  }
  
  # Return the annotated Seurat object
  return(seurat_object)
}
