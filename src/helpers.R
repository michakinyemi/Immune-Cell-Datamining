

# --== convertFeatures() ==--
# Input/Usage
#

convertToSymbol <- function(geneList, species="Hs") {
  #features <- data@Dimnames[[1]]
  
  if (species == "Mm") database = org.Mm.eg.db
  else database = org.Hs.eg.db
  
  print(typeof(database))
  
  geneTable <- mapIds(
    org.Mm.eg.db,
    keys = geneList, #rownames(data)
    column = "SYMBOL",
    keytype = "ENSEMBL"
  )

  inCommon <- geneList %in% names(geneTable)
  geneList[inCommon] <-
    unlist(geneTable[geneList[inCommon]])
  
  geneList <- toupper(geneList)
  
  return(geneList)
}


# --== annotateCellTypes() ==--
#

annotateCellTypes <- function(
    seuratObject,
    markerGenes=NULL) {
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
