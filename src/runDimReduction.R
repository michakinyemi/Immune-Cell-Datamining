
# --== runDimReduction() ==--
# a

runDimReduction <- function(seuratObject) {
  seuratObject <- FindVariableFeatures(seuratObject)
  
  seuratObject <- RunPCA(seuratObject, npcs = NUM_PCS)
  
  seuratObject <- FindNeighbors(seuratObject, dims = 3:20)
  
  seuratObject <- FindClusters(seuratObject, resolution = 1.0)
  # Can raise resolution value to produce more clusters
  # it's more tedious, but adds more precision to the cell type groupings
  
  seuratObject <-
    RunUMAP(
      seuratObject,
      reduction = "pca",
      dims = 3:20,
      reduction.name = "umap"
    )
  
  seuratObject <-
    RunTSNE(seuratObject,
            dims = 3:20,
            check_duplicates = FALSE)
  
  return(seuratObject)
}