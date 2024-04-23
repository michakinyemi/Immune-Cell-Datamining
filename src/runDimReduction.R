runDimReduction <- function(seuratObject, pcaFeats=NULL) {
  

  # seuratObject <- RunPCA(seuratObject,
  #                        npcs=2,
  #                        reduction.name="pca.treg",
  #                        features=c("FOXP3","CCL3"))
  

  seuratObject <- FindVariableFeatures(seuratObject)
  seuratObject <- RunPCA(seuratObject,
                         npcs=NUM_PCS)
  
  seuratObject <- FindVariableFeatures(seuratObject)
  seuratObject <- RunPCA(seuratObject, npcs = NUM_PCS, reduction.name = "orig.pca")
  seuratObject <-FindNeighbors(seuratObject, reduction = "orig.pca",
                             dims = 1:30)
  
  seuratObject <- FindClusters(seuratObject, resolution = 0.5)
  seuratObject <- RunUMAP(seuratObject,
                        dims = 1:20,
                        reduction = "orig.pca",
                        reduction.name = "orig.umap")
  seuratObject <- RunTSNE(seuratObject,
                          dims = 1:20,
                          check_duplicates = FALSE,
                          reduction.name = "orig.tsne")
  
  return(seuratObject)
}



