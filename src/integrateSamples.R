
integrateSamples <- function(seuratObject) {
  seuratObject[["RNA"]] <-
    split(seuratObject[["RNA"]], f = seuratObject$orig.sample)
  
  seuratObject <- IntegrateLayers(seuratObject,
                                method = CCAIntegration,
                                orig.reduction = "orig.pca",
                                new.reduction = "integrated.cca",
                                verbose = FALSE
  )

  seuratObject[["RNA"]] <- JoinLayers(seuratObject[["RNA"]])
  
  # Second Dimension Reduction Attempt
  seuratObject <-FindNeighbors(seuratObject, reduction = "integrated.cca", dims = 1:30)
  seuratObject <- FindClusters(seuratObject, resolution = 0.5)
  seuratObject <- RunUMAP(seuratObject,
                        dims = 1:20,
                        reduction = "integrated.cca",
                        reduction.name = "umap")
  seuratObject <- RunTSNE(seuratObject,
                          dims = 1:20,
                          check_duplicates = FALSE,
                          reduction = "integrated.cca",
                          reduction.name = "tsne")
  
  return(seuratObject)
}

