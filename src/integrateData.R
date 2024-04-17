
# --== integrateData() ==--
# test

integrateData <- function(sampleList) {
  sampleCount = length(sampleList)
  
  for (i in 1:sampleCount) {
    sampleList[[i]] <- filterSample(sampleList[[i]])
    
    sampleList[[i]] <- NormalizeData(sampleList[[i]])
    sampleList[[i]] <- ScaleData(sampleList[[i]])
  }
  
  mergedData <-
    merge(x = sampleList[[1]], y = sampleList[2:sampleCount], project = dataID)
  
  mergedData <- JoinLayers(mergedData)
  
  mergedData <- FindVariableFeatures(mergedData)
  mergedData <- RunPCA(mergedData, npcs = NUM_PCS)
  
  mergedData[["RNA"]] <-
    split(mergedData[["RNA"]], f = mergedData$sample)
  
  mergedData <-
    IntegrateLayers(
      object = mergedData,
      method = CCAIntegration,
      orig.reduction = "pca",
      new.reduction = "integrated.cca",
      verbose = FALSE
    )
  
  mergedData[["RNA"]] <- JoinLayers(mergedData[["RNA"]])
  
  mergedData <-
    FindNeighbors(mergedData, reduction = "integrated.cca", dims = 1:30)
  mergedData <- FindClusters(mergedData, resolution = 0.5)
  mergedData <-
    RunUMAP(
      mergedData,
      dims = 3:20,
      reduction = "integrated.cca",
      reduction.name = "umap"
    )
  
  return(mergedData)
}