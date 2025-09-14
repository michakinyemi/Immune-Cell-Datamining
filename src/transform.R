

#' --== Run Dimension Reduction ==--
#' 
#' 
run_dim_reduction <- function(seuratObject, reduction.target=NULL,
                              reduction.name="orig.pca",
                              numPCs=NULL, dimsUsed=NULL) {
  
  cat("Beginning Dimension Reduction Analysis.")
  
  if (is.null(numPCs)) { numPCs = CONFIG$numPCs }
  if (is.null(dimsUsed)) { dimsUsed = 1:CONFIG$numPCs }
  
  if (is.null(reduction.target)) {
    cat("No existing dimension reduction target set.")
    cat(sprintf("Identifying %d Principal Components.", CONFIG$numPCs))
    
    seuratObject <- RunPCA(seuratObject, npcs=numPCs,
                           features=VariableFeatures(seuratObject),
                           reduction.name=reduction.name)
  }
  
  # TODO: Potential interactive dimension selection mode where
  # ElbowPlot & JackStrawPlot/Score is displayed to user
  

  seuratObject <- FindNeighbors(seuratObject,
                                reduction=reduction.name,
                                dims=dimsUsed)
  
  seuratObject <- FindClusters(seuratObject,
                               resolution=CONFIG$clusterResolution,
                               cluster.name=paste(reduction.name, ".clusters", sep=""))

  seuratObject <- RunUMAP(seuratObject,
                          dims=dimsUsed,
                          reduction=reduction.name,
                          reduction.name=paste(reduction.name, ".umap", sep=""))
    
  seuratObject <- RunTSNE(seuratObject,
                          dims=dimsUsed,
                          check_duplicates=FALSE,
                          reduction=reduction.name,
                          reduction.name=paste(reduction.name, ".tsne", sep=""))
  
  return(seuratObject)
}



#' --== Perform Integration ==--
#' Integrate a list of Seurat objects using either CCA or SCT methods.
#' 
perform_integration <- function(sampleList, method="CCA",
                                numPCs=NULL, dimsUsed=NULL,
                                cellCycleRegression=TRUE) {
  
  sampleCount = length(sampleList)
  
  print(sprintf("Beginning to Combine Dataset Containing %d Samples", sampleCount))
  
  if (is.null(numPCs)) { numPCs = CONFIG$numPCs }
  if (is.null(dimsUsed)) { dimsUsed = 1:CONFIG$numPCs }
  
  if (method == "CCA") {
  
    seuratObject <- merge(x=sampleList[[1]], y=sampleList[2:sampleCount],
                          project=SAMPLE_DATA$datasetName)
    
    seuratObject <- NormalizeData(seuratObject)
    seuratObject <- FindVariableFeatures(seuratObject)
    seuratObject <- ScaleData(seuratObject)
    seuratObject <- RunPCA(seuratObject, npcs=numPCs, reduction.name="orig.pca")
    
    seuratObject <- IntegrateLayers(seuratObject,
                                    method = CCAIntegration,
                                    orig.reduction = "orig.pca",
                                    new.reduction = "integrated.cca",
                                    verbose = TRUE)
    
    seuratObject[["RNA"]] <- JoinLayers(seuratObject[["RNA"]])
    
    seuratObject <- ScaleData(seuratObject)
    
    # Second Dimension Reduction Attempt
    seuratObject <- FindNeighbors(seuratObject,
    				  reduction="integrated.cca",
                                  dims=dimsUsed)
    
    seuratObject <- FindClusters(seuratObject, resolution=CONFIG$clusterResolution,
                                 cluster.name="integrated.cca.clusters")
    
    seuratObject <- RunUMAP(seuratObject,
                            dims = dimsUsed,
                            reduction = "integrated.cca",
                            reduction.name = "integrated.cca.umap")
    
    seuratObject <- RunTSNE(seuratObject,
                            dims = dimsUsed,
                            check_duplicates = FALSE,
                            reduction = "integrated.cca",
                            reduction.name = "integrated.cca.tsne")
    
  } 
  else if (method == "SCT") {
    for (i in 1:length(sampleList)) {
    
      sampleList[[i]] <- SCTransform(sampleList[[i]], verbose = TRUE)
      
      if (cellCycleRegression) {
        
        sampleList[[i]] <- CellCycleScoring(sampleList[[i]],
                                            s.features = CONFIG[["s.genes"]],
                                            g2m.features = CONFIG[["g2m.genes"]],
                                            set.ident = TRUE, 
                                            assay = "SCT")
      }
    }
   
    features <- SelectIntegrationFeatures(object.list = sampleList, nfeatures = 3000)
    
    sampleList <- PrepSCTIntegration(object.list = sampleList, anchor.features = features)
    
    anchors <- FindIntegrationAnchors(object.list = sampleList,
                                      normalization.method = "SCT",
                                      anchor.features = features)  # or choose variableFeatures(n)
    
    seuratObject <- IntegrateData(anchorset = anchors,
                                  normalization.method = "SCT",
                                  new.assay.name = "integrated.sct")
    
    DefaultAssay(seuratObject) <- "integrated.sct"
    
    
    if (cellCycleRegression) {
      seuratObject <- ScaleData(seuratObject, 
                                features=rownames(seuratObject),
                                vars.to.regress = c("S.Score", "G2M.Score"))
    } else {
      seuratObject <- ScaleData(seuratObject, features=rownames(seuratObject))
    }
    
    seuratObject <- RunPCA(seuratObject, npcs = numPCs, verbose = TRUE,
                           reduction.name="integrated.sct.pca")
    
    seuratObject <- FindNeighbors(seuratObject, reduction = "integrated.sct.pca", dims = dimsUsed)
    seuratObject <- FindClusters(seuratObject, resolution = CONFIG$clusterResolution,
                                 cluster.name = "integrated.sct.clusters")
    
    seuratObject <- RunUMAP(seuratObject, dims = dimsUsed, reduction = "integrated.sct.pca",
                            reduction.name = "integrated.sct.umap")
    
    seuratObject <- RunTSNE(seuratObject, dims = dimsUsed, reduction = "integrated.sct.pca",
                            reduction.name = "integrated.sct.tsne")
    
  } 
  else if (method == "Merge") {
    seuratObject <- merge(x=sampleList[[1]], y=sampleList[2:sampleCount], project=SAMPLE_DATA$dataID)
    
    seuratObject <- NormalizeData(seuratObject)
    
    seuratObject <- FindVariableFeatures(seuratObject)
    
    seuratObject <- ScaleData(seuratObject, features=rownames(seuratObject))
    
    if (joinLayers) {
      seuratObject[["RNA"]] <- JoinLayers(seuratObject[["RNA"]])
    }
    
  }

  return(seuratObject)
}

#' --== Regress Cell Cycle Genes ==--
#' 
#' 
regress_cell_cycle <- function(seuratObject, reduction.target="integrated.cca",
                               reduction.name="regressed.cca",
                               numPCs=NULL, dimsUsed=NULL) {
  
  if (is.null(numPCs)) { numPCs = CONFIG$numPCs }
  if (is.null(dimsUsed)) { dimsUsed = 1:CONFIG$numPCs }
  
  seuratObject <- Seurat::CellCycleScoring(seuratObject,
                                           s.features = CONFIG[["s.genes"]],
                                           g2m.features = CONFIG[["g2m.genes"]],
                                           set.ident = FALSE)
  
  # seuratObject <- Seurat::ScaleData(seuratObject,
  #                                   features=rownames(seuratObject),
  #                                   vars.to.regress=c("S.Score", "G2M.Score"))
  
  embedMatrix <- Embeddings(seuratObject, reduction=reduction.target)
  
  X <- as.matrix(seuratObject@meta.data[, c("S.Score", "G2M.Score")])
  
  # Regress out contributions of S.Score and
  # G2M.Score from the reduction's embeddings
  residualEmbed <- apply(embedMatrix, 2, function(pc) {
    fit <- lm(pc ~ X)
    return(resid(fit))
  })
  
  rownames(residualEmbed) <- rownames(embedMatrix)
  
  seuratObject[[reduction.name]] <- CreateDimReducObject(
    embeddings = residualEmbed,
    key="RegressPC_",
    assay=DefaultAssay(seuratObject)
  )
  

  seuratObject <- Seurat::FindNeighbors(seuratObject, reduction=reduction.name,
                                        dims=dimsUsed)
  
  seuratObject <- Seurat::FindClusters(seuratObject, resolution=CONFIG$clusterRes,
                                       cluster.name=paste(reduction.name, ".clusters", sep=""))

  seuratObject <- Seurat::RunUMAP(seuratObject,
                                  dims=dimsUsed,
                                  reduction=reduction.name,
                                  reduction.name=paste(reduction.name, ".umap", sep=""))

  seuratObject <- Seurat::RunTSNE(seuratObject,
                                  dims=dimsUsed,
                                  check_duplicates=FALSE,
                                  reduction=reduction.name,
                                  reduction.name=paste(reduction.name, ".tsne", sep=""))
  
  return(seuratObject)

}
