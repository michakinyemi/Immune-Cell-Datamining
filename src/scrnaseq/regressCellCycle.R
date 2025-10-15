#' --== Regress Cell Cycle Genes ==--
#' 
#' 
regressCellCycle <- function(seuratObject, reduction.target="integrated.cca",
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



# # store mitochondrial percentage in object meta data
# pbmc <- PercentageFeatureSet(pbmc, pattern = "^MT-", col.name = "percent.mt")

# # run sctransform
# pbmc <- SCTransform(pbmc, vars.to.regress = "percent.mt", verbose = FALSE)




