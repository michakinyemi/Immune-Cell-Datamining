# Input/Usage
#   datasetID             Directory containing samples (10X Genomics Format)from a dataset
#   samples (optional)    List of samples to load (defaults to all samples)

# Output
#   aaa

library(Seurat)
library(yaml)
library(ggThemeAssist)
library(dplyr)
library(data.table)
library(biomaRt)

library(org.Hs.eg.db)
library(org.Mm.eg.db)


# Update to use the appropriate GSE identifier for your analysis
dataID = "Treg"

# TODO: Support integration of separate datasets
#(while retaining sample characteristics/identifiers)

genesOfInterest = c("TFAM",
                    "FOXK1",
                    "FOXK2",
                    "TFEB",
                    "LDHA",
                    "PCSK9",
                    "TFAM",
                    "GPR84")


# Filtering Metrics
FILT_METRICS = list(minFeatRNA = 500,
                    maxFeatRNA = 800,
                    mtExprLimit = 10) # The highest % of MT gene reads allowed

NUM_PCS = 50 # The # of principal components to generate



CELL_TYPES = c("metastatic", "primary")



# --== Convert Gene Names ==--
# ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
# gene_info <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
#                    filters = "ensembl_gene_id",
#                    values = features.tsv[1],
#                    mart = ensembl,
#                    uniqueRows = FALSE)





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




# --== mergeData() ==--
mergeData <- function(sampleList) {
  sampleCount = length(sampleList)
  
  #idList = list()
  
  for (i in 1:sampleCount) {
    sampleList[[i]] <- filterSample(sampleList[[i]])
    
    sampleList[[i]] <- NormalizeData(sampleList[[i]])
    sampleList[[i]] <- ScaleData(sampleList[[i]])
    #idList <- append(idList, paste("s", i, sep=""))
  }
  
  mergedData <-
    merge(x = sampleList[[1]], y = sampleList[2:sampleCount], project = dataID)
  return(mergedData)
  mergedData <- JoinLayers(mergedData)
  return(mergedData)
}


# --== integrateData() ==--
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






# --== runDimReduction() ==--
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





# --== annotateCellTypes() ==--
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
