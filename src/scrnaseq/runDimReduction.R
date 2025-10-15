#' --== Run Dimension Reduction ==--
#'
#' @param seuratObject Seurat object to analyze
#' @param reduction.target Existing reduction to use (e.g., "pca"). If NULL, runs PCA
#' @param assay Assay to use. If NULL, uses DefaultAssay
#' @param numPCs Number of PCs to compute. If NULL, uses CONFIG$numPCs
#' @param dimsUsed Dimensions to use for UMAP/tSNE/clustering. If NULL, uses 1:numPCs
#' @param reduction.output Base name for output reductions (e.g., "integrated" creates "integrated", "integrated.umap", etc.)
#' @param verbose Print progress messages
#'
runDimReduction <- function(seuratObject, reduction.target = NULL, assay = NULL,
                            numPCs = NULL, dimsUsed = NULL,
                            reduction.output = "pca", verbose = TRUE) {
  
  if (verbose) cat("Beginning Dimension Reduction Analysis.\n")
  
  # Set defaults
  if (is.null(assay)) assay <- DefaultAssay(seuratObject)
  if (is.null(numPCs)) numPCs <- CONFIG$numPCs
  if (is.null(dimsUsed)) dimsUsed <- 1:numPCs
  
  if (verbose) cat(sprintf("Using assay: %s\n", assay))
  
  # Run PCA if no existing dimension reduction target is provided
  if (is.null(reduction.target)) {
    if (verbose) cat("No existing dimension reduction target set.\n")
    
    npcs <- numPCs
    min_pcs <- 5
    pca_success <- FALSE
    
    for (attempt in 1:5) {
      if (verbose) cat(sprintf("Attempt %d: Identifying %d Principal Components.\n", 
                               attempt, npcs))
      
      result <- tryCatch({
        seuratObject <- RunPCA(
          seuratObject, 
          npcs = npcs, 
          assay = assay,
          features = VariableFeatures(seuratObject, assay = assay),
          reduction.name = paste0("pca.", reduction.output),
          verbose = FALSE
        )
        if (verbose) cat("RunPCA succeeded.\n")
        TRUE
      }, error = function(e) {
        if (verbose) {
          cat(sprintf("RunPCA failed with %d PCs: %s\n", npcs, conditionMessage(e)))
        }
        npcs <<- max(min_pcs, floor(npcs * 0.75))
        if (verbose) cat(sprintf("Reducing to %d PCs...\n", npcs))
        FALSE
      })
      
      if (result) {
        pca_success <- TRUE
        # Update dimsUsed if PCs were reduced
        if (npcs < numPCs) {
          dimsUsed <- 1:npcs
          if (verbose) cat(sprintf("Updated dimensions used to 1:%d\n", npcs))
        }
        break
      }
    }
    
    if (!pca_success) {
      stop("RunPCA failed after 5 attempts")
    }
    
    reduction.target <- paste0("pca.", reduction.output)
  }
  
  # Run downstream analyses
  if (verbose) cat("Running UMAP...\n")
  seuratObject <- RunUMAP(
    seuratObject, 
    dims = dimsUsed,
    reduction = reduction.target,
    reduction.name = paste0("umap.", reduction.output),
    verbose = FALSE
  )
  
  if (verbose) cat("Running tSNE...\n")
  seuratObject <- RunTSNE(
    seuratObject, 
    dims = dimsUsed, 
    check_duplicates = FALSE,
    reduction = reduction.target,
    reduction.name = paste0("tsne.", reduction.output),
    verbose = FALSE
  )
  
  if (verbose) cat("Finding neighbors...\n")
  seuratObject <- FindNeighbors(
    seuratObject, 
    dims = dimsUsed,
    reduction = reduction.target,
    verbose = FALSE
  )
  
  if (verbose) cat("Finding clusters...\n")
  seuratObject <- FindClusters(
    seuratObject,
    resolution = CONFIG$clusterResolution,
    cluster.name = paste0("clusters.", reduction.output),
    verbose = FALSE
  )
  
  if (verbose) cat("Dimension reduction complete.\n")
  return(seuratObject)
}

