#' Integrate Seurat objects to remove batch effects and regress cell cycle variations.
#' Recommended for large datasets (~150k cells): Use method="RPCA".
#' Regresses cell cycle by default; handles mitochondrial if enabled.
integrateSamples <- function(sampleList, method = "RPCA",  # RPCA for scalability/large data
                             numPCs = NULL, dimsUsed = NULL,
                             cellCycleRegression = TRUE,  # Enabled for cell cycle regression
                             mtRegression = FALSE) {               # Optional: Regress mitochondrial % 
  
  sampleCount <- length(sampleList)
  if (sampleCount == 0) stop("No samples provided.")
  
  regressionVars <- c()
  if (cellCycleRegression) regressionVars <- c(regressionVars, "S.Score", "G2M.Score")
  if (mtRegression) regressionVars <- c(regressionVars, "percent.mt")
  
  print(sprintf("Beginning integration of %d samples (%s method; ~%dk cells total)", sampleCount, method, round(sum(sapply(sampleList, ncol))/1e3, 0)))
  if (length(regressionVars) > 0) print(sprintf("Regressing variables: %s", paste(regressionVars, collapse = ", ")))
  
  if (is.null(numPCs)) numPCs <- CONFIG$numPCs
  if (is.null(dimsUsed)) dimsUsed <- 1:numPCs
  
  if (method == "RPCA") {
    
    # Perform per-sample normalization, cell cycle scoring, regression, and PCA
    print("Running per-sample normalization, PCA, and cell cycle scoring for RPCA...")
    for (i in seq_along(sampleList)) {
      seuratObject <- sampleList[[i]]
      cat(sprintf("Processing Sample %d/%d: %s\n", i, sampleCount, seuratObject@project.name))
      
      seuratObject <- NormalizeData(seuratObject)
      seuratObject <- FindVariableFeatures(seuratObject)
      
      # Perform cell cycle scoring if regressing it (needed before ScaleData)
      if (cellCycleRegression) {
        seuratObject <- CellCycleScoring(seuratObject,
                                         s.features = CONFIG[["s.genes"]],
                                         g2m.features = CONFIG[["g2m.genes"]],
                                         set.ident = FALSE)
      }
      
      # Regress in ScaleData (includes cell cycle if enabled)
      seuratObject <- ScaleData(seuratObject, vars.to.regress = regressionVars, verbose = TRUE)
      seuratObject <- RunPCA(seuratObject, npcs = numPCs, reduction.name = "pca", verbose = TRUE)
      
      sampleList[[i]] <- seuratObject
    }

    print("Performing batch correction with RPCA...")
    # Integrate directly from the list (no pre-merging, to preserve per-sample "pca" reductions)
    seuratObject <- IntegrateLayers(sampleList,
                                    method = RPCAIntegration,
                                    orig.reduction = "pca",
                                    new.reduction = "integrated.rpca",
                                    k.anchor = kAnchor,
                                    k.score = kScore,
                                    verbose = TRUE)

    # Post-integration processing
    seuratObject[["RNA"]] <- JoinLayers(seuratObject[["RNA"]])
    
    seuratObject <- runDimReduction(seuratObject = seuratObject,
                                    reduction.target = "integrated.rpca",
                                    reduction.output = "integrated")
    
  } else if (method == "SCT") {
    # SCT mode with cell cycle regression and RPCA tweaks for large data
    for (i in seq_along(sampleList)) {
      seuratObject <- sampleList[[i]]
      cat(sprintf("Processing Sample %d/%d: %s (SCT & cell cycle scoring)\n", i, sampleCount, seuratObject@project.name))

      seuratObject <- NormalizeData(seuratObject, verbose = TRUE)  # Log-normalize first for % mito

      seuratObject <- CellCycleScoring(seuratObject,
                                  s.features = CONFIG[["s.genes"]],
                                  g2m.features = CONFIG[["g2m.genes"]],
                                  set.ident = FALSE,
                                  nbin=16)
      
      # SCTransform with regressions
      seuratObject <- SCTransform(seuratObject, 
                                  method = "glmGamPoi",
                                  vars.to.regress = regressionVars,  # Regresses cell cycle/mito in SCT space
                                  verbose = TRUE)
      
      seuratObject <- RunPCA(seuratObject, npcs = numPCs, reduction.name = "pca", verbose = FALSE)

      DefaultAssay(seuratObject) <- "SCT"
    
      sampleList[[i]] <- seuratObject
    }
    
    # Integration prep
    features <- SelectIntegrationFeatures(sampleList, nfeatures = 3000, verbose = TRUE)
    sampleList <- PrepSCTIntegration(sampleList, anchor.features = features, verbose = TRUE)
    
    # Use RPCA reduction to avoid CCA fallbacks (for large data/logs)
    anchors <- FindIntegrationAnchors(sampleList, 
                                      normalization.method = "SCT",
                                      assay = NULL,
                                      reduction = "rpca",
                                      anchor.features = features,
                                      dims = dimsUsed,
                                      k.anchor = 5,
                                      k.filter = 200,
                                      k.score = 30,
                                      verbose = TRUE)
    
    seuratObject <- IntegrateData(anchors, normalization.method = "SCT", verbose = TRUE)
    VariableFeatures(seuratObject, assay = "SCT") <- features
    
    # For non-integrated SCT (run PCA on "SCT" assay)
    DefaultAssay(seuratObject) <- "SCT"
    seuratObject <- runDimReduction(seuratObject = seuratObject,
                                    assay = "SCT",
                                    reduction.output = "sct")
    
    # Analysis on integrated (IntegrateData creates "integrated" assay, so run PCA first)
    DefaultAssay(seuratObject) <- "integrated"
    seuratObject <- runDimReduction(seuratObject = seuratObject,
                                    assay = "integrated",
                                    reduction.output = "integrated")
    
    
  } else {
    stop("Invalid method. Use 'RPCA', or 'SCT'")
  }
  
  print("Integration complete. Seurat object ready for downstream analysis.")
  return(seuratObject)
}

