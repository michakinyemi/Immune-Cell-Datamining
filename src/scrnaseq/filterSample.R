#' --== Filter Samples ==--
#' 
#' 
filterSample <- function(seuratObject) {
  seuratObject$percent.mt <- PercentageFeatureSet(seuratObject,
                                                  pattern = "^MT-")
  
  origCount = length(Cells(seuratObject))
  
  seuratObject <- subset(
    seuratObject,
    subset = nFeature_RNA > CONFIG$minFeatRNA &
             percent.mt < CONFIG$mtExprLimit
  )
  
  newCount = length(Cells(seuratObject))
  filtCount = origCount - newCount
  
  cat("- Original Cell Count:", origCount, "\n")
  cat("- N Cells Removed:", filtCount, "\n")
  cat("- New Cell Count:", newCount, "\n")
  
  return(seuratObject)
}
