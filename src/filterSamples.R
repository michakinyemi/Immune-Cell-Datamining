# --== filterSample() ==--
# Input/Usage
#
filterSample <- function(seuratObject) {
  seuratObject$percent_MT <-
    PercentageFeatureSet(seuratObject, pattern = "^MT-")
  
  origCount = length(Cells(seuratObject))
  
  seuratObject <- subset(
    seuratObject,
    subset = nFeature_RNA > FILT_METRICS$minFeatRNA &
      nCount_RNA > FILT_METRICS$maxFeatRNA &
      percent_MT < FILT_METRICS$mtExprLimit
  )
  
  newCount = length(Cells(seuratObject))
  filtCount = origCount - newCount
  
  # Original Cell Count:
  # Filtered Cell Count:
  # (X Cells Removed)
  
  
  return(seuratObject)
}