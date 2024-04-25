
# --== filterSample() ==--
# Input/Usage

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
  
  
  # TODO: Is it preferred to normalize/scale data following the merge?
  # then again, we might not have enough memory anyways
  seuratObject <- NormalizeData(seuratObject)
  seuratObject <- ScaleData(seuratObject)
  
  return(seuratObject)
}
