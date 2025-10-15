#' --== Find All Markers ==--
#' 
#' 
findAllMarkers <- function(seuratObject, n=20, assay=NULL, group.by=NULL, fileName=NULL) {

  if (assay == "SCT") {
    seuratObject <- PrepSCTFindMarkers(seuratObject)
  }
  
  degMarkers <- FindAllMarkers(seuratObject, assay=assay,
                               min.pct=0.01, logfc.threshold=0.1,
                               group.by=group.by)
  
  topMarkers <- degMarkers %>% group_by(cluster) %>% top_n(n=n, wt=avg_log2FC)

  if (!is.null(fileName)) {
    write.csv(degMarkers, paste0(CONFIG$resultsHome, "/all_markers_", fileName, ".csv"),
              row.names = FALSE)
    
    write.csv(topMarkers, paste0(CONFIG$resultsHome, "/top_markers_", fileName, ".csv"),
              row.names = FALSE)
  }

  return(c(degMarkers, topMarkers))
}
