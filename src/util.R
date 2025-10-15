
#' --== Find Gene ==--
#' 
#' 
convert_to_symbol <- function(geneList, species=NULL) {
  #features <- data@Dimnames[[1]]
  
  if (species == "mouse") database = org.Mm.eg.db
  else if (species == "human") database = org.Hs.eg.db
  else {
    cat("Error: No species argument provided and no CONFIG detected.")
    return(-1)
  }
  
  print(typeof(database))
  
  geneTable <- mapIds(
    org.Mm.eg.db,
    keys = geneList, #rownames(data)
    column = "SYMBOL",
    keytype = "ENSEMBL"
  )
  
  inCommon <- geneList %in% names(geneTable)
  geneList[inCommon] <-
    unlist(geneTable[geneList[inCommon]])
  
  geneList <- toupper(geneList)
  
  return(geneList)
}


#' --== Find Gene ==--
#' 
#' 
find_gene <- function(seuratObject, gene) {
  grep(gene, SeuratObject::Features(seuratObject), value = TRUE)
  # Features(seuratObject)[Features(seuratObject) %like% gene]
}

#' Title: My custom function
#'
#' Description: Explains what the function does.
#' --== Subset By Metadata ==--
#'
#' @param x Numeric input.
#' @param y Numeric input.
#' @return Sum of x and y.
#' @examples
#' myfun(1, 2)
#' @export
subset_by_metadata <- function(seuratObject, metaCol, values) {
  if (!metaCol %in% colnames(seuratObject@meta.data)) {
    stop(paste("Column", metaCol, "not found in metadata"))
  }
  subset(seuratObject, cells = rownames(seuratObject@meta.data[
    seuratObject@meta.data[[metaCol]] %in% values, 
  ]))
}

#' --== Print All Markers ==--
#' 
#' 
printTopMarkers <- function(markers, n=20) {
  top_markers <- markers %>% group_by(cluster) %>% top_n(n=n, wt=avg_log2FC)
  
  top_markers %>%
    group_by(cluster) %>%
    summarise(genes = paste(gene, collapse = ", ")) %>%
    mutate(output = paste0("Cluster ", cluster, " Top Markers:\n[", genes, "]")) %>%
     pull(output) %>%
    cat(sep = "\n\n")
  
  return(top_markers)
}