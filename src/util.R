

#' --== Initialize Environment ==--
#' 
#' If the environment is already initialized, it will instead:
#' 1. Re-parse global config and marker files and update env variables
#' 2. Re-source all global R scripts (i.e. "/src")
#' 
.init_env <- function(remote=FALSE) {
  
  CONFIG$timestamp <<- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  
  CONFIG$resultsDir <<- paste(CONFIG$resultsHome, CONFIG$timestamp, sep="/")
  
  if (remote) {
    .libPaths("/home/mi542876/.conda/envs/immune_r/lib/R/library")
    CONFIG$dataHome <<- CONFIG$dataHomeRemote
    
    if (!require("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    
    BiocManager::install("celldex")
    
  } else {
    CONFIG$dataHome <<- CONFIG$dataHomeLocal 
  }
  
  # Load External Packages
  library(Seurat)
  library(dplyr)
  library(data.table)
  library(gridExtra)
  library(ggplot2)
  library(SingleR)
  library(celldex)
  library(DESeq2)


  # Source Custom Analysis Scripts
  for (f in list.files(path=CONFIG$srcDir, pattern="*.R")) {
    print(f)
    source(paste(CONFIG$srcDir, f, sep="/"))
  }
  rm(f)

  CONFIG$s.genes <<- Seurat::cc.genes$s.genes
  CONFIG$g2m.genes <<- Seurat::cc.genes$g2m.genes

  # Retrieve gene symbols from ensembl IDs
  # keyrefMm <- mapIds(org.Mm.eg.db, keys=counts@Dimnames[[1]], column="SYMBOL", keytype="ENSEMBL")
  # keyrefHs <- mapIds(org.Hs.eg.db, keys=counts@Dimnames[[1]], column="SYMBOL", keytype="ENSEMBL")
  
  # Import Markers
  MARKERS <<- yaml::read_yaml("markers.yml")
}


#' --== Validate Environment ==--
#' Ensure that all required packages are installed (and usable?)
#' 
#' 
.validate_env <- function() {
  if (!"dplyr" %in% loadedNamespaces()) library(dplyr)
  
  if (exists("METRICS") == FALSE) return(-1)
  
}



#' --== Prepare Output ==--
#' 
#' 
 .prepare_output_dir <- function(seuratObject) {
  
  if (is.null(seuratObject@misc$resultsPath)) {
    resultsPath <- paste0(CONFIG$resultsHome, seuratObject@project.name)
    
    if (file.exists(resultsPath)) {
      
    }
  }
  
}


#' --== Validate Save Path ==--
#' 
#' 
.validate_save_path <- function(path) {
  
  
  
  
}



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

#' --== Annotate Cells Using SingleR ==--
#' 
#' 
annotate_singler_labels <- function(seuratObject,
                                    species="human") {

  if (species == "human") {
    ref <- ref_human
  } else if (species == "mouse") {
    ref <- ref_mouse
  }

  sce <- as.SingleCellExperiment(seuratObject)

  # Main Cell Type Labeling
  singler_results_main <- SingleR(
    test = sce,
    ref = ref,
    labels = ref$label.main
  )

  seuratObject$SingleR_main <- singler_results_main$labels

  # Fine Cell Type Labeling
  singler_results_fine <- SingleR(
    test = sce,
    ref = ref,
    labels = ref$label.fine
  )
  
  seuratObject$SingleR_fine <- singler_results_fine$labels
  
  return(seuratObject)
  
}

#' --== Manually Annotate Clusters ==--
#' 
#' 
annotateClusters <- function(seuratObject, markerGenes=NULL, useSingleR=FALSE) {
  # Initialize a new column for cell type annotations
  seuratObject$cell_type <- NA
  
  # Loop through each cluster and annotate cells
  for (cell_type in names(markerGenes)) {

    # Identify cells expressing marker genes for the current cell type
    marker_genes_expr <- markerGenes[[cell_type]]
    cell_indices <-
      WhichCells(seuratObject, expression = marker_genes_expr)
    
    # Assign the cell type label to the corresponding cells
    seuratObject$cell_type[cell_indices] <- cell_type
  }
  
  # Return the annotated Seurat object
  Idents(seuratObject) <- "cell_type"

  return(seuratObject)
}




#' --== Find Gene ==--
#' 
#' 
find_gene <- function(seuratObject, gene) {
  grep(gene, SeuratObject::Features(seuratObject), value = TRUE)
  # Features(seuratObject)[Features(seuratObject) %like% gene]
}


#' --== Subset By Metadata ==--
#' 
#' 
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
print_top_markers <- function(markers, n=20) {
  top_markers <- markers %>% group_by(cluster) %>% top_n(n=n, wt=avg_log2FC)
  
  top_markers %>%
    group_by(cluster) %>%
    summarise(genes = paste(gene, collapse = ", ")) %>%
    mutate(output = paste0("Cluster ", cluster, " Top Markers:\n[", genes, "]")) %>%
     pull(output) %>%
    cat(sep = "\n\n")
  
  return(top_markers)
}

#' --== Find All Markers ==--
#' 
#' 
find_all_markers <- function(seuratObject, n=20) {
  
  clusters_used <- paste0(seuratObject@misc$mainReduction, ".clusters")
  
  de_markers <- FindAllMarkers(seuratObject, only.pos=TRUE,
                               min.pct=0.10, logfc.threshold=0.25,
                               group.by=clusters_used)
  
  write.csv(de_markers, "results/all_markers_sle.csv", row.names = FALSE)
  
  topMarkers <- de_markers %>% group_by(cluster) %>% top_n(n=n, wt=avg_log2FC)
  
  write.csv(topMarkers, "results/top_markers_sle.csv", row.names = FALSE)
}


