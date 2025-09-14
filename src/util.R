

.initialize_env <- function(remote=FALSE) {
  
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
  library(jsonlite)
  library(Seurat)
  library(dplyr)
  library(data.table)
  library(gridExtra)
  library(ggplot2)
  library(SingleR)
  
  library(celldex)


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
  MARKERS <<- fromJSON("markers.json")
}

.validate_env <- function() {
  if (!"dplyr" %in% loadedNamespaces()) library(dplyr)
  
  if (exists("METRICS") == FALSE) return(-1)
  
}



#' --== Prepare Output ==--
#' 
#' 
  <- function(seuratObject) {
  
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

annotate_cell_types <- function(
    seuratObject,
    markerGenes=NULL) {
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

find_gene <- function(seuratObject, gene) {
  grep(gene, SeuratObject::Features(seuratObject), value = TRUE)
  # Features(seuratObject)[Features(seuratObject) %like% gene]
}

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
  
  resultsDir <- 
  
  clusters_used <- paste0(seuratObject@misc$main_reduc, ".clusters")
  
  de_markers_pos <- FindAllMarkers(seuratObject, only.pos=TRUE,
                               min.pct=0.10, logfc.threshold=0.25,
                               group.by=clusters_used)
  
  write.csv(topMarkers, "results/topMarkers.csv", row.names = FALSE)
  
  topMarkers <- de_markers_epi %>% group_by(cluster) %>% top_n(n=20, wt=avg_log2FC)
  
  write.csv(topMarkers, "results/topMarkers.csv", row.names = FALSE)
  
  
  
}


#' --== Generate Main Figures ==--
#' 
#' 
generate_main_figures <- function(seuratObject, saveFolder=NULL) {
  
  visualize_batch_effect(seuratObject, saveFolder=saveFolder)
  
  visualize_umap_clusters(seuratObject, saveFolder=saveFolder)
  
}

#' --== Visualize Batch Effect ==--
#' 
#' 
visualize_batch_effect <- function(seuratObject, saveFolder="results/") {
  
  fig1a <- DimPlot(seuratObject, reduction="orig.pca.umap", 
                   group.by="orig.ident") + ggtitle('1. Initial Dimension Reduction')
  fig1b <- DimPlot(seuratObject, reduction="orig.pca.umap",
                   group.by="orig.ident", split.by="orig.ident")
  
  fig2a <- DimPlot(seuratObject, reduction="integrated.cca.umap",
                   group.by="orig.ident") + ggtitle('2. Post-Batch Effect Correction')
  fig2b <- DimPlot(seuratObject, reduction="integrated.cca.umap",
                   group.by="orig.ident", split.by="orig.ident")
  
  if (!is.null(seuratObject@reductions$regressed.umap)) {
    
    fig3a <- DimPlot(seuratObject, reduction="regressed.umap",
                     group.by="orig.ident") + ggtitle('3. Post-Cell Cycle Regression')
    fig3b <- DimPlot(seuratObject, reduction="regressed.umap",
                     group.by="orig.ident", split.by="orig.ident")
    
    fig <- (fig1a + fig1b) / (fig2a + fig2b) / (fig3a + fig3b)
    
    ggsave(paste0(saveFolder, "Batch Effect Analysis.png"), fig, width=15, height=9, dpi=300)
    
  } else {
    
    fig <- (fig1a + fig1b) / (fig2a + fig2b)
    
    ggsave(paste0(saveFolder, "Batch Effect Analysis.png"), fig, width=15, height=6, dpi=300)
  }
  
  
  return(fig)
}


#' --== Visualize UMAP Clusters ==--
#' 
#' 
visualize_umap_clusters <- function(seuratObject, resultsDir=NULL) {
  
  if (is.null(resultsDir)) { resultsDir = CONFIG$resultsDir }
  
  # TODO:
  # - Better figure title/axis labels
  # - Include cluster resolution used in subtitle
  
  fig <- DimPlot(seuratObject, reduction="orig.pca.umap",
                 label = TRUE, group.by="orig.pca.clusters")
  ggsave(paste(resultsDir, "[DimPlot] Original UMAP.png", sep="/"), fig,
         width=6, height=4, dpi=300)
  
  fig <- DimPlot(seuratObject, reduction="integrated.cca.umap",
                 label = TRUE, group.by="integrated.cca.clusters")
  ggsave(paste(resultsDir, "[DimPlot] CCA UMAP.png", sep="/"), fig,
         width=6, height=4, dpi=300)
  
  if (!is.null(seuratObject@reductions$regressed.cca.umap)) {
    fig <- DimPlot(seuratObject, reduction="regressed.cca.umap", label = TRUE,
                   group.by="regressed.cca.clusters")
    ggsave(paste(resultsDir, "[DimPlot] Regressed UMAP.png", sep="/"), fig,
           width=6, height=4, dpi=300)
  
  }
  
  return(fig)
}

#' --== Visualize Primary Heatmap ==--
#' 
#' 
visualize_primary_heatmap <- function(seuratObject) {
  
  fig <- pheatmap(
    mat,
    #annotation_col = anno_col,
    
    col = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    
    # breaks = seq(-1,1, length.out=(100 + 1)),
    
    cluster_cols = FALSE,
    cluster_rows = TRUE,
    angle_col = "45",
    gaps_col = seq(1, by = 1, length.out = 1)
  )
  
  ggsave("figures/Average Expression (via Scaled Clustering).png", fig, width=8, height=6, dpi=300)
  
  return(fig)
}


#' --== Visualize Cell Type Proportions ==--
#' 
#' 
visualize_cell_type_proportions <- function(seuratObject, metaCol="seurat_clusters")  {
  
  cell_metadata <- seuratObject@meta.data[, c(metaCol, "orig.ident")]
  
  # Calculate Individual Proportions
  proportions <- cell_metadata %>%
    group_by(orig.ident, seurat_clusters) %>%
    summarise(count = n()) %>%
    mutate(proportion = count / sum(count)) %>%
    ungroup()
  
  # Calculate average proportions
  average_proportions <- proportions %>%
    group_by(seurat_clusters) %>%
    summarise(average_proportion = mean(proportion)) %>%
    mutate(orig.ident = "Average")
  
  # Normalize average proportions
  total_average_proportion <- sum(average_proportions$average_proportion)
  average_proportions <- average_proportions %>%
    mutate(proportion = average_proportion / total_average_proportion) %>%
    dplyr::select(-average_proportion)
  
  # Combine proportions and normalized average proportions
  combined_proportions <- bind_rows(proportions, average_proportions)
  
  # Add a dummy column for spacing
  combined_proportions <- combined_proportions %>%
    mutate(orig.ident = factor(orig.ident, levels = c(unique(proportions$orig.ident), "", "Average")))
  
  # Plot proportions with a wider average bar and spacing
  fig <- ggplot() +
    # Plot individual sample proportions
    geom_bar(data = combined_proportions,
             aes(x = orig.ident, y = proportion, fill = seurat_clusters),
             stat = "identity") +
    theme_minimal() +
    labs(title = "Cell Type Proportions",
         x = "Original Sample",
         y = "Proportion") +
    scale_y_continuous(labels = scales::percent) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  return(fig)
}





