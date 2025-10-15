





visualizeVolcanoPlot <- function(markers) {
  library(ggplot2)

  ggplot(markers, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = ifelse(p_val_adj < 5e-250, "Significant", "Not Significant")), 
             alpha = 0.8, size = 2) +
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "black")) +
  labs(title = "Volcano Plot",
       x = "Log2(Fold Change)",
       y = "-log10(Adjusted P-value)",
       color = "Significance") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "top"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue")
}








#' --== Visualize UMAP Clusters ==--
#' 
#' 
visualizeClusters <- function(seuratObject, type="UMAP", resultsDir=NULL) {
  
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
#' Displays a heatmap comparing expression levels across cell types and conditions
#' 
visualizePrimaryHeatmap <- function(seuratObject, cellTypeCol=NULL, condCol="condition") {

  if (is.null(cellTypeCol)) {

  }


  # Set up annotations
  unique_combinations <- colnames(mat)

  anno_row <- data.frame(
    Pathway = factor(rep(c("Glycolysis", "T Cell Markers"), c(4, 2)))
  )

  anno_colors <- list(
    State = c("white", "firebrick"),
    CellType = c("white", "firebrick"),
    Pathway = c("Glycolysis" = "#7570B3", "T Cell Markers" = "#E7298A")
  )

  split_groups <- strsplit(unique_combinations, "_")

  anno_col <- data.frame(
    State = sapply(split_groups, "[[", 2),
    CellType = sapply(split_groups, "[[", 1)
  )

  rownames(anno_col) <- colnames(mat)



  
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

  figAvgExp <- pheatmap(
    mat,
    #annotation_col = anno_col,
    
    col = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    
    # breaks = seq(-1,1, length.out=(100 + 1)),
    
    cluster_cols = FALSE,
    cluster_rows = TRUE,
    angle_col = "45",
    gaps_col = seq(1, by = 1, length.out = 1)
  )

  ggsave("figures/Average Expression (via Scaled Clustering).png", figAvgExp, width=8, height=6, dpi=300)
  

  DotPlot(treg_singler, features = gene_list, cols = c("blue", "red"),
          dot.scale = 8, group.by = "orig.ident") +
          RotatedAxis()

  ggsave("figures/Dot Plot (via SingleR).png", figDotPlot, width=8, height=6, dpi=300)

  DotPlot(bm_samples, features=gene_list)
  
  return(fig)
}


#' --== Visualize Cell Type Proportions ==--
#' Displays the per sample 
#' 
visualizeCellTypeProportions <- function(seuratObject, metaCol="seurat_clusters")  {
  
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





#' --== Visualize Primary Heatmap ==--
#' Displays a heatmap comparing expression levels across cell types and conditions
#' 
visualizePrimaryHeatmap <- function(seuratObject, targets, cellTypeCol = NULL, conditionCol = "condition") {

  
  meta <- seuratObject@meta.data
  if (!conditionCol %in% colnames(meta)) stop("conditionCol not found in metadata.")
  if (!is.null(cellTypeCol) && !cellTypeCol %in% colnames(meta)) stop("cellTypeCol not found in metadata.")
  
  # Build group column
  if (!is.null(cellTypeCol)) {
    meta$group <- paste(meta[[cellTypeCol]], meta[[conditionCol]], sep = "_")
  } else {
    meta$group <- meta[[conditionCol]]
  }
  
  # Expression subset
  expr <- GetAssayData(seuratObject, slot = "data")[targets, , drop = FALSE]
  
  # Average per group
  avgExpr <- sapply(split(colnames(expr), meta$group), function(cells) {
    if (length(cells) == 1) expr[, cells, drop = FALSE]
    else rowMeans(expr[, cells, drop = FALSE])
  })
  avgExpr <- as.matrix(avgExpr)
  
  # Derive annotation from colnames
  if (!is.null(cellTypeCol)) {
    annot <- do.call(rbind, strsplit(colnames(avgExpr), "_"))
    colnames(annot) <- c("CellType", "Condition")
  } else {
    annot <- data.frame(Condition = colnames(avgExpr))
    rownames(annot) <- colnames(avgExpr)
  }
  annot <- as.data.frame(annot)
  rownames(annot) <- colnames(avgExpr)
  
  # Colors
  condLevels <- unique(annot$Condition)
  condCols <- setNames(brewer.pal(min(length(condLevels), 8), "Set1")[seq_along(condLevels)], condLevels)
  
  annotCols <- list(Condition = condCols)
  
  if (!is.null(cellTypeCol)) {
    cellLevels <- unique(annot$CellType)
    cellCols <- setNames(brewer.pal(min(length(cellLevels), 8), "Set3")[seq_along(cellLevels)], cellLevels)
    annotCols$CellType <- cellCols
  }
  
  # Plot
  pheatmap(
    avgExpr,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    scale = "row",
    show_colnames = FALSE,
    annotation_col = annot,
    annotation_colors = annotCols,
    cellwidth = 15,
    cellheight = 15,
    fontsize = 10
  )
}

