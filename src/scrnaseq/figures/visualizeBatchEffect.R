#' --== Visualize Batch Effect ==--
#' 
#' 
visualizeBatchEffect <- function(seuratObject, saveFolder="results/") {
  
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
    
    ggsave(paste0(saveFolder, "Batch Effect Analysis.png"), fig, create.dir=TRUE,
           width=15, height=9, dpi=300)
    
  } else {
    fig <- (fig1a + fig1b) / (fig2a + fig2b)
    
    ggsave(paste0(saveFolder, "Batch Effect Analysis.png"), fig, create.dir=TRUE,
           width=15, height=6, dpi=300)
  }
  
  
  return(fig)
}
