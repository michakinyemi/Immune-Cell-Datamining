

#' --== Calculate Gene Module Score ==--
#' 
#' TODO: Get better metrics than the simple mean
#' 
computeModuleScore <- function(seuratObject, genes, score_name = "ModuleScore") {
  # Add module score
  seuratObject <- AddModuleScore(seuratObject, features = list(genes), name = score_name)
  
  # The score column gets a number suffix (e.g. ModuleScore1)
  score_col <- paste0(score_name, "1")
  
  # Compute average per cluster
  avg_scores <- seuratObject@meta.data %>%
    dplyr::group_by(cluster = Idents(seuratObject)) %>%
    dplyr::summarise(mean_score = mean(.data[[score_col]], na.rm = TRUE))
  
  return(list(seuratObject = seuratObject, cluster_scores = avg_scores))
}



runDESeq2 <- function(seuratObject, 
                                conditionCol = "condition", 
                                cellTypeCol = "cell_type",
                                sampleCol   = "orig.ident",
                                contrast    = c("condition", "Healthy", "SLE"),
                                targetCellType = NULL) {

  
  # raw counts
  counts <- GetAssayData(seuratObject, slot = "counts")
  meta   <- seuratObject@meta.data
  
  # restrict to one cell type if requested
  if (!is.null(targetCellType)) {
    keep_cells <- rownames(meta)[meta[[cellTypeCol]] == targetCellType]
    counts <- counts[, keep_cells, drop = FALSE]
    meta   <- meta[keep_cells, , drop = FALSE]
  }
  
  # make group ID (sample × cell type)
  meta$group <- paste(meta[[sampleCol]], meta[[cellTypeCol]], sep = "_")
  
  # aggregate counts to pseudobulk
  pseudobulk_counts <- do.call(cbind, lapply(unique(meta$group), function(g) {
    cells <- rownames(meta[meta$group == g, ])
    if (length(cells) == 1) counts[, cells, drop=FALSE]
    else rowSums(counts[, cells, drop=FALSE])
  }))
  colnames(pseudobulk_counts) <- unique(meta$group)
  storage.mode(pseudobulk_counts) <- "integer"

  
  # build pseudobulk metadata
  pb_meta <- data.frame(
    sample    = sub("_(.*)$", "", colnames(pseudobulk_counts)),
    cell_type = sub(".*_", "", colnames(pseudobulk_counts))
  )
  pb_meta[[conditionCol]] <- meta[[conditionCol]][match(pb_meta$sample, meta[[sampleCol]])]
  rownames(pb_meta) <- colnames(pseudobulk_counts)

  keep <- !is.na(pb_meta[[conditionCol]])
  pseudobulk_counts <- pseudobulk_counts[, keep]
  pb_meta <- pb_meta[keep, ]
  
  # construct DESeq2 object
  dds <- DESeqDataSetFromMatrix(
    countData = pseudobulk_counts,
    colData   = pb_meta,
    design    = ~ cell_type + get(conditionCol)
  )
  
  dds <- DESeq(dds)
  res <- results(dds, contrast = contrast)
  
  return(list(dds = dds, results = res))
}

