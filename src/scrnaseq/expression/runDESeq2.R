

runDESeq2 <- function(
  seuratObject,
  assays = "RNA",
  conditionCol = "condition",
  cellTypeCol = "cell_type",
  sampleCol   = "orig.ident",
  contrast    = c("condition", "Healthy", "Disease"),
  targetCellType = NULL,
  includeCellType = TRUE  # TRUE → ~ cell_type + condition, FALSE → ~ condition
) {
  # Optionally restrict to one cell type
  if (!is.null(targetCellType)) {
    seuratObject <- subset(
      seuratObject,
      subset = !!as.name(cellTypeCol) == targetCellType
    )
    # when a subset is specified, force design to condition only
    includeCellType <- FALSE
  }

  # Choose grouping variables for pseudobulk
  groupVars <- if (includeCellType) {
    c(conditionCol, sampleCol, cellTypeCol)
  } else {
    c(conditionCol, sampleCol)
  }

  # Aggregate pseudobulk expression
  agg <- AggregateExpression(
    seuratObject,
    assays = assays,
    return.seurat = TRUE,
    group.by = groupVars
  )

  counts <- GetAssayData(agg, assay = assays, slot = "counts")
  meta <- data.frame(cell = colnames(counts))

  # Rebuild metadata columns from combined group names
  parts <- strsplit(meta$cell, "_")
  meta[[conditionCol]] <- sapply(parts, `[`, 1)
  meta[[sampleCol]]    <- sapply(parts, `[`, 2)
  if (includeCellType) {
    meta[[cellTypeCol]] <- sapply(parts, `[`, 3)
  } else if (!is.null(targetCellType)) {
    meta[[cellTypeCol]] <- targetCellType
  }
  rownames(meta) <- meta$cell

  # Construct design formula
  designFormula <- if (includeCellType) {
    as.formula(paste("~", cellTypeCol, "+", conditionCol))
  } else {
    as.formula(paste("~", conditionCol))
  }

  # Build and run DESeq2 model
  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData   = meta,
    design    = designFormula
  )

  dds <- DESeq(dds)
  res <- results(dds, contrast = contrast)

  list(dds = dds, results = res)
}
