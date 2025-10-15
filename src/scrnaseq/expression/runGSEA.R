
#' --== Run Gene Set Enrichment Analysis (GSEA) ==--
#' 
runGSEA <- function(markers, database = "KEGG", species = "human") {

  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  library(ReactomePA)
  library(msigdbr)
  library(fgsea)

  # Prepare ranked gene list
  if (!("gene" %in% colnames(markers))) {
    markers$gene <- rownames(markers)
  }

  gene_list <- markers$avg_log2FC
  names(gene_list) <- markers$gene
  gene_list <- sort(gene_list, decreasing = TRUE)

  # Species reference
  reference <- if (species == "human") org.Hs.eg.db else org.Mm.eg.db

  gene_ids <- bitr(names(gene_list), fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = reference)

  gene_list <- gene_list[gene_ids$SYMBOL]
  names(gene_list) <- gene_ids$ENTREZID

  # Database selection
  if (database == "GO") {
    gsea_results <- gseGO(
      geneList = gene_list,
      OrgDb = reference,
      ont = "ALL",
      keyType = "ENTREZID",
      minGSSize = 10,
      maxGSSize = 500,
      pvalueCutoff = 0.05
    )
    return(gsea_results)

  } else if (database == "KEGG") {
    gsea_results <- gseKEGG(
      geneList = gene_list,
      organism = ifelse(species == "human", "hsa", "mmu"),
      minGSSize = 10,
      maxGSSize = 500,
      pvalueCutoff = 0.05
    )
    return(gsea_results)

  } else if (database == "REACTOME") {
    gsea_results <- gsePathway(
      geneList = gene_list,
      organism = ifelse(species == "human", "human", "mouse"),
      pvalueCutoff = 0.05,
      minGSSize = 10,
      maxGSSize = 500
    )
    return(gsea_results)

  } else if (database == "METABOLIC") {
    msig <- msigdbr(species = species, collection = "C2")

    # Keep only metabolism-related sets
    metab_sets <- msig[grep("METABOLISM", msig$gs_name), ]

    # Convert to ENTREZ IDs
    metab_ids <- bitr(metab_sets$gene_symbol, fromType = "SYMBOL",
                      toType = "ENTREZID", OrgDb = reference)
    metab_sets <- merge(metab_sets, metab_ids, by.x = "gene_symbol", by.y = "SYMBOL")
    metab_sets <- split(metab_sets$ENTREZID, metab_sets$gs_name)

    fgsea_res <- fgsea(
      pathways = metab_sets,
      stats = gene_list,
      minSize = 5,
      maxSize = 2000
    )
    return(list(fgsea_res, gene_list, metab_sets))
  }

}
