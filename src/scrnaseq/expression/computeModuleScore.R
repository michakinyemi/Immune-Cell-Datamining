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


