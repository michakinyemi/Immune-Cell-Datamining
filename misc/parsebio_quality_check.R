# Parse Biosciences Testing

# sample <- "12_TUMOR_HFD_ML"

# filteredDir <- sprintf("/media/michael/Nguyen-Lab/sequencing-data/nguyen/GE8070/output_combined/combined_mouse/%s/DGE_unfiltered", sample)

library(ggplot2)
library(dplyr)

METRICS = list(
  "minFeatRNA" = 800,
  "mtExprLimit" = 20, # The highest % of MT gene reads allowed
  "numPCs" = 30, # Amount of principal components to generate
  "clusterResolution" = 1
)

results <- data.frame(sampleID = character(),
                      mean_mt_percent = numeric(),
                      post_mt_filter_count = integer(),
                      post_n200_filter_count = integer(),
                      post_n800_filter_count = integer(),
                      stringsAsFactors = FALSE)

for (f in list.dirs(path="/media/michael/Nguyen-Lab/sequencing-data/nguyen/GE8070/Deep/output/combined_mouse",
                    recursive = FALSE)) {
  
  sampleID <- basename(f)

  if (length(grep(sampleID, list("all-sample", "process"))) == 0) {
    print(sampleID)
    
    data_unfiltered <- Seurat::ReadParseBio(paste(f, "/DGE_unfiltered", sep=""))
    pb_unfiltered <- Seurat::CreateSeuratObject(counts = data_unfiltered)
    pb_unfiltered[["percent.mt"]] <- Seurat::PercentageFeatureSet(pb_unfiltered, pattern = "^mt-")
    pb_unfiltered$filter_status <- "Unfiltered"
    
    data_filtered <- Seurat::ReadParseBio(paste(f, "/DGE_filtered", sep=""))
    pb_filtered <- Seurat::CreateSeuratObject(counts = data_filtered)
    pb_filtered[["percent.mt"]] <- Seurat::PercentageFeatureSet(pb_filtered, pattern = "^mt-")
    pb_filtered$filter_status <- "Filtered"
    
    pb_combined <- merge(pb_unfiltered, pb_filtered)
    
    figure <- Seurat::VlnPlot(pb_combined, features = "percent.mt", group.by = "filter_status") &
      geom_hline(yintercept=20, linetype="dashed", color = "red") &
      ylim(0, 100)
    
    ggsave(sprintf("/media/michael/Nguyen-Lab/sequencing-data/nguyen/GE8070/Deep/quality/[MT Filter] %s.png", sampleID), figure, width=6, height=8, dpi=300)
    
    mean_mt <- mean(pb_unfiltered$percent.mt, na.rm=TRUE)
    mean_nfeature <- mean(pb_unfiltered$nFeature_RNA)

    mt_filter_count <- nrow(pb_filtered@meta.data %>% filter(percent.mt < METRICS$mtExprLimit))
    n200_filter_count <- nrow(pb_filtered@meta.data %>% filter(percent.mt < METRICS$mtExprLimit & nFeature_RNA > 200))
    n800_filter_count <- nrow(pb_filtered@meta.data %>% filter(percent.mt < METRICS$mtExprLimit & nFeature_RNA > 800))
    
    print(sprintf("Post MT Filter Count: %s", mt_filter_count))
    print(sprintf("Post >200 nFeature Filter Count: %s", n200_filter_count))
    print(sprintf("Post >800 nFeature Filter Count: %s", n800_filter_count))
    
    results <- rbind(results, data.frame(sampleID = sampleID,
                                         mean_mt_percent = mean_mt,
                                         post_mt_filter_count = mt_filter_count,
                                         post_n200_filter_count = n200_filter_count,
                                         post_n800_filter_count = n800_filter_count,
                                         stringsAsFactors = FALSE))
  }
}
rm(f,sampleID)


