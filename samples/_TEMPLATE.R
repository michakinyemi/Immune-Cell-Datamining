
# --== Initialize Pipeline ==--

source("init.R")

SAMPLE_DATA <- list(
  "datasetName" = "GSE175604",
  "sampleDir" = "/nguyen/GE8070/Deep/output/combined_human/",
  "species" = "human", # Options: ("human", "mouse")
  "dataFormat" = "parsebio" # Options: ("10X", "counts", "rdata", "parsebio")
)


include_samples <- c("14_GUT_DSS_FK1FL_2", "15_GUT_DSS_FK1FL_3",
                     "16_GUT_DSS_FK1-FOXP3_4", "17_GUT_DSS_FK1-FOXP3_5")

# --== Run Pipeline ==--

samples <- load_dataset(SAMPLE_DATA)

samples <- perform_integration(samples)

samples <- run_dim_reduction(samples, reduction.target="orig.pca")

samples <- annotate_singler_labels(samples, species = "human")

visualize_batch_effect(samples)


