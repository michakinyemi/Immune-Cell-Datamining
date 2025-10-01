
# --== Initialize Pipeline ==--

source("init.R")

SAMPLE_DATA <- list(
  "datasetName" = "GUT_DSS_FK1",
  "sampleDir" = "/nguyen/GE8070/Deep/output/combined_mouse/",
  "species" = "mouse", # Options: ("human", "mouse")
  "dataFormat" = "parsebio" # Options: ("10X", "counts", "rdata", "parsebio")
)

include_samples <- c("14_GUT_DSS_FK1FL_2", "15_GUT_DSS_FK1FL_3",
                     "16_GUT_DSS_FK1-FOXP3_4", "17_GUT_DSS_FK1-FOXP3_5")

# --== Run Pipeline ==--

samples_mouse <- load_dataset(SAMPLE_DATA, includeList = include_samples)

samples_mouse <- perform_integration(samples_mouse)

samples_mouse <- run_dim_reduction(samples_mouse, reduction.target="orig.pca")

visualize_batch_effect(samples_mouse)
