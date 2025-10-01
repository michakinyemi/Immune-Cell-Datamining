
# --== Initialize Pipeline ==--

source("init.R")

SAMPLE_DATA <- list(
  "datasetName" = "GSE175604",
  "sampleDir" = "/nguyen/GE8070/Deep/output/combined_human/",
  "species" = "human", # Options: ("human", "mouse")
  "dataFormat" = "parsebio" # Options: ("10X", "counts", "rdata", "parsebio")
)

# --== Run Pipeline ==--

samples_human <- load_dataset(SAMPLE_DATA)

samples_human <- perform_integration(samples_human)

samples_human <- run_dim_reduction(samples_human, reduction.target="orig.pca")

visualize_batch_effect(samples_human)

