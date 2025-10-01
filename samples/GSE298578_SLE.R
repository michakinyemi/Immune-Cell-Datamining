
# --== Initialize Pipeline ==--
cat(getwd())

source("init.R")

SAMPLE_DATA <- list(
  "datasetName" = "GSE298578",
  "sampleDir" = "public/GSE298578/",
  "species" = "human", # Options: ("human", "mouse")
  "dataFormat" = "10X" # Options: ("10X", "counts", "rdata", "parsebio")
)

# --== Run Pipeline ==--

samples_sle <- loadDataset(SAMPLE_DATA)

samples_sle <- integrateSamples(samples_sle)

samples_sle <- runDimReduction(samples_sle, reduction.target="orig.pca")

visualize_batch_effect(samples_sle)

saveRDS(samples_sle, "data/samples_sle.rds")

library(presto)

samples_sle <- readRDS("data/samples_sle.rds")

samples_sle <- annotate_singler_labels(samples_sle)

saveRDS(samples_sle, "data/samples_sle.rds")

# samples_sle@misc$mainReduction <- "integrated.cca"

# find_all_markers(samples_sle)




