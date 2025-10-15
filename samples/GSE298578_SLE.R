
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

# samples_sle <- loadDataset(SAMPLE_DATA)

# samples_sle <- integrateSamples(samples_sle, method="SCT")

# saveRDS(samples_sle, "data/samples_sle.rds")

samples_sle <- readRDS("data/samples_sle.rds")

# samples_sle <- annotateSinglerLabels(samples_sle)

# saveRDS(samples_sle, "data/samples_sle.rds")

findAllMarkers(samples_human, assay="SCT", group.by="clusters.integrated", fileName="sle")
