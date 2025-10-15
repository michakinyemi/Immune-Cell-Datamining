
# --== Initialize Pipeline ==--

source("init.R")

SAMPLE_DATA <- list(
  "datasetName" = "GSE175604",
  "sampleDir" = "/nguyen/GE8070/Deep/output/combined_human/",
  "species" = "human", # Options: ("human", "mouse")
  "dataFormat" = "parsebio" # Options: ("10X", "counts", "rdata", "parsebio")
)

# --== Run Pipeline ==--

samples_human <- loadDataset(SAMPLE_DATA)

samples_human <- integrateSamples(samples_human, method="SCT")

samples_human <- annotateSinglerLabels(samples_human, species="human")

saveRDS(samples_human, "data/samples_human.rds")

findAllMarkers(samples_human, assay="SCT", group.by="clusters.integrated", fileName="human")