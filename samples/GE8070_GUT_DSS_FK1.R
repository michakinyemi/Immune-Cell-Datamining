
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

samples_mouse <- loadDataset(SAMPLE_DATA, includeList=include_samples)

samples_mouse <- integrateSamples(samples_mouse, method="SCT")

samples_mouse <- annotateSinglerLabels(samples_mouse, species="mouse")

saveRDS(samples_mouse, "data/samples_mouse.rds")

findAllMarkers(samples_mouse, assay="SCT", group.by="clusters.integrated", fileName="mouse")
