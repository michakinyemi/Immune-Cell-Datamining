# Load External Packages
library(Seurat)
library(WGCNA)
library(yaml)
library(ggplot2)
library(ggThemeAssist)
library(dplyr)
library(data.table)
library(biomaRt)
library(jsonlite)
library(RColorBrewer)

library(org.Hs.eg.db)
library(org.Mm.eg.db)

# Source Custom Analysis Scripts
for (f in list.files(path="src/", pattern="*.R")) {
  print(f)
  source(paste("src/", f, sep=""))
}
rm(f)


# --== Script Config ==-- 

dataID = "transplant"

FILT_METRICS = list(
  "minFeatRNA" = 500,
  "maxFeatRNA" = 800,
  "mtExprLimit" = 10 # The highest % of MT gene reads allowed
) 

NUM_PCS = 50 # The # of principal components to generate

CELL_TYPES = c("metastatic", "primary")

# Retrieve gene symbols from ensembl IDs
#keyrefMm <- mapIds(org.Mm.eg.db, keys=counts@Dimnames[[1]], column="SYMBOL", keytype="ENSEMBL")
#keyrefHs <- mapIds(org.Hs.eg.db, keys=counts@Dimnames[[1]], column="SYMBOL", keytype="ENSEMBL")

CTM = fromJSON("markers.json")




samples <- parseDataset(dataID)
samples <- mergeSamples(samples)
samples <- runDimReduction(samples)
sample <- integrateSamples(samples)


