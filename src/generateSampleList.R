
# --== generateSampleList() ==--
# Input/Usage
#   dataID              Name of the folder containing your dataset
#   ignoreList          List of sample substrings to ignore
#                       ("meta" included by default)
#   convertFeatures     Enable if gene IDs must be converted to gene symbols

generateSampleList <- function(dataID,
                               ignoreList = vector(),
                               convertFeatures = FALSE) {
  
  dataDirList <-
    list.dirs(path = paste("datasets/", dataID, sep = ""),
              recursive = FALSE)
  
  sampleList <-  vector()
  ignoreList <- append(ignoreList, "meta")
  tempList <- dataDirList
  
  if (length(dataDirList)==0) single=TRUE
  else single=FALSE
  
  if (single) {
    dataDirList = append(dataDirList, paste("datasets/", dataID, sep = ""))
    sampleCount = 1
  } else {
    # Remove meta folder and ignored samples
    sampleCount = 0
    for (folder in tempList) {
      include = TRUE
      for (word in ignoreList) {
        if (grepl(word, folder)) {
          dataDirList <- dataDirList[dataDirList != folder]
          include = FALSE
        }
      }
      if (include)
        sampleCount = sampleCount + 1
    }
  }
  
  i = 0
  # Read 10X folder, convert to Seurat, and filter
  for (folder in dataDirList) {
    i = i + 1
    if (single)
      name = dataID
    else
      name = strsplit(folder, split = "/")[[1]][3]
    
    # TODO: By default take final suffix (i.e. "_UMM061") as sample names
    # and provide user with optional manual input/override aaa
    
    sprintf("\nLoading Sample [$d/$d]")
    
    cat("\nLoading Sample [", i, "/", sampleCount, "]: ", name, "\n", sep = "")
    
    data <- Read10X(folder)
    
    if (convertFeatures) {
      data <- convertToSymbol(data)
      duplicate_rows <- duplicated(rownames(data))
      data <- data[!duplicate_rows, ]
    }
    
    data <- CreateSeuratObject(counts = data)
    data$orig.sample = name
    
    cat("Filtering", name)
    data <- filterSample(data)
    
    # TODO: Handle interactive addition of samples
    # if (isInteractive) {}
    
    # Return single Seurat object for singular samples
    if (single)
      return(data)
    
    sampleList <- append(sampleList, data)
  }
  
  return(sampleList)
}


# TODO Features:
# - Auto/manually classify as mouse/human model



