
#' --== Load Dataset ==--
#' Parse the provided sample directory and generate
#' a Seurat Object for each identified sample.
#' 
loadDataset <- function(SAMPLE_DATA, excludeList=vector(), includeList=vector(),
                         convertFeatures=FALSE, sampleMetaCol=NULL) {
  
  # Validate environment and argument setup
  if (exists("CONFIG") == FALSE) {
    cat("Error: Could not find CONFIG within environment.")
    return(-1)
  }
  
  exclusionMode = FALSE
  inclusionMode = FALSE
  if (length(excludeList) > 0 & length(includeList) > 0) {
    cat("Error: Please only specify one include or exclude list, not both.")
    return(-1)
  } else {
    if (length(includeList) != 0) {
      inclusionMode = TRUE
    
    # Enter exclusion mode by default (to ensure exclusion of meta folders)
    } else {
      exclusionMode = TRUE
      excludeList <- append(excludeList, CONFIG$defaultExclude)
    }
  }

  final_dir <- file.path(CONFIG$dataHome, SAMPLE_DATA$sampleDir)
  
  cat("Beginning to Load Dataset:", SAMPLE_DATA$datasetName, "\n", sep="")
  cat("Parsing Directory:\n", final_dir, "\n", sep="")
  
  sampleDirList <- list.dirs(path=final_dir,
                             recursive=FALSE)
  sampleList <- list() # The final list of accepted samples
  tempList <- sampleDirList
  
  # Single Sample Setup
  if (length(sampleDirList)==0) {
    cat("No subdirectories found. Setting sample count to 1.\n")
    sampleDirList = append(sampleDirList, final_dir)
    sampleCount = 1
    
  # Multi Sample Setup
  } else {
    sampleCount = 0
    for (folder in tempList) {
      include = TRUE
      
      # Inclusion Filter
      if (inclusionMode) {
        if (basename(folder) %in% includeList) {
          sampleCount = sampleCount + 1
        } else {
          sampleDirList <- sampleDirList[sampleDirList != folder]
          include = FALSE
        }
      }
      
      # Exclusion Filter
      else if (exclusionMode) {
        for (word in excludeList) {
          if (grepl(word, folder)) {
            sampleDirList <- sampleDirList[sampleDirList != folder]
            include = FALSE
          }
        }
        if (include) { sampleCount = sampleCount + 1 }
      }
      
    }
  }
  
  cat(sprintf("Total Samples Identified: %d\n", sampleCount))
  
  format <- tolower(SAMPLE_DATA$dataFormat)
  
  i = 0
  for (sampleDir in sampleDirList) {
    i = i + 1
    sampleName = basename(sampleDir)
    cat(sprintf("\nLoading Sample [%d/%d]: %s\n", i, sampleCount, sampleName))
    
    if (format == "10x") {
      
      data <- Seurat::Read10X(sampleDir)
      
      if (convertFeatures) {
        data <- convert_to_symbol(data)
        duplicate_rows <- duplicated(rownames(data))
        data <- data[!duplicate_rows, ]
      }
      
      sampleObj <- CreateSeuratObject(counts=data, project=sampleName)
      sampleObj$orig.sample <- sampleName
      
      cat("Filtering Sample:", sampleName, "\n")
      sampleObj <- filterSample(sampleObj)
      
    }
    
    else if (format == "counts") {
      counts <- read.csv(file.path(sampleDir, "gene_counts.csv.gz"), row.names="X")
      counts <- as(as.matrix(counts), "dgCMatrix")
      
      sampleObj <- CreateSeuratObject(counts, project = sampleName)
      metadata <- read.csv(file.path(sampleDir, "barcodes_samples_info.csv.gz"))
      rownames(metadata) <- metadata$barcode
      sampleObj <- AddMetaData(sampleObj, metadata)
      sampleObj$orig.sample <- sampleName
      
      cat("Filtering Sample:", sampleName, "\n")
      sampleObj <- filterSample(sampleObj)
    }

    else if (format == "parsebio") {
      sampleDir <- paste0(sampleDir, "/DGE_filtered")
      
      mat <- ReadParseBio(sampleDir)
      
      rownames(mat)[rownames(mat) == ""] <- "unknown"
      
      cell_meta <- read.csv(paste0(sampleDir, "/cell_metadata.csv"), row.names = 1)
      
      sampleObj <- CreateSeuratObject(mat, meta.data = cell_meta)
      
      cat("Filtering Sample:", sampleName, "\n")
      sampleObj <- filterSample(sampleObj)
    }
    
    else if (format == "rdata") {
      load(file.path(sampleDir, paste0(sampleName, ".Rdata")))
      sampleObj <- data
      sampleObj$orig.sample <- sampleName
    }
    
    # Handle splitting samples by metadata column
    if (!is.null(sampleMetaCol)) {
      cat("Splitting Sample by Metadata Column:", sampleMetaCol, "\n")
      splitSamples <- SplitObject(sampleObj, split.by = sampleMetaCol)
      for (subName in names(splitSamples)) {
        sub <- splitSamples[[subName]]
        sub$orig.sample <- subName
        sub@project.name <- subName
        sampleList <- append(sampleList, list(sub))
      }
    }
    else {
      sampleList <- append(sampleList, list(sampleObj))
    }
  }
  
  if (length(sampleList) == 1) { 
    return(sampleList[1])
  } else {
    return(sampleList)
  }
  
}


#' --== Filter Samples ==--
#' 
#' 
filterSample <- function(seuratObject) {
  seuratObject$percent.mt <- PercentageFeatureSet(seuratObject,
                                                  pattern = "^MT-")
  
  origCount = length(Cells(seuratObject))
  
  seuratObject <- subset(
    seuratObject,
    subset = nFeature_RNA > CONFIG$minFeatRNA &
             percent.mt < CONFIG$mtExprLimit
  )
  
  newCount = length(Cells(seuratObject))
  filtCount = origCount - newCount
  
  cat("- Original Cell Count:", origCount, "\n")
  cat("- N Cells Removed:", filtCount, "\n")
  cat("- New Cell Count:", newCount, "\n")
  
  return(seuratObject)
}



 



