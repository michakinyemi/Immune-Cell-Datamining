#' --== Initialize Environment ==--
#' 
#' If the environment is already initialized, it will instead:
#' 1. Re-parse global config and marker files and update env variables
#' 2. Re-source all global R scripts (i.e. "/src")
#' 
.init_env <- function(remote=FALSE) {
  
  CONFIG$timestamp <<- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  
  CONFIG$resultsDir <<- paste(CONFIG$resultsHome, CONFIG$timestamp, sep="/")
  
  if (remote) {
    .libPaths("/home/mi542876/.conda/envs/immune_r/lib/R/library")
    CONFIG$dataHome <<- CONFIG$dataHomeRemote
    
    if (!require("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    
    BiocManager::install("celldex")
    
  } else {
    CONFIG$dataHome <<- CONFIG$dataHomeLocal 
  }
  
  # Load External Packages
  library(Seurat)
  library(dplyr)
  library(data.table)
  library(gridExtra)
  library(ggplot2)
  library(SingleR)
  library(celldex)
  # library(DESeq2)
  require(pheatmap)
  require(RColorBrewer)
  library(presto)


  # Source Custom Analysis Scripts
  for (f in list.files(path=CONFIG$srcDir, pattern="\\.R$", recursive=TRUE)) {
    print(f)
    source(paste(CONFIG$srcDir, f, sep="/"))
  }
  rm(f)

  CONFIG$s.genes <<- Seurat::cc.genes$s.genes
  CONFIG$g2m.genes <<- Seurat::cc.genes$g2m.genes

  REFERENCES <<- list(
    "celldex.human" = celldex::HumanPrimaryCellAtlasData(),
    "celldex.mouse" = celldex::MouseRNAseqData()
  )

  # Retrieve gene symbols from ensembl IDs
  # keyrefMm <- mapIds(org.Mm.eg.db, keys=counts@Dimnames[[1]], column="SYMBOL", keytype="ENSEMBL")
  # keyrefHs <- mapIds(org.Hs.eg.db, keys=counts@Dimnames[[1]], column="SYMBOL", keytype="ENSEMBL")
  
  # Import Markers
  MARKERS <<- yaml::read_yaml("markers.yml")
}


#' --== Validate Environment ==--
#' Ensure that all required packages are installed (and usable?)
#' 
#' 
.validate_env <- function() {
  if (!"dplyr" %in% loadedNamespaces()) library(dplyr)
  
  if (exists("METRICS") == FALSE) return(-1)
  
}



#' --== Prepare Output ==--
#' 
#' 
 .prepare_output_dir <- function(seuratObject) {
  
  if (is.null(seuratObject@misc$resultsPath)) {
    resultsPath <- paste0(CONFIG$resultsHome, seuratObject@project.name)
    
    if (file.exists(resultsPath)) {
      
    }
  }
  
}


#' --== Validate Save Path ==--
#' 
#' 
.validate_save_path <- function(path) {
  
  
  
  
}
