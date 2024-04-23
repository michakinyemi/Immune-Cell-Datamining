generateDataVis <- function(genesOfInterest){
  
  VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)
  
  out1 <- x + y
  out2 <- t + v
  
  
  if (is.null(genesOfInterest)) {
    cat("target genes not provided")
  }
  
}


# fun <- function(){
#   x <- readline("What is the value of x? ")
#   y <- readline("What is the value of y? ")
#   t <- readline("What are the T values? ")
#   v <- readline("What are the V values? ")
# 
#   x <- as.numeric(unlist(strsplit(x, ",")))
#   y <- as.numeric(unlist(strsplit(y, ",")))
#   t <- as.numeric(unlist(strsplit(t, ",")))
#   v <- as.numeric(unlist(strsplit(v, ",")))
# 
#   out1 <- x + y
#   out2 <- t + v
# 
#   return(list(out1, out2))
# 
# }



