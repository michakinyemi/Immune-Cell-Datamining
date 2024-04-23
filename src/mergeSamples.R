
mergeSamples <- function(sampleList) {
  sampleCount = length(sampleList)
  
  #idList = list()
  
  for (i in 1:sampleCount) {
    sampleList[[i]] <- filterSample(sampleList[[i]])
    
    #idList <- append(idList, paste("s", i, sep=""))
  }
  
  mergedData <-
    merge(x = sampleList[[1]], y = sampleList[2:sampleCount], project = dataID)
  mergedData <- JoinLayers(mergedData)
  return(mergedData)
}
