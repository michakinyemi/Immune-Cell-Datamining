

samples <- c("LesionB_Dermis", "LesionB_Epi", "Normal_Dermis", "Normal_Epi")


for (dataID in samples) {
  dataDir <- sprintf("/media/michael/Nguyen-Lab/sequencing-data/public/GSE191335/%s", dataID)
  
  data <- Read10X(folder)
  data <- CreateSeuratObject(counts = data, project=dataID)
  
  save(data, file=sprintf("%s.Rdata", dataID))
}
