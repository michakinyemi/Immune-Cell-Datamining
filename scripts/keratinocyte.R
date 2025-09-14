

sampleList <- load_dataset(dataDir, format=dataFormat,
                           exclude=c("LesionB_Epi", "Normal_Epi"))

Dermis_samples <- merge_samples(sampleList)

Dermis_samples <- run_dim_reduction(Dermis_samples)

Dermis_samples <- perform_integration(Dermis_samples)

Dermis_samples <- regress_cell_cycle(Dermis_samples)

Dermis_samples <- FindNeighbors(Dermis_samples, reduction="regressed.cca",
                                dims=1:METRICS$numPCs)

Dermis_samples <- FindClusters(Dermis_samples, resolution=METRICS$clusterRes,
                               cluster.name="regressed.clusters")

Dermis_samples <- RunUMAP(Dermis_samples,
                          dims = 1:METRICS$numPCs,
                          reduction = "regressed.cca",
                          reduction.name = "regressed.umap")

Dermis_samples <- RunTSNE(Dermis_samples,
                          dims = 1:METRICS$numPCs,
                          check_duplicates = FALSE,
                          reduction = "regressed.cca",
                          reduction.name = "regressed.tsne")

Dermis_samples[["RNA"]] <- JoinLayers(Dermis_samples[["RNA"]])

SaveSeuratRds(Dermis_samples, file="Dermis_samples.rds" )

de_markers_dermis <- FindAllMarkers(Dermis_samples, only.pos=TRUE, min.pct=0.10,
                                    logfc.threshold=0.25, group.by="seurat_clusters")

write.csv(de_markers_dermis, "de_markers_dermis.csv")




# sampleList <- load_dataset(dataDir, format=dataFormat,
#                            exclude=c("LesionB_Dermis", "Normal_Dermis"))
# 
# Epi_samples <- merge_samples(sampleList)
# 
# Epi_samples <- run_dim_reduction(Epi_samples)
# 
# Epi_samples <- perform_integration(Epi_samples)
# 
# Epi_samples <- regress_cell_cycle(Epi_samples)
# 
# Epi_samples <- FindNeighbors(Epi_samples, reduction="regressed.cca",
#                               dims=1:METRICS$numPCs)
# 
# Epi_samples <- FindClusters(Epi_samples, resolution=METRICS$clusterRes,
#                              cluster.name="regressed.clusters")
# 
# Epi_samples <- RunUMAP(Epi_samples,
#                         dims = 1:METRICS$numPCs,
#                         reduction = "regressed.cca",
#                         reduction.name = "regressed.umap")
# 
# Epi_samples <- RunTSNE(Epi_samples,
#                         dims = 1:METRICS$numPCs,
#                         check_duplicates = FALSE,
#                         reduction = "regressed.cca",
#                         reduction.name = "regressed.tsne")
# 
# Epi_samples[["RNA"]] <- JoinLayers(Epi_samples[["RNA"]])
# 
# SaveSeuratRds(Epi_samples, file="Epi_samples.rds" )
# 
# de_markers_epi <- FindAllMarkers(Epi_samples, only.pos=TRUE, min.pct=0.10,
#                              logfc.threshold=0.25, group.by="seurat_clusters")
# 
# write.csv(de_markers_epi, "de_markers_epi.csv")
