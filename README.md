# Immune Cell Datamining

------------------------------------------------------------------------

This pipeline currently only supports scRNA-Seq datasets in the 10X genomics format.

All pre-processing steps are fully automated, with hands-on analysis steps being included in the `DataVisualization.Rmd` file.

**Additional Notes:**

-   Supports automatic conversion of ENSEMBL IDs to gene symbols.

-   Filtering metrics can be configured in the `DataVisualization.Rmd` file.

## Instructions

------------------------------------------------------------------------

1.  

2.  Source all R script files and load packages/config variables through the `DataVisualization.Rmd` file.

3.  Place 10X Genomics formatted scRNA-Seq data files into the `datasets/` directory.

4.  Read and store data into individual seurat object variables

    -   `samples <- generateSampleList(dataID)`

5.  Perform batch correction through `integrateData()` function or simply merge the samples normally through `mergeData()`.

    -   *Not necessary if only one sample is present in the dataset.*

    -   TODO: Create metrics for determining if batch correction is necessary

6.  Use `runDimReduction()` function to perform dimension reduction analysis on the processed dataset.

7.  Use code blocks in `DataVisualization.Rmd` file to generate figures and visualize clustering/markers.

    -   TODO: Better support image file generation for results

## Function Documentation

------------------------------------------------------------------------

(TODO) Usage instructions for all custom functions created for the pipeline:

| Function               | Usage |
|------------------------|-------|
| `generateSampleList()` |       |
| `filterSample()`       |       |
| `runDimReduction()`    |       |
| `annotateCellTypes()`  |       |
| `mergeData()`          |       |
| `integrateData()`      |       |

# Pipeline Planning

------------------------------------------------------------------------

## Implementation/Design

------------------------------------------------------------------------

-   Upon generation of a SeuratObject, the majority of sample-specific information will be imprinted into the `@misc` slot

## Planned Features

------------------------------------------------------------------------

**Analysis Techniques:**

-   Weighted Gene Co-Expression Network Analysis (WGCNA)

-   Clonality Trees

-   Trajectory Analysis/Pseudotime

-   Copy-Number Variations (CNVs)

**Pipeline Improvements:**

-   Bring expression level violin plots more in-line with past lab papers
-   Change symbol conversion function to update dataset files rather than being used during pipeline.
-   Support integration of separate datasets, while retaining sample characteristics/identifiers
-   Automatically classify as mouse vs human model (with manually override) and update gene references if necessary

### ATAC-Seq

------------------------------------------------------------------------

`filter_gtf.py` - Generates a subset GTF file containing only the annotations for the genes of interest. This prevents non-target genes from being included in the tracks.

## Additional Information

------------------------------------------------------------------------

### Quick Notes

It is unlikely that the information each database stores for entries will be the same

-   We can use marker genes (those that are highly associated with specific cell types) to differentiate

-   Papers involving sequencing data analysis often require you to include the steps you had taken along the way (filtering process, etc.)

    -   Seurat does a lot of the heavy lifting by storing data transformations
    -   We could still make this a part of the script's job by saving the most relevant parameters & results (Ex: \# of clusters) and generating an associated output file

-   Global Assignment (within functions): `variable <<- data`

### Future Plans

-   Sub Cluster Discovery Pipeline

    -   Support "zooming in" on clusters we had generated to further investigate

**Questions:**

-   Most documentation/papers suggest using high [% mitochondrial gene expression as a filtering metric]{.underline}. How does this impact the study of genes such as **TFAM** (Mitochondrial Transcription Factor A)?

### Resources

-   [(Seurat) Guided Clustering Tutorial](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html)

-   [(Elixir Excelerate) Single RNA-seq data analysis with R (Elixir Excelerate)](https://nbisweden.github.io/excelerate-scRNAseq/session-qc/Quality_control.html)

-   [(Harvard) Introduction to Bioinformatics and Computational Biology](https://liulab-dfci.github.io/bioinfo-combio/)
