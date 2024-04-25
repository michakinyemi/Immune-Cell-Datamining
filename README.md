# Immune Cell Datamining

------------------------------------------------------------------------

This pipeline currently only supports scRNA-Seq datasets in the 10X genomics format.

All pre-processing steps are fully automated, with hands-on analysis steps being included in the `DataVisualization.Rmd` file.

**Additional Notes:**

-   Supports automatic conversion of ENSEMBL IDs to gene symbols.

-   Filtering metrics can be configured in the `DataVisualization.Rmd` file.

## Instructions

------------------------------------------------------------------------

1.  Source all R script files and load packages/config variables through the `DataVisualization.Rmd` file.
2.  Place 10X Genomics formatted scRNA-Seq data files into the `datasets/` directory.
3.  Read and store data into individual seurat object variables
    -   `samples <- generateSampleList(dataID)`
4.  Perform batch correction through `integrateData()` function or simply merge the samples normally through `mergeData()`.
    -   *Not necessary if only one sample is present in the dataset.*

    -   TODO: Create metrics for determining if batch correction is necessary
5.  Use `runDimReduction()` function to perform dimension reduction analysis on the processed dataset.
6.  Use code blocks in `DataVisualization.Rmd` file to generate figures and visualize clustering/markers.
    -   TODO: Better support image file generation for results

## Documentation

------------------------------------------------------------------------

(TODO) Usage instructions for all custom functions created for the pipeline:

`generateSampleList()`

`filterSample()`

`runDimReduction()`

`annotateCellTypes()`

`mergeData()`

`integrateData()`

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
