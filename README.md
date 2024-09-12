# scRNA-seq Data Analysis Tutorial

In this tutorial you will use the Kallisto | bustools workflow to perform pseudo-alignment of scRNA-seq reads to a reference transcriptome and generate count matrices. Then you will analyze the count data in R. 

------------------------------------------------------------------------
## Overview 

- **Part1**: Get a count matrix from fastq files. 
  - Use the `prefetch` and `fasterq dump` functions from the  [SRA-toolkit](https://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump) to download fastq files from the SRA.   
  - Use `kb ref` from the [kallisto | bustools](https://www.kallistobus.tools) workflow to download the pre-made mouse reference index.    
  - Use `kb-count`  from the [kallisto | bustools](https://www.kallistobus.tools) to get cell by gene count data.  
  - Use `wget` to download processed data from GEO.   
  - Part1 [Tutorial](https://msubioinformaticscore.github.io/scRNAseq_training/html/countsFromFastq.html)
- **Part2**: Analyze the count data in R. 
  - Make a `SingleCellExperiment` object from count data derived from.   `Kaslisto-Bustools` or processed counts downloaded from GEO.    
  - Detect empty droplets with `DropletUtils`.    
  - Detect barcodes that correspond to 'doublets'.   
  - Identify and remove low quality cells from your data.   
  - Normalize and transform the raw counts for downstream analysis    
  - Select highly variable genes in the data.   
  - Perform dimesnionality reductions using PCA, tSNE and UMAP.    
  - Cluster the cells based on their gene expression profiles.    
  - Identify cluster enriched genes.   
  - Assign cell labels from gene sets.    
  - Calculate mean gene expression per clusters or cell type.   
  - Visualize your results.  
  - Part2 [Tutorial](https://msubioinformaticscore.github.io/scRNAseq_training/html/bioconductor_scRNAseq_analysis.html)

------------------------------------------------------------------------

## Set up

### Download workshop data

The data used in this workshop are all publicly available and the download links are included throughout the document. For efficiency's sake, you can dowload the data and code from Google Drive, [here](https://drive.google.com/drive/folders/1gfrdzceslLBtHPVnfB8sSIbV70rdvbkh?usp=drive_link). Log into [ondemand.hpcc.msu.edu](https://ondemand.hpcc.msu.edu) to upload the `scRNAseq_training` directory into your hpcc space.

### Install Anaconda

If you haven't yet follow these instructions to [install Anaconda](https://docs.icer.msu.edu/Using_conda/).

### Install `kb-python` if you haven't

From these [instructions](https://bustools.github.io/download).
```
module purge
module load Conda/3
pip install kb-python
```

### Install R packages

Open the R environment on HPCC
```
module purge
module load R-bundle-CRAN/2023.12-foss-2023a
R --vanilla
```

In the R environment:
```
#install cran packages
install.packages(c("tidyverse", "Matrix", "patchwork",
                   "pheatmap", "RColorBrewer", "readxl"))


#install bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("SingleCellExperiment", "scater",
                     "scran", "DropletUtils", "bluster",
                     "scDblFinder", "AUCell", "PCAtools"))
```

------------------------------------------------------------------------

## Part1: Get a count matrix from fastq files

 - [html](https://msubioinformaticscore.github.io/scRNAseq_training/html/countsFromFastq.html)
 - [Rmarkdown](https://github.com/MSUBioinformaticsCore/scRNAseq_training/blob/main/src/countsFromFastq.Rmd)

------------------------------------------------------------------------

## Part2: Analyze the count data in R

 - [html](https://msubioinformaticscore.github.io/scRNAseq_training/html/bioconductor_scRNAseq_analysis.html)
 - [Rmarkdown](https://github.com/MSUBioinformaticsCore/scRNAseq_training/blob/main/src/bioconductor_scRNAseq_analysis.Rmd)

------------------------------------------------------------------------

## Credits

Part2 of this tutorial is based off of the e-book [*Orchestrating Single-Cell Analysis with Bioconductor*](https://bioconductor.org/books/release/OSCA/index.html), a comprehensive resource designed to guide users through the process of analyzing single-cell RNA sequencing (scRNA-seq) data using the Bioconductor ecosystem in R. 
