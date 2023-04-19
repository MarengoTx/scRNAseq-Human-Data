## Table of Contents

1. [Main Description](#main-description)
2. [File Descriptions](#file-descriptions)
3. [Linked Files](#linked-files)
4. [Installation and Instructions](#installation-and-instructions)

&nbsp;
&nbsp;
&nbsp;
&nbsp;
&nbsp;
&nbsp;



### **Main Description**
------------------------

This is the GitHub repository for the single cell RNA sequencing data analysis for the human manuscript titled.

The following essential libraries are required for script execution:

```
Seurat
scReportoire
ggplot2
dplyr
ggridges
ggrepel
ComplexHeatmap

```
&nbsp;
&nbsp;
&nbsp;




### **File Descriptions**
---------------------------

The **code** can be downloaded and opened in **RStudios**. </br>


&nbsp;
&ensp;
&ensp;

### **Linked Files**
---------------------

This repository contains code for the analysis of single cell RNA-seq dataset. The dataset contains raw FASTQ files, as well as, the aligned files that were deposited in GEO. The **Rdata** or `.Rds` file was deposited in Zenodo. Provided below are descriptions of the linked datasets:

1. Gene Expression Omnibus (GEO) ID: GSE229626 
   - **Title**: Gene expression profile at single cell level of human T cells stimulated via antibodies against the T Cell Receptor (TCR)
   - **Description**: This submission contains the `matrix.mtx`, `barcodes.tsv`, and `genes.tsv` files for each replicate and condition, corresponding to the aligned files for single cell sequencing data. 
   - **Submission type**: Private. In order to gain access to the repository, you must use a [reviewer token](https://www.ncbi.nlm.nih.gov/geo/info/reviewer.html).

&ensp;

2. Sequence read archive (SRA) repository
    - **Title**: Gene expression profile at single cell level of human T cells stimulated via antibodies against the T Cell Receptor (TCR)
   - **Description**:  This submission contains the **raw sequencing** or `.fastq.gz` files, which are tab delimited text files. 
   - **Submission type**: Private. In order to gain access to the repository, you must use a [reviewer token](https://www.ncbi.nlm.nih.gov/geo/info/reviewer.html).
   
&ensp;


&nbsp;
&ensp;
&nbsp;
&ensp;




### **Installation and Instructions**
--------------------------------------
The code included in this submission requires several essential packages, as listed above. Please follow these instructions for installation:

&nbsp;

> Ensure you have R version 4.1.2 or higher for compatibility. 

> Although it is not essential, you can use R-Studios (Version 2022.12.0+353 (2022.12.0+353)) for accessing and executing the code. 

The following code can be used to set working directory in R:

> setwd(directory)

&nbsp;

**Steps**:
1. Download the **Human_code_April2023.R** and **Install_Packages.R** R scripts, and the processed data from **GSE229626**.
2. Open [R-Studios](https://www.rstudio.com/tags/rstudio-ide/) or a similar integrated development environment (IDE) for R. 
3. Set your working directory to where the following files are located:
   - `Human_code_April2023.R`
   - `Install_Packages.R`
4. Open the file titled `Install_Packages.R` and execute it in R IDE. This script will attempt to install all the necessary pacakges, and its dependencies.
&nbsp;
5. Open the `Human_code_April2023.R` R script and execute commands as necessary. 




