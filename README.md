# Single Cell Transcriptomics Tutorial 

**01.¿What is RNA-seq?¿What is scRNA-seq?¿spatial transcriptomics?** 🤔

**RNA-seq** are RNA sequence molecules from specific cells, tissues, or species, this allows differential expression of genes between different samples, identifies transcript variants, or captures all RNA in a given species. **Single-cell** scRNAseq differs from “traditional” bulk RNAseq. focuses on studying individual cells rather than averaging data across a population of cells. This approach reveals variations in gene expression and cellular states that can be masked when only bulk tissue samples are examined.scRNAseq is ideal for identifying cell types, understanding the development process, exploring disease mechanisms, or mapping cell lineages. On the other hand, **spatial transcriptomics** is a technique that maps gene expression to specific locations within tissue sections, preserving spatial information about the tissue architecture.

**02. How Does Single-Cell RNA-Seq Work?** 👩🏽‍🔬

1. Isolate single cells (techniques like microfluidics, fluorescence-activated cell sorting (FACS), or manual pipetting).
2. Cell Lysis and RNA Extraction
3. Reverse Transcription and cDNA Library Construction
4. Amplification and Sequencing
5. Data Analysis 

**03. Retrieve data** 🗂️ 	

We are going to navigate to https://www.10xgenomics.com/datasets. Filter products following your interests. 

**04. Clustering and Biomarkers** 🧬
Run the scRNAseq_tutorial.R script for doing the following steps:

1. Load a count matrix
2. QC and filtering
3. Normalization
4. Identify highly variable genes
5. Scale data
6. PCA
7. Clustering
8. Nonlinear-dimensionality reduction
9. Biomarkers

