# Brown CBC's Bioconductor Workshop 2 Materials

[![Docs](https://img.shields.io/badge/docs-stable-blue.svg?style=flat-square)](https://compbiocore.github.io/bioconductor-workshop-2)
[![nbviewer](https://img.shields.io/badge/jupyter_notebooks-nbviewer-purple.svg?style=flat-square)](http://nbviewer.jupyter.org/github/compbiocore/bioconductor-workshop-2/tree/master/docs/src/notebooks)
[![License](https://img.shields.io/github/license/compbiocore/bioconductor-workshop-2.svg)](https://raw.githubusercontent.com/compbiocore/bioconductor-workshop-2/master/LICENSE)


## Goals

This workshop, 'Analyses of RNA-seq and ChiP-seq with R/Bioconductor', is primarily designed to introduce attendees DEseq and its analysis pipeline; it is an intermediate-level lecture that presupposes some knowledge of R like that found in our 'Fundamentals of R' workshop.  This workshop is **not** a direct sequel to 'R/Bioconductor Workshop for Genomic Data Analysis', and that material is not required to understand this content.

Here, users will learn how to perform exploratory data analysis on RNA-seq data, then use that data (if appropriate) to conduct a complete differential expression and gene-set analysis.  The workshop will discuss how to interpret the results of the exploratory analyses it presents and what they say about the data's quality, as well as offer advice on how to make the data more suitable for analysis.  It takes the user from a raw table of counts all the way to a set of differentially expressed genes and their pathways - the traditional endpoint for an RNA-seq analysis.

Finally, this workshop also includes a brief introduction to ChIP-seq analysis.  This topic exists in its own self-contained notebook (Notebook 4), for which the RNA-seq material is not a prerequisite.


## Workshop Topic Overview

This workshop, 'Analyses of RNA-seq and ChiP-seq with R/Bioconductor', is designed to complement (rather than succeed) the Bioconductor training provided in the first Bioconductor workshop; it is an intermediate-level lecture that presupposes some knowledge of R like that found in our 'Fundamentals of R' workshop.  The topics presented herein include the following:

* Standard analyses for RNA-seq Data
    * Approaches to processing raw data from RNA seq experiments
    * Annotating and importing expression counts
    * Basic exploratory data analysis
    * Assessing model fit and dispersion

* Gene-wise Testing in Depth, Including Multiple Testing
    * Differential expression analysis with the 'edgeR' package
    * Comparison of the 'edgeR' and 'DESeq2' packages for differential expression analysis
    * Hypothesis testing for count data
    * Approaches to multiple hypothesis testing result correction (e.g. false discovery rate calculation)
    * Visualizing differentially expressed genes/transcripts

* Gene Set Analysis
    * Over-representation analysis
    * Gene set enrichment analysis
    * Network-based enrichment analysis
    * Region-based enrichment analysis
    * Resources for these types of analyses (e.g. GO, KEGG, MSigDb)

* Analyses of ChIP-seq Data
    * Introduction to standard approaches used in ChIP Seq
    * Analyses for peak calling
    * Calculation of differential binding across samples
    * The annotation and plotting of peaks


### **[Notebook](http://nbviewer.jupyter.org/github/compbiocore/bioconductor-workshop-2/tree/master/docs/src/notebooks)**

To access a read-only version of the workshop notebooks (complete with full output), please click the above link; you will be redirected to Jupyter's NBViewer utility.  Simply click the name of the notebook and it will be rendered as HTML.

### Authors

Levi Waldron

### Contact

For assistance, contact cbc-help@brown.edu - this is our general help line, so please specify that your issue is with this site's contents

### Original Time and Location

> Date: 2/19/2018

> Location: Campus Center
