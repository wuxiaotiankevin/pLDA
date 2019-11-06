# Introduction
Single cell RNA sequencing (scRNA-seq) is a recently developed technology that allows quantification of RNA transcripts at individual cell level, providing cellular level resolution of gene expression variation. We develop the penalized Latent Dirichlet Allocation (pLDA) model to extract robust and interpretable biological information from single cell mRNA sequencing (scRNA-seq) data. The method is adapted from the generative probabilistic model LDA originated in natural language processing. pLDA models the scRNA-seq data by considering genes as words, cells as documents, and latent biological functions as topics, with a penalty that increases the robustness of the estimation. The topics identified by pLDA are interpretable with biological functions.

For details, see Xiaotian Wu, Hao Wu, Zhijin Wu. “Penalized Latent Dirichlet Allocation Model in Single Cell RNA Sequencing”. In preperation.

The package is maintained by Xiaotian (Kevin) Wu. Contact xiaotian_wu at brown dot edu for questions.

# Getting Started

## Install
```
library(devtools)
install_github("wuxiaotiankevin/pLDA", build_vignettes=TRUE)
```

## Quick Start
```
library(pLDA)
```
For a simulation example, see
```
browseVignettes('pLDA')
```
