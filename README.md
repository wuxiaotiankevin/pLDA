# Introduction
We develop the penalized Latent Dirichlet Allocation (pLDA) model to extract robust and interpretable biological information from single cell mRNA sequencing (scRNA-seq) data. The method is adapted from the generative probabilistic model LDA originated in natural language processing. pLDA models the scRNA-seq data by considering genes as words, cells as documents, and latent biological functions as topics, with a penalty that increases the robustness of the estimation. We apply pLDA to scRNA-seq datasets from both Drop-seq and SMARTer v1 technologies, and demonstrate improved performances in cell type classification. The topics identified by pLDA are interpretable with biological functions.

For details, see Xiaotian Wu, Hao Wu, Zhijin Wu. “Penalized Latent Dirichlet Allocation Model in Single Cell RNA Sequencing”. In preperation.

The package is maintained by Xiaotian (Kevin) Wu. Contact xiaotian_wu at brown dot edu for questions.

# Getting Started

## Install
```
git clone https://github.com/wuxiaotiankevin/pLDA
R CMD build pLDA
R CMD INSTALL pLDA_0.1.1.tar.gz
```

## Quick Start
```
library(pLDA)
browseVignettes('pLDA')
fit = plda(x = CELL_BY_GENE_EXPRESSION_COUNT_MATRIX, k = NUMBER_OF_TOPICS, lambda = YOUR_PENALTY)
```
