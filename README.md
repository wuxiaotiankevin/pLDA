# pLDA

This is the implementation of Penalized Latent Dirichlet Allocation (pLDA) as described in Xiaotian Wu, Hao Wu, Zhijin Wu. “Penalized Latent Dirichlet Allocation Model in Single Cell RNA Sequencing”. In preperation.

# Quick Start

## Download and Install Package
```
git clone https://github.com/wuxiaotiankevin/pLDA
R CMD build pLDA
R CMD INSTALL pLDA_0.1.1.tar.gz
```

## Analyze Single Cell RNA Sequencing Data with pLDA
```
library(pLDA)
fit = plda(x=CELL_BY_GENE_EXPRESSION_COUNT_MATRIX, k=NUMBER_OF_TOPICS, lambda=3)
```
