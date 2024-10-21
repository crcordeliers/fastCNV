## Introduction

`fastCNV` is a package that helps you detect, plot and analyse the putative Copy Number Variations (CNVs) in single cell (scRNA-seq) data or Spatial Transcriptomics (ST) data. Built on `SeuratObject`, it is easily integrated into scRNA-seq or ST pipelines.

*WARNING:* Project is still under construction and function usage may change.

## Installation

To install fastCNV, run in R:

```         
remotes::install_github("must-bioinfo/fastCNV")
```

## Capabilities

`fastCNV` can plot a heatmap of inferred CNVs: ![](https://must-bioinfo.github.io/fastCNV/articles/fastCNV_files/figure-html/run1-1.png)

It also calculates a `cnv_fraction`, which can be plotted with `Seurat` standard plotting functions:

![](https://must-bioinfo.github.io/fastCNV/articles/fastCNV_files/figure-html/plot_cnv_umap-1.png)

`cnv_fractions` can also be visualized spatially:

![](https://must-bioinfo.github.io/fastCNV/articles/fastCNV_files/figure-html/plot_cnv_fraction2-1.png)

And `cnv_fractions` can be used to obtain clonal clusters (`cnv_clusters`):

![](https://must-bioinfo.github.io/fastCNV/articles/fastCNV_files/figure-html/plot_cnv_clusters-1.png)

## Usage

An extensive tutorial is available to get started [here](https://must-bioinfo.github.io/fastCNV/articles/fastCNV.html).
