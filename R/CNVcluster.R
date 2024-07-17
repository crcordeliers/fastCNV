#' CNV Clustering within a Seurat Object
#'
#' Performs clustering based on Copy Number Variations (CNVs) calculated by `CNVcalling`.
#'
#' @param seuratObj A Seurat object containing single-cell genomic data.
#' @param nPCs The number of principal components to use for the clustering. Default is 10.
#' @param resolution Resolution parameter for the clustering algorithm. Affects the granularity of the resulting clusters. Default is 0.8.
#' @param k.param Number of nearest neighbors to use in the construction of the SNN graph. Default is 20.
#'
#' @return A Seurat object with CNV-based clusters added
#' @import Seurat SeuratObject
#' @export

CNVcluster <- function(seuratObj, nPCs = 10, resolution = 0.8, k.param = 20) {
  dft_a <- DefaultAssay(seuratObj)
  exp_clusters <- seuratObj[["seurat_clusters"]]
  DefaultAssay(seuratObj) <- "genomicScores"

  seuratObj <- ScaleData(seuratObj, assay = "genomicScores", features = Features(seuratObj, assay = "genomicScores"))
  seuratObj <- RunPCA(seuratObj, assay = "genomicScores", reduction.name = "cnv_pca", reduction.key = "cnvPC_",
                      features = Features(seuratObj, assay = "genomicScores"))
  seuratObj <- FindNeighbors(seuratObj, reduction = "cnv_pca", dims = 1:nPCs)
  seuratObj <- FindClusters(seuratObj, resolution = resolution, k.param = k.param)

  DefaultAssay(seuratObj) <- dft_a
  seuratObj[["cnv_clusters"]] <- seuratObj[["seurat_clusters"]]
  seuratObj[["seurat_clusters"]] <- exp_clusters

  return(seuratObj)
}
