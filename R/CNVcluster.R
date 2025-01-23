#' Perform CNV Clustering with Seurat Object
#'
#' The `CNVcluster` function performs hierarchical clustering on a genomic score matrix extracted from a Seurat object.
#' It provides options for plotting a dendrogram, an elbow plot for optimal cluster determination,
#' and cluster visualization on the dendrogram. The resulting cluster assignments are stored in the Seurat object.
#'
#' @param seuratObj A Seurat object containing a "genomicScores" assay with a matrix of genomic scores for clustering.
#' @param k Optional. Number of clusters to cut the dendrogram into. If `NULL`, the optimal number of clusters is determined automatically using the elbow method.
#' @param h Optional. The height at which to cut the dendrogram for clustering. If both `k` and `h` are provided, `k` takes precedence.
#' @param plotDendrogram Logical. Whether to plot the dendrogram. Defaults to `FALSE`.
#' @param plotClustersOnDendrogram Logical. Whether to highlight clusters on the dendrogram. Defaults to `FALSE`.
#' @param plotElbowPlot Logical. Whether to plot the elbow plot used for determining the optimal number of clusters. Defaults to `FALSE`.
#'
#' @details
#' The function computes a Manhattan distance matrix and performs hierarchical clustering using the Ward.D2 method.
#' If `k` is not provided, the elbow method is applied to determine the optimal number of clusters based on within-cluster sum of squares (WSS).
#'
#' The clusters are assigned to the Seurat object under the metadata column `cnv_clusters`.
#'
#' @return A Seurat object with an additional metadata column, `cnv_clusters`, containing cluster assignments.
#'
#' @import proxy
#' @import utils
#' @import graphics
#'
#' @export


CNVcluster <- function(seuratObj,
                         k = NULL,
                         h = NULL,
                         plotDendrogram = F,
                         plotClustersOnDendrogram = F,
                         plotElbowPlot = F) {

  if (is.null(k)){kDetection = "automatic"}
  if (!is.null(k)){kDetection = "manual"}

  genomicMatrix <- t(as.matrix(Seurat::GetAssay(seuratObj, assay = "genomicScores")$data))

  dist_cos <- proxy::dist(genomicMatrix, method = "Manhattan")

  dist_matrix <- proxy::as.dist(dist_cos)
  hc <- hclust(dist_matrix, method = "ward.D2")


  if (plotDendrogram){plot(hc, main = "Dendrogram",
                           xlab = "Observations", ylab = "Distance", labels = F)}

  dist_matrix_full <- as.matrix(dist_matrix)


  find_elbow <- function(x, y, sensitivity = 0.05) {
    x_norm <- (x - min(x)) / (max(x) - min(x))
    y_norm <- (y - min(y)) / (max(y) - min(y))

    slopes <- diff(y_norm) / diff(x_norm)
    slope_changes <- diff(slopes)

    significant_changes <- which(abs(slope_changes) > sensitivity * max(abs(slope_changes)))

    if (length(significant_changes) == 0) {
      return(x[ceiling(length(x)/2)])
    }

    first_quarter <- ceiling(length(x) / 4)
    elbow_point <- significant_changes[significant_changes > first_quarter][1]

    if (is.na(elbow_point)) {
      elbow_point <- utils::tail(significant_changes, 1)
    }

    return(x[elbow_point + 1])
  }



  k_values <- 1:20
  wss <- sapply(k_values, function(k) {
    clusters <- cutree(hc, k = k)
    sum(sapply(unique(clusters), function(cl) {
      sum(dist_matrix_full[which(clusters == cl), which(clusters == cl)]^2) / (2 * sum(clusters == cl))
    }))
  })

  if(kDetection == "automatic") {k <- find_elbow(k_values, wss)}

  if(plotElbowPlot){
    plot(k_values, wss, type = "b", pch = 19, frame = FALSE,
         xlab = "Number of Clusters (k)",
         ylab = "WSS",
         main = "Elbow Plot ")

    graphics::abline(v = k, col = "red", lty = 2)
    if(kDetection == "automatic") {graphics::text(k, max(wss), paste("optimal k =", k), pos = 4, col = "red")}
  }


  clusters <- cutree(hc, k = k, h = h)

  if (plotClustersOnDendrogram){rect.hclust(hc, k = k, h = h, border = "red")}

  seuratObj$cnv_clusters = as.factor(clusters)

  return(seuratObj)
}





# CNV Clustering within a Seurat Object
#
# Performs clustering based on Copy Number Variations (CNVs) calculated by `CNVcalling`.
#
# @param seuratObj A Seurat object containing single-cell genomic data.
# @param nPCs The number of principal components to use for the clustering. Default is 10.
# @param resolution Resolution parameter for the clustering algorithm. Affects the granularity of the resulting clusters. Default is 0.8.
# @param k.param Number of nearest neighbors to use in the construction of the SNN graph. Default is 20.
#
# @return A Seurat object with CNV-based clusters added
# @import Seurat SeuratObject
# @export
#
# CNVcluster <- function(seuratObj, nPCs = 10, resolution = 0.8, k.param = 20) {
#   dft_a <- DefaultAssay(seuratObj)
#   exp_clusters <- seuratObj[["seurat_clusters"]]
#   DefaultAssay(seuratObj) <- "genomicScores"
#
#   seuratObj <- ScaleData(seuratObj, assay = "genomicScores", features = Features(seuratObj, assay = "genomicScores"))
#   seuratObj <- RunPCA(seuratObj, assay = "genomicScores", reduction.name = "cnv_pca", reduction.key = "cnvPC_",
#                       features = Features(seuratObj, assay = "genomicScores"))
#   seuratObj <- FindNeighbors(seuratObj, reduction = "cnv_pca", dims = 1:nPCs)
#   seuratObj <- FindClusters(seuratObj, resolution = resolution, k.param = k.param)
#
#   DefaultAssay(seuratObj) <- dft_a
#   seuratObj[["cnv_clusters"]] <- seuratObj[["seurat_clusters"]]
#   seuratObj[["seurat_clusters"]] <- exp_clusters
#
#   return(seuratObj)
# }
