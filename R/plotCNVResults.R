#' A function to plot the CNV results into a heatmap
#'
#' @param seuratObj Seurat object containing the data used to get the genomicScores.
#' @param referenceVar The variable name of the annotations in the Seurat object metadata
#' @param clustersVar The variable name of the clusters in the Seurat object metadata
#' @param splitPlotOnVar The variable name on which to split the heatmap rows.
#' @param savePath Path to save the pdf heatmap. If `NULL`, plot won't be saved (default = `.`).
#' @param printPlot If the heatmap should be printed in the console.
#' @param referencePalette The color palette that should be used for `referenceVar`.
#' @param clusters_palette The color palette that should be used for `clustersVar`.
#' @param outputType "png" or "pdf".
#'
#' @import ComplexHeatmap
#' @import Seurat
#' @import circlize
#' @import qpdf
#' @import grDevices
#' @import ggplot2
#' @import patchwork
#' @import grid
#' @import paletteer
#' @import magick
#' @import scales
#'
#' @return This function builds a heatmap and saves it as a .pdf file in your working directory
#'
#' @export
#'

plotCNVResults <- function(seuratObj, referenceVar = NULL, clustersVar = "cnv_clusters",
                           splitPlotOnVar = clustersVar, savePath = ".",
                           printPlot = FALSE, referencePalette = as.character(paletteer::paletteer_d("pals::glasbey")),
                           clusters_palette = scales::hue_pal(),
                           outputType = "png"){
  if (outputType != "png" && outputType != "pdf"){
    message("Warning : outputType not valid, should be 'pdf' or 'png'. Setting outputType to 'png'")
    outputType = "png"
  }

  M <- t(as.matrix(Seurat::GetAssay(seuratObj, "genomicScores")["data"]))

  if (!is.null(referenceVar)) {
    annotation_df <- as.data.frame(seuratObj@meta.data[[referenceVar]])
    colnames(annotation_df) <- "Annotations"
    annot_colors <- setNames(referencePalette[1:length(unique(annotation_df$Annotations))], unique(annotation_df$Annotations))
    if (!is.null(splitPlotOnVar)) {
      split_df <- as.data.frame(Seurat::FetchData(seuratObj, vars = splitPlotOnVar))
      colnames(split_df) <- "Split"
    } else {
      split_df <- NULL
    }
  } else if (is.null(splitPlotOnVar)) {
    split_df <- NULL
  }
  if (!is.null(clustersVar)) {
    clusters_df <- as.data.frame(seuratObj@meta.data[[clustersVar]])
    colnames(clusters_df) <- "Clusters"
    clusters_colors <- setNames(clusters_palette(length(unique(clusters_df$Clusters))), unique(clusters_df$Clusters))
  }

  if (!is.null(referenceVar) && is.null(clustersVar)) {
    annotation_heatmap <- ComplexHeatmap::rowAnnotation(
      Annotations = annotation_df$Annotations,
      col = list(Annotations = annot_colors),
      annotation_legend_param = list(
        title = "Annotations",
        title_gp = grid::gpar(fontsize = 11),
        labels_gp = grid::gpar(fontsize = 8),
        legend_height = grid::unit(3, "cm"),
        legend_width = grid::unit(1.5, "cm"),
        grid_height = grid::unit(0.6, "cm"),
        grid_width = grid::unit(0.6, "cm")
      )
    )
  }

  if (is.null(referenceVar) && !is.null(clustersVar)) {
    annotation_heatmap <- ComplexHeatmap::rowAnnotation(
      Clusters = clusters_df$Clusters,
      col = clusters_colors,
      annotation_legend_param = list(
        title = "Clusters",
        title_gp = grid::gpar(fontsize = 11),
        labels_gp = grid::gpar(fontsize = 8),
        legend_height = grid::unit(3, "cm"),
        legend_width = grid::unit(1.5, "cm"),
        grid_height = grid::unit(0.6, "cm"),
        grid_width = grid::unit(0.6, "cm")
      )
    )
  }

  if (is.null(referenceVar) && is.null(clustersVar)) {
    annotation_heatmap <- NULL
  }

  if (!is.null(referenceVar) && !is.null(clustersVar)) {
    annotation_heatmap <- ComplexHeatmap::rowAnnotation(
      Annotations = annotation_df$Annotations,
      Clusters = clusters_df$Clusters,
      col = list(
        Annotations = annot_colors,
        Clusters = clusters_colors
      ),
      annotation_legend_param = list(
        Annotations = list(
          title = "Annotations",
          title_gp = grid::gpar(fontsize = 11),
          labels_gp = grid::gpar(fontsize = 8)
        ),
        Clusters = list(
          title = "CNV cluster",
          title_gp = grid::gpar(fontsize = 11),
          labels_gp = grid::gpar(fontsize = 8)
        )
      )
    )
  }

  hm <-  ComplexHeatmap::Heatmap(
    M,
    right_annotation = annotation_heatmap,
    border = TRUE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    use_raster = TRUE,
    clustering_distance_rows = "euclidean",
    clustering_method_rows = "ward.D",
    column_split = as.numeric(sapply(strsplit(rownames(as.matrix(Seurat::GetAssay(seuratObj, assay = "genomicScores")["data"])), ".", fixed = TRUE), function(z) z[1])),
    column_title_gp = grid::gpar(fontsize = 8),
    column_title = c(1:22,"X"),
    row_split = as.factor(split_df[[1]]),
    row_title = NULL,
    col = circlize::colorRamp2(c(-1, -0.6, -0.3, 0, 0.3, 0.6, 1), c("#0B2F7EFF", "#2A4D9EFF", "#A0A0FFFF", "white", "#E3807D", "#A4161A","#7A0A0D")),
    heatmap_legend_param = list(
      title = "CNV",
      title_gp = grid::gpar(fontsize = 11),
      labels_gp = grid::gpar(fontsize = 8),
      grid_height = grid::unit(1, "cm"),
      grid_width = grid::unit(0.6, "cm"))
    )

  if(printPlot == TRUE) {
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow = 2, heights = grid::unit.c(grid::unit(1, "cm"), grid::unit(1, "null")))))
    grid::pushViewport(grid::viewport(layout.pos.row = 1))
    grid::grid.text(paste0("CNV heatmap for sample ", seuratObj@project.name), gp = grid::gpar(fontsize = 20))
    grid::popViewport()
    grid::pushViewport(grid::viewport(layout.pos.row = 2))
    ComplexHeatmap::draw(hm, newpage = FALSE)
    grid::popViewport()
  }

  if(!is.null(savePath)) {
    if (outputType == "png") {
      fname <- file.path(savePath, paste0("heatmap.fastCNV_",seuratObj@project.name,".png"))
      png(width = 3500, height = 2200, filename = fname, res = 300)
    }
    if (outputType == "pdf"){
      fname <- file.path(savePath, paste0("heatmap.fastCNV_",seuratObj@project.name,".pdf"))
      pdf(width = 12, height = 7, file = fname)
    }
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow = 2, heights = grid::unit.c(grid::unit(1, "cm"),grid::unit(1, "null")))))
    grid::pushViewport(grid::viewport(layout.pos.row = 1))
    grid::grid.text(paste0("CNV heatmap for sample ", seuratObj@project.name), gp = grid::gpar(fontsize = 20))
    grid::popViewport()
    grid::pushViewport(grid::viewport(layout.pos.row = 2))
    ComplexHeatmap::draw(hm, newpage = FALSE)
    grid::popViewport()
    Sys.sleep(3)
    dev.off()
    message("CNV plot for sample ",seuratObj@project.name, " saved at ", fname)
  }
  #return(hm)
}
