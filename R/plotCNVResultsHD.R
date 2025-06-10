#' Plot Visium HD CNV Results into a Heatmap
#' Builds a heatmap to visualize the Visium HD CNV results based on genomic scores.
#'
#' @param seuratObjHD A Seurat object containing the genomic scores computed previously.
#' @param referenceVar The name of the metadata column in the Seurat object containing reference annotations.
#' @param clustersVar The name of the metadata column containing cluster information (default = `"cnv_clusters"`).
#' @param splitPlotOnVar The name of the metadata column used to split the heatmap rows (e.g., cell type or cluster) (default = `clustersVar`).
#' @param savePath The path where the heatmap will be saved. If `NULL`, the plot will not be saved (default = `"."`).
#' @param printPlot Logical. If `TRUE`, prints the heatmap to the console.
#' @param referencePalette A color palette for `referenceVar`.
#' You can provide a custom palette as a vector of color codes (e.g., `c("#FF0000", "#00FF00")`).
#' @param clusters_palette A color palette for `clustersVar`.
#' You can provide a custom palette as a vector of color codes (e.g., `c("#F8766D", "#A3A500", "#00BF7D")`).
#' @param outputType Character. Specifies the file format for saving the plot, either `"png"` or `"pdf"`.
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
#' @return This function generates a heatmap and saves it as a `.pdf` or `.png` file in the specified path (default = working directory).
#'
#' @export


plotCNVResultsHD <- function(seuratObjHD,
                             referenceVar = NULL,
                             clustersVar = "cnv_clusters",
                             splitPlotOnVar = clustersVar,
                             savePath = ".",
                             printPlot = FALSE,
                             referencePalette = "default",
                             clusters_palette = "default",
                             outputType = "png"){
  if (outputType != "png" && outputType != "pdf"){
    message("Warning : outputType not valid, should be 'pdf' or 'png'. Setting outputType to 'png'")
    outputType = "png"
  }

  if (clustersVar == "cnv_clusters"){
    if (!(clustersVar %in% names(seuratObjHD@meta.data))) {
      if(splitPlotOnVar == clustersVar){
        splitPlotOnVar = referenceVar
      }
      clustersVar = NULL
    }
  }

  M <- t(as.matrix(Seurat::GetAssay(seuratObjHD, "genomicScores")["data"]))
  if (any(referencePalette == "default")) {
    referencePalette = as.character(paletteer::paletteer_d("pals::glasbey"))
  }
  if (!is.null(referenceVar)) {
    annotation_df <- as.data.frame(seuratObjHD@meta.data[[referenceVar]])
    colnames(annotation_df) <- "Annotations"
    annot_colors <- setNames(referencePalette[1:length(unique(annotation_df$Annotations))], unique(annotation_df$Annotations))
    if (!is.null(splitPlotOnVar)) {
      split_df <- as.data.frame(Seurat::FetchData(seuratObjHD, vars = splitPlotOnVar))
      colnames(split_df) <- "Split"
    } else {
      split_df <- NULL
    }
  } else if (is.null(splitPlotOnVar)) {
    split_df <- NULL
  }
  if (!is.null(clustersVar)) {
    clusters_df <- as.data.frame(seuratObjHD@meta.data[[clustersVar]])
    colnames(clusters_df) <- "Clusters"
    if (any(clusters_palette == "default")){
      clusters_palette = scales::hue_pal()(length(unique(clusters_df$Clusters)))
    }
    clusters_colors <- setNames(clusters_palette, unique(clusters_df$Clusters))
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
    cluster_rows = F,
    show_row_names = FALSE,
    show_column_names = FALSE,
    use_raster = TRUE,
    clustering_distance_rows = "euclidean",
    clustering_method_rows = "ward.D",
    column_split = as.numeric(sapply(strsplit(rownames(as.matrix(Seurat::GetAssay(seuratObjHD, assay = "genomicScores")["data"])), ".", fixed = TRUE), function(z) z[1])),
    column_title_gp = grid::gpar(fontsize = 8),
    column_title = c(1:22,"X"),
    row_split = as.factor(split_df[[1]]),
    row_title = NULL,
    col = circlize::colorRamp2(c(-0.2, -0.1, -0.05, 0, 0.05, 0.1, 0.2), c("#0B2F7EFF", "#2A4D9EFF", "#A0A0FFFF", "white", "#E3807D", "#A4161A","#7A0A0D")),
    heatmap_legend_param = list(
      title = "CNV Score",
      title_gp = grid::gpar(fontsize = 11),
      labels_gp = grid::gpar(fontsize = 8),
      grid_height = grid::unit(1, "cm"),
      grid_width = grid::unit(0.6, "cm"))
  )

  if(printPlot == TRUE) {
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow = 2, heights = grid::unit.c(grid::unit(1, "cm"), grid::unit(1, "null")))))
    grid::pushViewport(grid::viewport(layout.pos.row = 1))
    grid::grid.text(paste0("CNV heatmap for sample ", seuratObjHD@project.name), gp = grid::gpar(fontsize = 20))
    grid::popViewport()
    grid::pushViewport(grid::viewport(layout.pos.row = 2))
    ComplexHeatmap::draw(hm, newpage = FALSE)
    grid::popViewport()
  }

  if(!is.null(savePath)) {
    if (outputType == "png") {
      fname <- file.path(savePath, paste0("heatmap.fastCNV_",seuratObjHD@project.name,".png"))
      png(width = 3500, height = 2200, filename = fname, res = 300)
    }
    if (outputType == "pdf"){
      fname <- file.path(savePath, paste0("heatmap.fastCNV_",seuratObjHD@project.name,".pdf"))
      pdf(width = 12, height = 7, file = fname)
    }
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow = 2, heights = grid::unit.c(grid::unit(1, "cm"),grid::unit(1, "null")))))
    grid::pushViewport(grid::viewport(layout.pos.row = 1,gp = grid::gpar(fill = "white")))
    grid::grid.text(paste0("CNV heatmap for sample ", seuratObjHD@project.name), gp = grid::gpar(fontsize = 20))
    grid::popViewport()
    grid::pushViewport(grid::viewport(layout.pos.row = 2))
    ComplexHeatmap::draw(hm, newpage = FALSE)
    grid::popViewport()
    Sys.sleep(3)
    dev.off()
    message("CNV plot for sample ",seuratObjHD@project.name, " saved at ", fname)
  }
}
