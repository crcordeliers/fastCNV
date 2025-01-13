#' A function to plot the CNV results into a heatmap
#'
#' @param seuratObj Seurat object containing the data used to get the genomicScores.
#' @param referenceVar The variable name of the annotations in the Seurat metadata
#' @param splitPlotOnVar The variable name on which to split the heatmap rows.
#' @param savePath Path to save the pdf heatmap. If `NULL`, plot won't be saved (default = `.`).
#' @param printPlot If the heatmap should be printed in the console.
#' @param downsizePlot Subset the observations to speed up the plotting process (default = `FALSE`).
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
#'
#' @return This function builds a heatmap and saves it as a .pdf file in your working directory
#'
#' @export
#'

plotCNVResults <- function(seuratObj, referenceVar = NULL,
                           splitPlotOnVar = referenceVar, savePath = ".",
                           printPlot = FALSE, downsizePlot = FALSE){
  M <- t(as.matrix(Seurat::GetAssay(seuratObj, "genomicScores")["data"]))
  if (downsizePlot == TRUE && !is.null(referenceVar)) {
    # K-means clustering
    LClust <- kmeans(M, centers = 40)
    LClust <- split(Seurat::Cells(seuratObj), LClust$cluster)
    LClust <- lapply(LClust, function(x) split(x, Seurat::FetchData(seuratObj, vars = referenceVar, cells = x)))

    centroid_list <- list()
    annotation_list <- list()
    split_list <- list()

    for (cluster in names(LClust)) {
      for (annotation in names(LClust[[cluster]])) {
        rows <- intersect(rownames(M), LClust[[cluster]][[annotation]])
        if (length(rows) > 0) {
          subset_matrix <- M[rows, , drop=FALSE]
          if (nrow(subset_matrix) > 1) {
            centroid <- colMeans(subset_matrix)
            centroid_list[[paste0(cluster, "_", annotation)]] <- centroid
            annotation_list[[paste0(cluster, "_", annotation)]] <- annotation
            if (!is.null(splitPlotOnVar)) {
              split_values <- seuratObj@meta.data[LClust[[cluster]][[annotation]], splitPlotOnVar, drop = TRUE]
              split_list[[paste0(cluster, "_", annotation)]] <- names(sort(table(split_values), decreasing = TRUE))[1]  # Take the most common value
            }
          } else if (nrow(subset_matrix) == 1) {
            centroid_list[[paste0(cluster, "_", annotation)]] <- subset_matrix
            annotation_list[[paste0(cluster, "_", annotation)]] <- annotation
            if (!is.null(splitPlotOnVar)) {
              split_list[[paste0(cluster, "_", annotation)]] <- seuratObj@meta.data[rownames(subset_matrix), splitPlotOnVar]
            }
          }
        }
      }
    }
    annotation_df <- data.frame(Annotation = unlist(annotation_list))
    M <- as.matrix(do.call(rbind, centroid_list))
    rownames(M) <- names(centroid_list)
    if (!is.null(splitPlotOnVar)) {
      split_df <- data.frame(Split = unlist(split_list))
    } else {
      split_df <- NULL
    }
  } else if (downsizePlot == FALSE && !is.null(referenceVar)) {
    annotation_df <- as.data.frame(seuratObj@meta.data[[referenceVar]])
    colnames(annotation_df) <- "Annotations"
    if (!is.null(splitPlotOnVar)) {
      split_df <- as.data.frame(Seurat::FetchData(seuratObj, vars = splitPlotOnVar))
      colnames(split_df) <- "Split"
    } else {
      split_df <- NULL
    }
  } else if (downsizePlot == TRUE && is.null(referenceVar)) {
    LClust <- kmeans(M, centers = 150)
    LClust <- split(Seurat::Cells(seuratObj), LClust$cluster)
    centroid_list <- list()
    split_list <- list()
    for (cluster in names(LClust)) {
      cluster_rows <- unlist(LClust[[cluster]])
      rows <- intersect(rownames(M), cluster_rows)
      if (length(rows) > 0) {
        subset_matrix <- M[rows, , drop = FALSE]
        if (nrow(subset_matrix) > 1) {
          centroid <- colMeans(subset_matrix)
          centroid_list[[cluster]] <- centroid
          if (!is.null(splitPlotOnVar)) {
            split_values <- seuratObj@meta.data[cluster_rows, splitPlotOnVar, drop = TRUE]
            split_list[[cluster]] <- names(sort(table(split_values), decreasing = TRUE))[1]  # Take the most common value
          }
        } else if (nrow(subset_matrix) == 1) {
          centroid_list[[cluster]] <- as.numeric(subset_matrix)
          if (!is.null(splitPlotOnVar)) {
            split_list[[cluster]] <- seuratObj@meta.data[rownames(subset_matrix), splitPlotOnVar]
          }
        }
      }
    }
    M <- as.matrix(do.call(rbind, centroid_list))
    if (!is.null(splitPlotOnVar)) {
      split_df <- data.frame(Split = unlist(split_list))
    } else {
      split_df <- NULL
    }
  } else if (downsizePlot == FALSE && is.null(splitPlotOnVar)) {
    split_df <- NULL
  }

  if(is.null(referenceVar)) {
    if (is.null(splitPlotOnVar)) {splitPlot = NULL}else{splitPlot = as.factor(split_df[[1]])}
    hm <-  ComplexHeatmap::Heatmap(
      M,
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
      row_split = splitPlot,
      col = circlize::colorRamp2(c(-1, -0.6, -0.3, 0, 0.3, 0.6, 1), c("#0B2F7EFF", "#2A4D9EFF", "#A0A0FFFF", "white", "#E3807D", "#A4161A","#7A0A0D")),
      heatmap_legend_param = list(
        title = "CNV",
        title_gp = grid::gpar(fontsize = 11),  # Font size for the title
        labels_gp = grid::gpar(fontsize = 8),  # Font size for the labels
        grid_height = grid::unit(1, "cm"),  # Size of the grid cells in the legend
        grid_width = grid::unit(0.6, "cm"))
    )
  } else {
    if (referenceVar == splitPlotOnVar) {
      rown = NULL
    } else {
      rown = character(0)
    }
    palette <- as.character(paletteer::paletteer_d("pals::glasbey"))
    if (length(unique(annotation_df$Annotation)) > length(palette)) {
      hm <-  ComplexHeatmap::Heatmap(
        M,
        right_annotation = ComplexHeatmap::rowAnnotation(
          Annotations = annotation_df$Annotations,
          annotation_legend_param = list(
            title = "Annotations",  # Set the title of the annotation legend
            title_gp = grid::gpar(fontsize = 11),  # Font size for the title
            labels_gp = grid::gpar(fontsize = 8),  # Font size for the labels
            legend_height = grid::unit(3, "cm"),  # Height of the annotation legend
            legend_width = grid::unit(1.5, "cm"),  # Width of the annotation legend
            grid_height = grid::unit(0.6, "cm"),  # Size of the grid cells in the legend
            grid_width = grid::unit(0.6, "cm")  # Size of the grid cells in the legend
          )
        ),
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
        row_title = rown,
        col = circlize::colorRamp2(c(-1, -0.6, -0.3, 0, 0.3, 0.6, 1), c("#0B2F7EFF", "#2A4D9EFF", "#A0A0FFFF", "white", "#E3807D", "#A4161A","#7A0A0D")),
        heatmap_legend_param = list(
          title = "CNV",
          title_gp = grid::gpar(fontsize = 11),  # Font size for the title
          labels_gp = grid::gpar(fontsize = 8),  # Font size for the labels
          grid_height = grid::unit(1, "cm"),  # Size of the grid cells in the legend
          grid_width = grid::unit(0.6, "cm"))
      )
    } else {
      annot_colors <- setNames(palette[1:length(unique(annotation_df$Annotation))], unique(annotation_df$Annotation))
      hm <-  ComplexHeatmap::Heatmap(
        M,
        right_annotation = ComplexHeatmap::rowAnnotation(
          Annotations = annotation_df$Annotations,
          col = list(Annotations = annot_colors),
          annotation_legend_param = list(
            title = "Annotations",  # Set the title of the annotation legend
            title_gp = grid::gpar(fontsize = 11),  # Font size for the title
            labels_gp = grid::gpar(fontsize = 8),  # Font size for the labels
            legend_height = grid::unit(3, "cm"),  # Height of the annotation legend
            legend_width = grid::unit(1.5, "cm"),  # Width of the annotation legend
            grid_height = grid::unit(0.6, "cm"),  # Size of the grid cells in the legend
            grid_width = grid::unit(0.6, "cm")  # Size of the grid cells in the legend
          )
        ),
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
        row_title = rown,
        col = circlize::colorRamp2(c(-1, -0.6, -0.3, 0, 0.3, 0.6, 1), c("#0B2F7EFF", "#2A4D9EFF", "#A0A0FFFF", "white", "#E3807D", "#A4161A","#7A0A0D")),
        heatmap_legend_param = list(
          title = "CNV",
          title_gp = grid::gpar(fontsize = 11),  # Font size for the title
          labels_gp = grid::gpar(fontsize = 8),  # Font size for the labels
          grid_height = grid::unit(1, "cm"),  # Size of the grid cells in the legend
          grid_width = grid::unit(0.6, "cm"))
      )
    }
  }
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
    fname <- file.path(savePath, paste0("heatmap.fastCNV_", splitPlotOnVar,"_",seuratObj@project.name,".pdf"))
    pdf(width = 12, height = 7, file = fname)
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
