#' CNVanalysis
#' Runs the CNV functions given a seurat object or a list of seurat objects
#'
#' @param object The Seurat object or list of seurat objects containing the data to analyze
#' @param assay Name of the assay to run the CNV on. Takes the results of prepareCountsForCNVAnalysis by default if available
#' @param referenceVar The variable name of the annotations in the Seurat metadata
#' @param referenceLabel The label given to the observations wanted as reference (can be any type of annotation)
#' @param scaleOnReferenceLabel If you want to scale the results depending on the normal observations
#' @param thresholdPercentile Which quantiles to take (if 0.01 it will take 0.01-0.99). Background noise appears with higher numbers.
#' @param genes List of genes from ensembl
#' @param windowSize Size of the genomic windows
#' @param windowStep Step between the genomic windows
#' @param topNGenes Number of top expressed genes to keep
# @param doRecapPlot Default `TRUE`. Will output the CNV heatmaps by annotation if `TRUE`.
#' @param pooledReference Default `TRUE`. Will build a pooled reference across all samples if `TRUE`.
#'
#' @return This function returns the genomic scores per genomic window per seurat object. If given a list with more than one seurat object, and their annotations, it will output a heatmap per cell type given
#' @export
#'

CNVanalysis <- function(object,
                       referenceVar = NULL,
                       referenceLabel = NULL,
                       #doRecapPlot = TRUE,
                       pooledReference = TRUE,
                       scaleOnReferenceLabel = TRUE,
                       assay = NULL,
                       thresholdPercentile = 0.01,
                       genes=getGenes(),
                       windowSize=100,
                       windowStep=20,
                       topNGenes=7000) {

    if (!is.list(object)) {
      object <- CNVcalling(object,
                           assay = assay,
                           referenceVar = referenceVar,
                           referenceLabel = referenceLabel,
                           scaleOnReferenceLabel = scaleOnReferenceLabel,
                           thresholdPercentile = thresholdPercentile,
                           genes=genes,
                           windowSize=windowSize,
                           windowStep=windowStep,
                           topNGenes=topNGenes)
    } else {
      if (length(object) == 1) {
        object <- list(CNVcalling(object[1],
                                  assay = assay,
                                  referenceVar = referenceVar,
                                  referenceLabel = referenceLabel,
                                  scaleOnReferenceLabel = scaleOnReferenceLabel,
                                  thresholdPercentile = thresholdPercentile,
                                  genes=genes,
                                  windowSize=windowSize,
                                  windowStep=windowStep,
                                  topNGenes=topNGenes))
      } else {
        if (pooledReference == TRUE) {
          object <- CNVcallingList(object,
                             assay = assay,
                             referenceVar = referenceVar,
                             referenceLabel = referenceLabel,
                             scaleOnReferenceLabel = scaleOnReferenceLabel,
                             thresholdPercentile = thresholdPercentile,
                             genes=genes,
                             windowSize=windowSize,
                             windowStep=windowStep,
                             topNGenes=topNGenes)
        } else {
          object <- lapply(object, function(x) {
                                    CNVcalling(x,
                                    assay = assay,
                                    referenceVar = referenceVar,
                                    referenceLabel = referenceLabel,
                                    scaleOnReferenceLabel = scaleOnReferenceLabel,
                                    thresholdPercentile = thresholdPercentile,
                                    genes=genes,
                                    windowSize=windowSize,
                                    windowStep=windowStep,
                                    topNGenes=topNGenes) } )
        }

        # if (doRecapPlot == TRUE) {
        #   print("Plotting the CNV recap heatmap per category accross samples. This could take some time.")
        #   if (!is.null(referenceVar)) {
        #   palette <- as.character(paletteer::paletteer_d("ggsci::default_igv"))
        #   annot_colors <- setNames(palette[1:length(names(object))], names(object))
        #   uni <- list()
        #   for (i in names(object)){
        #     uni[[i]] <- unique(object[[i]]@meta.data[[referenceVar]])
        #   }
        #   cellTypes <- unique(unlist(uni))
        #   cellTypes <- cellTypes[nzchar(cellTypes)]
        #
        #   Lmat <- list()
        #   for (i in cellTypes){
        #     Lmat[[i]] <- lapply(object, function(x)
        #       Seurat::Cells(x)[which(x@meta.data[[referenceVar]] == i)])
        #     Lmat[[i]] <- Lmat[[i]][which(sapply(Lmat[[i]], length)>=5)]
        #     for (j in names(Lmat[[i]])) {
        #       Lmat[[i]][[j]] <- as.matrix(Seurat::GetAssay(object[[j]], assay = "genomicScores")$counts)[,Lmat[[i]][[j]]]
        #     }
        #   }
        #   Lhm <- list()
        #   for (x in names(Lmat)) {
        #     modified <- list()
        #     sample_info <- data.frame(col_name = character(), sample = character(), stringsAsFactors = FALSE)
        #     for (sample in names(Lmat[[x]])) {
        #       matrix <- Lmat[[x]][[sample]]
        #       col_names <- colnames(matrix)
        #       new_col_names <- paste0(sample, "_", col_names)
        #       colnames(matrix) <- new_col_names
        #       modified[[sample]] <- rowMeans(matrix)
        #
        #       sample_info <- rbind(sample_info, data.frame(col_name = new_col_names, sample = sample))
        #     }
        #     rownames(sample_info) <- sample_info[,1]
        #     sample_info[,1] <- NULL
        #     combined <- t(do.call(cbind, modified))
        #
        #     sample_df <- as.data.frame(unique(sample_info))
        #     rownames(sample_df) <- sample_df[,1]
        #
        #     Lhm[[x]] <- ComplexHeatmap::pheatmap(combined,
        #                                          border=T,
        #                                          border_color = NA,
        #                                          use_raster = F,
        #                                          cluster_cols = F,
        #                                          show_rownames = F,
        #                                          show_colnames = F,
        #                                          clustering_distance_rows = "euclidean",
        #                                          clustering_method = "ward.D",
        #                                          column_split = as.numeric(sapply(strsplit(colnames(combined),".",fixed=T),function(z)z[1])),
        #                                          row_split = sample_df,
        #                                          annotation_row = sample_df,
        #                                          annotation_colors = list(sample = annot_colors),
        #                                          col=circlize::colorRamp2(c(-1, 0, 1), c("darkblue", "white", "darkred")))
        #
        #   }
        #   output_file <- "heatmaps_per_annotation.pdf"
        #
        #   pdf(output_file, width = 20, height = 10)
        #   for (i in seq_along(Lhm)){
        #     grid::grid.newpage()
        #     grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow = 2, heights = grid::unit.c(grid::unit(1, "cm"),grid::unit(1, "null")))))
        #     grid::pushViewport(grid::viewport(layout.pos.row = 1))
        #     grid::grid.text(paste0("CNV heatmap for annotation : ", names(Lhm)[i]), gp = grid::gpar(fontsize = 20))
        #     grid::popViewport()
        #     grid::pushViewport(grid::viewport(layout.pos.row = 2))
        #     ComplexHeatmap::draw(Lhm[[i]], newpage = FALSE)
        #     grid::popViewport()
        #
        #   }
        #   dev.off()
        # }
        # }

      }
    }
  return (object)
}
