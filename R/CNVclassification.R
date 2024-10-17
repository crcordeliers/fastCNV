#' CNVclassification
#'
#' @param seuratObj A seurat object containing the results of fastCNV
#' @param tumorLabels A list containing the tumor labels
#'
#' @return A seurat object containing the classification (loss gain no_alteration) for
#' each observation and each chromosome arm.
#'
#' @import scales
#'
#' @export
#'

CNVclassification <- function(seuratObj, tumorLabels) {
  chromosomes <- c(paste0(1:22, "p_CNV"), paste0(1:22, "q_CNV"))
  chromosomes <- c(chromosomes, "X_CNV")

  metadata <- Seurat::FetchData(seuratObj, vars = colnames(seuratObj@meta.data))

  cnvVectorAll <- c()

  for (chrom in chromosomes) {
    tumor_cells <- metadata$annot %in% tumorLabels
    cnvVectorAll <- c(cnvVectorAll, metadata[[chrom]][tumor_cells])
  }

  peaks <- cit.peaks(cnvVectorAll, maxNbPeaks = 3, percentHighestPeak = 0.02)

  for (chrom in chromosomes) {
    cnvVector <- metadata[[chrom]]
    classification <- classify_cnv(cnvVector, peaks)
    metadata[[paste0(chrom, "_classification")]] <- classification
  }

  seuratObj@meta.data <- metadata

  return (seuratObj)
}
