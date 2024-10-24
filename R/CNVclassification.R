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

CNVclassification <- function(seuratObj, peaks = c(-0.1,0,0.1)) {
  chromosomes <- c(paste0(1:22, "p_CNV"), paste0(1:22, "q_CNV"))
  chromosomes <- c(chromosomes, "X_CNV")

  metadata <- Seurat::FetchData(seuratObj, vars = colnames(seuratObj@meta.data))


  for (chrom in chromosomes) {
    cnvVector <- metadata[[chrom]]
    classification <- classify_cnv(cnvVector, peaks)
    metadata[[paste0(chrom, "_classification")]] <- classification
  }

  seuratObj@meta.data <- metadata

  return (seuratObj)
}
