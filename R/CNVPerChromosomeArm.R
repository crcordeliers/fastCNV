#' CNVPerChromosomeArm
#' Computes the CNV fraction of each spot/cell per chromosome arm, then stocks the results into the metadata
#'
#' @param seuratObj The output of fastCNV()
#'
#' @return This function returns the same seurat object with the CNV for each chromosome arm in the metadata
#' @export
#'

CNVPerChromosomeArm <- function(seuratObj) {
  genomicScores <- as.matrix(Seurat::GetAssay(seuratObj, assay = "genomicScores")["data"])
  window_names <- rownames(genomicScores)

  extract_chrom_arm <- function(window_name) {
    chrom_arm <- sub("(\\d+\\.\\w).*", "\\1", window_name)
    chrom_arm <- sub("^23\\.(p|q)$", "X.\\1", chrom_arm)
    return(chrom_arm)
  }

  window_info <- data.frame(
    window = window_names,
    chrom_arm = sapply(window_names, extract_chrom_arm)
  )

  chrom_arms_all <- unique(c(
    paste0(1:22, ".p"), paste0(1:22, ".q"),
    "X.p", "X.q"
  ))

  arm_averages <- list()
  for (chrom_arm in chrom_arms_all) {
    if (chrom_arm %in% unique(window_info$chrom_arm)) {
      windows <- window_info$window[window_info$chrom_arm == chrom_arm]
      subset_scores <- genomicScores[windows, , drop = FALSE]
      avg_values <- colMeans(subset_scores, na.rm = TRUE)
    } else {
      avg_values <- rep(0, ncol(genomicScores))
    }
    arm_averages[[chrom_arm]] <- avg_values
  }

  arm_averages_df <- do.call(rbind, lapply(names(arm_averages), function(chrom_arm) {
    data.frame(
      chrom_arm = chrom_arm,
      value = arm_averages[[chrom_arm]],
      barcode = colnames(genomicScores),
      check.names = FALSE
    )
  }))

  for (chrom_arm in unique(arm_averages_df$chrom_arm)) {
    column_name <- paste0(chrom_arm, "_CNV")
    seuratObj[[column_name]] <- arm_averages_df$value[arm_averages_df$chrom_arm == chrom_arm]
  }

  return(seuratObj)
}
