#' CNVPerChromosomeArm
#' Computes the CNV fraction of each spot/cell per chromosome arm, then stocks the results into the metadata
#'
#' @param seuratObject The output of fastCNV()
#'
#' @return This function returns the same seurat object with the CNV for each chromosome arm in the metadata
#' @export
#'

CNVPerChromosomeArm <- function(seuratObject) {
  genomicScores <- as.matrix(Seurat::GetAssay(seuratObject, assay = "genomicScores")["data"])

  proportions <- list(
    '1' = c(0.48, 0.52),
    '2' = c(0.42, 0.58),
    '3' = c(0.38, 0.62),
    '4' = c(0.36, 0.64),
    '5' = c(0.37, 0.63),
    '6' = c(0.40, 0.60),
    '7' = c(0.43, 0.57),
    '8' = c(0.45, 0.55),
    '9' = c(0.47, 0.53),
    '10' = c(0.50, 0.50),
    '11' = c(0.47, 0.53),
    '12' = c(0.46, 0.54),
    '13' = c(0.15, 0.85),
    '14' = c(0.15, 0.85),
    '15' = c(0.15, 0.85),
    '16' = c(0.49, 0.51),
    '17' = c(0.51, 0.49),
    '18' = c(0.25, 0.75),
    '19' = c(0.53, 0.47),
    '20' = c(0.49, 0.51),
    '21' = c(0.33, 0.67),
    '22' = c(0.33, 0.67))

  for (chrom in names(proportions)) {
    chrom_rows <- grep(paste0("^", chrom, "\\."), rownames(genomicScores), value = TRUE)
    total_windows <- length(chrom_rows)

    p_size <- floor(proportions[[chrom]][1] * total_windows)
    q_size <- total_windows - p_size


    if (p_size == 0) {
      p_size <- 1
      q_size <- total_windows - p_size
    }
    if (q_size == 0) {
      q_size <- 1
      p_size <- total_windows - q_size
    }


    if (p_size > 1) {
      seuratObject[[paste0(chrom, "p_CNV")]] <- apply(genomicScores[chrom_rows[1:p_size], ], 2, function(x) {
        non_zero_values <- x
        if (length(non_zero_values) > 0) {
          mean(non_zero_values)
        } else {
          0
        }
      })
    } else {
      seuratObject[[paste0(chrom, "p_CNV")]] <- genomicScores[chrom_rows[1:p_size],]
    }

    if (q_size > 1) {
      seuratObject[[paste0(chrom, "q_CNV")]] <- apply(genomicScores[chrom_rows[(p_size + 1):total_windows], ], 2, function(x) {
        non_zero_values <- x
        if (length(non_zero_values) > 0) {
          mean(non_zero_values)
        } else {
          0
        }
      })
    } else {
      seuratObject[[paste0(chrom, "q_CNV")]] <- genomicScores[chrom_rows[(p_size + 1):total_windows],]
    }
  }

  chrom_X <- grep(paste0("^", "23", "\\."), rownames(genomicScores), value = TRUE)
  if (length(chrom_X) == 1) {
    seuratObject[["X_CNV"]] <- genomicScores[chrom_X,]
  } else {
    seuratObject[["X_CNV"]] <- apply(genomicScores[chrom_X,], 2, function(x) {
      non_zero_values <- x
      if (length(non_zero_values) > 0) {
        mean(non_zero_values)
      } else {
        0
      }
    })
  }


  return(seuratObject)
}
