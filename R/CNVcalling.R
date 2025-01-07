#' CNVcalling
#' Runs a CNV analysis given a seurat object
#'
#' @param seuratObj The Seurat object containing the data to analyze
#' @param assay Name of the assay to run the CNV on. Takes the results of `prepareCountsForCNVAnalysis` by default if available
#' @param referenceVar The variable name of the annotations in the Seurat metadata
#' @param referenceLabel The label given to the observations wanted as reference (can be any type of annotation)
#' @param scaleOnReferenceLabel If you want to scale the results depending on the reference observations
#' @param denoise If the data needs to be denoised (default = `TRUE`)
#' @param thresholdPercentile which quantiles to take (if 0.01 it will take 0.01-0.99). Background noise appears with higher numbers.
#' @param geneMetadata Dataframe of genes and their metadata from ensembl
#' @param genesToForce List of genes to force into the analysis (by default, fastCNV only takes the `topNgenes` most expressed).
#' @param windowSize Size of the genomic windows
#' @param windowStep Step between the genomic windows
#' @param saveGenomicWindows If the information of the genomic windows need to be saved in the current directory (default = `FALSE`).
#' @param topNGenes Number of top expressed genes to keep
#'
#' @return This function returns the genomic scores per genomic window
#'
#' @import stats
#' @import Seurat
#' @import SeuratObject
#'
#' @export
#'
#'
CNVcalling <- function(seuratObj,
                       assay = NULL,
                       referenceVar = NULL,
                       referenceLabel = NULL,
                       scaleOnReferenceLabel = TRUE,
                       denoise = TRUE,
                       thresholdPercentile = 0.01,
                       geneMetadata=getGenes(),
                       genesToForce = NULL,
                       windowSize=100,
                       windowStep=20,
                       saveGenomicWindows = FALSE,
                       topNGenes=7000) {
  if (dim(seuratObj)[1] < topNGenes) {topNGenes = dim(seuratObj)[1]}
  # getting reference cells / spots
  if (is.null(referenceVar) || is.null(referenceLabel)){
    message(paste0("referenceVar and/or referenceLabel parameters not found. Computing the CNV without a reference."))
    # unable to scale the results on reference data if we don't know what the reference data is
    scaleOnReferenceLabel = FALSE
  } else {
    if (length(referenceLabel) == 1) {
      referenceCells <- Seurat::Cells(seuratObj)[which(Seurat::FetchData(seuratObj, vars = referenceVar) == referenceLabel)]
      if (length(referenceCells) == 0) {
        message(paste0("CNVcalling : there is no annotation called ", referenceLabel," in the ", referenceVar," metadata slot of your seurat object.
Computing the CNV without a reference."))
        scaleOnReferenceLabel = FALSE
        rm(referenceCells)
      }
    } else {
      referenceCells <- list()
      for (i in referenceLabel){
        if (length(Seurat::Cells(seuratObj)[which(Seurat::FetchData(seuratObj, vars = referenceVar) == i)]) >= 5) {
                  referenceCells[[i]] <- Seurat::Cells(seuratObj)[which(Seurat::FetchData(seuratObj, vars = referenceVar) == i)]
        }
      }
      numberOfReferenceCells <- sum(sapply(referenceCells, length))
      if (numberOfReferenceCells == 0) {
        message("Couldn't find any cells annotated as ",referenceLabel," in the ",referenceVar," metadata slot of the seurat object.
                Computing the CNV without a reference.")
        scaleOnReferenceLabel = FALSE
        rm(numberOfReferenceCells,referenceCells)
      }
    }
  }

  # preparation of gene information
  geneMetadata <- geneMetadata[which(geneMetadata$gene_biotype %in% c("protein_coding","lncRNA") & geneMetadata$chromosome_name %in% c(1:22,"X") & geneMetadata$hgnc_symbol !=""),]
  geneMetadata$chromosome_num <- geneMetadata$chromosome_name
  geneMetadata$chromosome_num[which(geneMetadata$chromosome_num=="X")]<- 23
  geneMetadata$chromosome_num <- as.numeric(geneMetadata$chromosome_num)
  geneMetadata2 <- unique(geneMetadata[,c("hgnc_symbol","chromosome_num","start_position", "chr_arm")])

  # internal functions
  funTrim <- function(normcounts,lo=-3,up=3){
    t(apply(normcounts,1, function(z) {
      z[which(z < lo)] <- lo ; z[which(z > up)]<-up;z } ))
  }


  # analysis
  if (is.null(assay)) {
    if ("AggregatedCounts" %in% Seurat::Assays(seuratObj)) {
      assay = "AggregatedCounts"
    } else {
      assay = Seurat::Assays(seuratObj)[1]
    }
  }

  rawCounts <- as.matrix(Seurat::GetAssay(seuratObj, assay = assay)["counts"])
  commonGenes <- intersect(rownames(rawCounts),geneMetadata2$hgnc_symbol)
  rawCounts <- rawCounts[commonGenes,]

  if (!is.null(referenceVar) && !is.null(referenceLabel)){
    averageExpression <- rowMeans(rawCounts[,unlist(referenceCells)])
  } else {
    averageExpression <- rowMeans(rawCounts)
  }
  topExprGenes <- commonGenes[order(averageExpression, decreasing = T)[1:topNGenes]]

  topExprGenes_metadata <- geneMetadata2[geneMetadata2$hgnc_symbol %in% topExprGenes, ]
  topExprGenes_metadata$chr_arm_full <- paste0(topExprGenes_metadata$chromosome_num, topExprGenes_metadata$chr_arm)

  genes_by_arm <- split(
    topExprGenes_metadata$hgnc_symbol,
    topExprGenes_metadata$chr_arm_full
  )

  for (arm in unique(geneMetadata2$chr_arm)) {
    if (!(arm %in% names(genes_by_arm))) {
      genes_by_arm[[arm]] <- character(0)
    }
    if (length(genes_by_arm[[arm]]) < 200) {
      remaining_genes <- commonGenes[!commonGenes %in% genes_by_arm[[arm]]]
      top_arm_genes <- remaining_genes[order(averageExpression[commonGenes %in% remaining_genes], decreasing = TRUE)[1:200]]
      genes_by_arm[[arm]] <- unique(c(genes_by_arm[[arm]], top_arm_genes))
    }
  }
  final_selected_genes <- unlist(genes_by_arm)

  rawCounts <- rawCounts[topExprGenes,]

  normCounts <- log2(1+rawCounts)
  normCounts <- scale(normCounts, scale = FALSE)

  if (scaleOnReferenceLabel) {
    if (length(referenceLabel) == 1) {
      scaleFactor <- rowMeans(normCounts[,referenceCells])
    } else {
      scaleFactor <- list()
      for (i in referenceLabel) {
        scaleFactor[[i]] <- rowMeans(normCounts[,referenceCells[[i]]])
      }
      scaleFactor <- do.call(rbind,scaleFactor)
      scaleFactor <- na.omit(scaleFactor)
      scaleFactor <- apply(scaleFactor, 2, median)
    }
  } else {
    scaleFactor <- rowMeans(normCounts)
  }

  normCounts <- normCounts - scaleFactor
  normCounts <- funTrim(normCounts, lo = -3, up = 3)

  geneMetadata2 <- geneMetadata2[which(geneMetadata2$hgnc_symbol %in% topExprGenes),]
  geneMetadata2 <- geneMetadata2[order(geneMetadata2$chromosome_num,geneMetadata2$start_position),]

  # Preparation of genomic windows
  genomicWindows <- lapply(c(1:23), function(chrom) {
    genesC <- geneMetadata2[which(geneMetadata2$chromosome_num == chrom),]
    chr_arms <- unique(genesC$chr_arm)
    chrom_windows <- list()
    for (arm in chr_arms) {
      genesArm <- genesC[which(genesC$chr_arm == arm),]
      N <- nrow(genesArm)
      iter <- round(windowSize / 2)
      if (N > windowSize) {
        gw <- lapply(seq(iter + 1, N - iter, by = windowStep), function(i) {
          genesArm[(i - iter):(i + iter), "hgnc_symbol"] |> unlist() |> as.character()
        })
        names(gw) <- paste0(chrom, ".", arm, 1:length(gw))
      } else {
        gw <- list(as.character(unlist(genesArm$"hgnc_symbol")))
        names(gw) <- paste0(chrom, ".", arm, 1)
      }
      chrom_windows <- c(chrom_windows, gw)
    }
    return(chrom_windows)
  })

  genomicWindows <- unlist(genomicWindows,recursive=F)
  genomicScores <- sapply(genomicWindows, function(g) colMeans(normCounts[g,]) )

  if (denoise) {
    if (scaleOnReferenceLabel) {
      if (length(referenceLabel) == 1) {
        genomicScoresReferenceLabel<- genomicScores[referenceCells,]
      } else {
        genomicScoresReferenceLabel <- list()
        for (i in referenceLabel){
          genomicScoresReferenceLabel[[i]] <- genomicScores[referenceCells[[i]],]
        }
        genomicScoresReferenceLabel <- do.call(rbind, genomicScoresReferenceLabel)
      }

    Q01Q99 <- apply(genomicScoresReferenceLabel, 2, stats::quantile, "probs"=c(0+thresholdPercentile,1-thresholdPercentile))
    genomicScoresTrimmed <- apply(genomicScores, 1,function(v)
      {v[which(v >= Q01Q99[1,] & v <= Q01Q99[2,])] <- 0 ; v })

    } else {
      if (is.null(referenceVar)) {
        Q01Q99 <- apply(genomicScores,2,stats::quantile,"probs"=c(0+thresholdPercentile,1-thresholdPercentile))
        genomicScoresTrimmed <- apply(genomicScores,1,function(v)
          {v[which(v >= Q01Q99[1,] & v <= Q01Q99[2,])] <- 0;v})
      } else {
        cellLines <- split(Seurat::Cells(seuratObj), Seurat::FetchData(seuratObj, vars = referenceVar))
        high_threshold <- median(unlist(sapply(cellLines, function(z) apply (genomicScores[z,], 2, function(v) quantile(v, probs = c(0.99))))))
        low_threshold <- median(unlist(sapply(cellLines, function(z) apply (genomicScores[z,], 2, function(v) quantile(v, probs = c(0.01))))))
        genomicScoresTrimmed <- t(apply(genomicScores, 2, function(z){
          z[which( z>low_threshold & z<high_threshold)] = 0; z}))
      }
    }
  } else {
    genomicScoresTrimmed <- t(genomicScores)
  }

  if (saveGenomicWindows){
    save(genomicWindows, file = paste0("genomicWindows_size",windowSize,"_step",windowStep,".RData"))
  }

  genomicAssay <- Seurat::CreateAssayObject(data = as.matrix(genomicScoresTrimmed))
  seuratObj[["genomicScores"]] <- genomicAssay
  seuratObj[["cnv_fraction"]] <- colMeans(abs(genomicScoresTrimmed) > 0)

  return (seuratObj)
}
