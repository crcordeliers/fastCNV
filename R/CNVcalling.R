#' CNVcalling
#' Runs a CNV analysis given a seurat object
#'
#' @param seuratObj The Seurat object containing the data to analyze
#' @param assay Name of the assay to run the CNV on. Takes the results of `prepareCountsForCNVAnalysis` by default if available
#' @param referenceVar The variable name of the annotations in the Seurat metadata
#' @param referenceLabel The label given to the observations wanted as reference (can be any type of annotation)
#' @param scaleOnReferenceLabel If you want to scale the results depending on the reference observations
#' @param thresholdPercentile which quantiles to take (if 0.01 it will take 0.01-0.99). Background noise appears with higher numbers.
#' @param genes List of genes from ensembl
#' @param windowSize Size of the genomic windows
#' @param windowStep Step between the genomic windows
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
                       thresholdPercentile = 0.01,
                       genes=getGenes(),
                       windowSize=100,
                       windowStep=20,
                       topNGenes=7000) {
  if (dim(seuratObj)[1] < topNGenes) {topNGenes = dim(seuratObj)[1]}
  # getting reference cells / spots
  if (is.null(referenceVar) || is.null(referenceLabel)){
    print(paste0("referenceVar and/or referenceLabel parameters not found. Computing the CNV without a reference."))
    # unable to scale the results on reference data if we don't know what the reference data is
    scaleOnReferenceLabel = FALSE
  } else {
    if (length(referenceLabel) == 1) {
      referenceCells <- Seurat::Cells(seuratObj)[which(Seurat::FetchData(seuratObj, vars = referenceVar) == referenceLabel)]
      if (length(referenceCells) == 0) {
        print(paste0("CNVcalling : there is no annotation called ", referenceLabel," in the ", referenceVar," metadata slot of your seurat object.
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
        print("Couldn't find any cells annotated as ",referenceLabel," in the ",referenceVar," metadata slot of the seurat object.
                Computing the CNV without a reference.")
        scaleOnReferenceLabel = FALSE
        rm(numberOfReferenceCells,referenceCells)
      }
    }
  }

  # preparation of gene information
  genes <- genes[which(genes$gene_biotype %in% c("protein_coding","lncRNA") & genes$chromosome_name %in% c(1:22,"X") & genes$hgnc_symbol !=""),]
  genes$chromosome_num <- genes$chromosome_name
  genes$chromosome_num[which(genes$chromosome_num=="X")]<- 23
  genes$chromosome_num <- as.numeric(genes$chromosome_num)
  genes2 <- unique(genes[,c("hgnc_symbol","chromosome_num","start_position")])

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

  rawCounts <- as.matrix(Seurat::GetAssay(seuratObj, assay = assay)$counts)

  commonGenes <- intersect(rownames(rawCounts),genes2$hgnc_symbol)
  rawCounts <- rawCounts[commonGenes,]
  averageExpression <- rowMeans(rawCounts)
  topExprGenes <- commonGenes[order(averageExpression, decreasing = T)[1:topNGenes]]
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

  genes2 <- genes2[which(genes2$hgnc_symbol %in% topExprGenes),]
  genes2 <- genes2[order(genes2$chromosome_num,genes2$start_position),]

  # Preparation of genomic windows
  genomicWindows <- lapply(c(1:23),function(chrom){
    genesC <- genes2[which(genes2$chromosome_num == chrom),]
    N <- nrow(genesC)
    iter <- round(windowSize/2)
    if(nrow(genesC)>windowSize){
      gw <- lapply(seq(iter+1,N-iter,by=windowStep),function(i)genesC[(i-iter):(i+iter),"hgnc_symbol"])
      names(gw) <- paste0(chrom,".",1:length(gw))
    }else{
      gw <- list(genesC$"hgnc_symbol")
      names(gw) <- paste0(chrom,".",1)
    }
    return(gw)
  })
  genomicWindows <- unlist(genomicWindows,recursive=F)

  genomicScores <- sapply(genomicWindows, function(g) colMeans(normCounts[g,]) )

  if (scaleOnReferenceLabel) {
    if (length(referenceLabel) == 1) {
      Q01Q99 <- apply(genomicScores[referenceCells,],2,stats::quantile,"probs"=c(0+thresholdPercentile,1-thresholdPercentile))
      genomicScoresTrimmed <- apply(genomicScores,1,function(v)
        {v[which(v >= Q01Q99[1,] & v <= Q01Q99[2,])] <- 0; return(v)})
    } else {
      genomicScoresReferenceLabel <- list()
      for (i in referenceLabel){
        genomicScoresReferenceLabel[[i]] <- genomicScores[referenceCells[[i]],]
      }
      genomicScoresReferenceLabel <- do.call(rbind, genomicScoresReferenceLabel)
      Q01Q99 <- apply(genomicScoresReferenceLabel, 2, stats::quantile, "probs"=c(0+thresholdPercentile,1-thresholdPercentile))
      genomicScoresTrimmed <- apply(genomicScores, 1,function(v)
        {v[which(v >= Q01Q99[1,] & v <= Q01Q99[2,])] <- 0 ; v })
    }

  } else {
    Q01Q99 <- apply(genomicScores,2,stats::quantile,"probs"=c(0+thresholdPercentile,1-thresholdPercentile))
    genomicScoresTrimmed <- apply(genomicScores,1,function(v)
    {v[which(v >= Q01Q99[1,] & v <= Q01Q99[2,])] <- 0;v})
  }

  genomicAssay <- Seurat::CreateAssayObject(data = as.matrix(genomicScoresTrimmed))
  seuratObj[["genomicScores"]] <- genomicAssay
  seuratObj[["cnv_fraction"]] <- colMeans(abs(genomicScoresTrimmed) > 0)

  return (seuratObj)
}
