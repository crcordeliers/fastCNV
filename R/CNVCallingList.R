#' CNVcalling but for a list of seurat objects
#'
#' @param seuratList List of seurat objects
#' @param assay Name of the assay to run the CNV on. Takes the results of `prepareCountsForCNVAnalysis` by default if available
#' @param referenceVar The variable name of the annotations in the Seurat metadata
#' @param referenceLabel  The label given to the observations wanted as reference (can be any type of annotation)
#' @param scaleOnReferenceLabel If you want to scale the results depending on the reference observations
#' @param thresholdPercentile Which quantiles to take (if 0.01 it will take 0.01-0.99). Background noise appears with higher numbers.
#' @param genes List of genes from ensembl
#' @param windowSize Size of the genomic windows
#' @param windowStep Step between the genomic windows
#' @param topNGenes Number of top expressed genes
#'
#' @return This function returns a list of the genomic scores per genomic window
#'
#' @import stats
#' @import Seurat
#'
#' @export
#'


CNVcallingList <- function(seuratList,
                           assay = NULL,
                           referenceVar = NULL,
                           referenceLabel = NULL,
                           scaleOnReferenceLabel= TRUE,
                           thresholdPercentile = 0.01,
                           genes=getGenes(),
                           windowSize=100,
                           windowStep=20,
                           topNGenes=7000
){
  LrawcountsByPatient <- lapply(seuratList, function(x) {
    if (!is.null(assay)) {
      as.matrix(Seurat::GetAssay(x, assay = assay)$counts)
    } else if ("AggregatedCounts" %in% Seurat::Assays(x)) {
      as.matrix(Seurat::GetAssay(x, assay = "AggregatedCounts")$counts)
    } else {
      as.matrix(Seurat::GetAssay(x, assay = Seurat::Assays(x)[1])$counts)
    } } )

  names(LrawcountsByPatient) <- lapply(seuratList, function(x) Seurat::Project(x))
  Lannot <- lapply(seuratList, function(x) Seurat::FetchData(x, vars = referenceVar))
  names(Lannot) <- lapply(seuratList, function(x) Seurat::Project(x))


  # getting reference items per patient
  if (is.null(referenceVar) || is.null(referenceLabel)){
    # unable to scale the results on reference data if we don't know what the reference data is
    warning("referenceVar and/or referenceLabel parameters not found. Computing the CNV without a reference.")
    scaleOnReferenceLabel = FALSE
  } else if ( as.numeric(sum(sapply(lapply(Lannot, function(annot) rownames(annot)[which(annot == referenceLabel)]), function(x) length(x)))) == 0 ) {
    warning(paste0("No observations found for annotation labels ", referenceLabel,". Computing the CNV without a reference."))
    scaleOnReferenceLabel = FALSE
  } else if (length(referenceLabel) == 1){
    LN <- lapply(Lannot, function(annot) rownames(annot)[which(annot == referenceLabel)])
    names(LN) <- names(Lannot)
    LN <- LN[which(sapply(LN,length) >= 5)]
    LN <- Filter(function(x) length(x) > 0, LN)
  } else {
    LN <- list()
    for (i in referenceLabel){
      LN[[i]] <- lapply(Lannot, function(annot) rownames(annot)[which(annot == i)])
      LN[[i]] <- LN[[i]][which(sapply(LN[[i]],length) >= 5)]
      LN <- Filter(function(x) length(x) > 0, LN)
    }
  }


  # preparation of gene info
  genes <- genes[which(genes$gene_biotype %in% c("protein_coding","lncRNA") & genes$chromosome_name %in% c(1:22,"X") & genes$hgnc_symbol !=""),]
  genes$chromosome_num <- genes$chromosome_name
  genes$chromosome_num[which(genes$chromosome_num=="X")]<- 23
  genes$chromosome_num <- as.numeric(genes$chromosome_num)

  genes2 <- unique(genes[,c("hgnc_symbol","chromosome_num","start_position")])


  # internal functions

  funTrim <- function(normcounts,lo=-3,up=3){
    t(apply(normcounts, 1, function(z) {z[which(z < lo)]<-lo ; z[which(z > up)]<-up ; z}))
  }

  funGenomicScore <- function(normcounts,GW=genomicWindows){
    res <- sapply(GW, function(g) colMeans(normcounts[g,]))
    res
  }

  # Analyses

  commonGenes <- Reduce(intersect, c(lapply(LrawcountsByPatient, rownames), list(genes2$hgnc_symbol)))
  LrawcountsByPatient <- lapply(LrawcountsByPatient, function(x) x[commonGenes,])

  aveExpr <- rowMeans(sapply(LrawcountsByPatient, rowMeans))
  topExprGenes <- commonGenes[order(aveExpr, decreasing = T)[1:topNGenes]]
  LrawcountsByPatient <- lapply(LrawcountsByPatient, function(x) x[topExprGenes,])

  LnormcountsByPatient <- lapply(LrawcountsByPatient, function(d) log2(1+d))



  LnormcountsByPatient <- lapply(LnormcountsByPatient, scale, "scale"=F)


  if(scaleOnReferenceLabel){
    if (length(referenceLabel) == 1) {
      SF <- rowMeans(sapply(names(LN), function(patient) rowMeans(LnormcountsByPatient[[patient]][,LN[[patient]]])))
    } else {
      SF <- list()
      for (i in referenceLabel) {
        SF[[i]] <- rowMeans(sapply(names(LN[[i]]), function(patient) rowMeans(LnormcountsByPatient[[patient]][,LN[[i]][[patient]]])))
      }
      SF <- do.call(rbind,SF)
      SF <- apply(SF, 2, median)
    }
  }else{
    SF <- rowMeans(sapply(LnormcountsByPatient,rowMeans))
  }


  LnormcountsByPatient <- lapply(LnormcountsByPatient, function(d) d-SF)

  LnormcountsByPatient <- lapply(LnormcountsByPatient, function(d) funTrim(d, lo=-3, up=3))

  genes2 <- genes2[which(genes2$hgnc_symbol %in% topExprGenes),]
  genes2 <- genes2[order(genes2$chromosome_num,genes2$start_position),]

  # preparation of genomic windows
  genomicWindows <- lapply(c(1:23), function(chrom){
    genesC <- genes2[which(genes2$chromosome_num == chrom),]
    N <- nrow(genesC)
    iter <- round(windowSize/2)
    if(nrow(genesC)>windowSize){
      gw <- lapply(seq(iter+1,N-iter,by=windowStep), function(i) genesC[(i-iter):(i+iter),"hgnc_symbol"])
      names(gw) <- paste0(chrom, ".", 1:length(gw))
    }else{
      gw <- list(genesC$"hgnc_symbol")
      names(gw) <- paste0(chrom,".",1)
    }
    gw
  })
  genomicWindows <- unlist(genomicWindows,recursive=F)

  LgenomicScores <- lapply(LnormcountsByPatient,funGenomicScore,"GW"=genomicWindows)

  if(scaleOnReferenceLabel){
    if (length(referenceLabel) == 1){
      genomicScoresReferenceLabel <- do.call(rbind, lapply(names(LN), function(patient) LgenomicScores[[patient]][LN[[patient]],]))
      Q01Q99 <- apply(genomicScoresReferenceLabel, 2, stats::quantile, "probs"=c(0+thresholdPercentile,1-thresholdPercentile))
      LgenomicScoresTrimmed <- lapply(LgenomicScores, function(gs){
          apply(gs, 1, function(v) { v[which(v >= Q01Q99[1,] & v <= Q01Q99[2,])] <- 0 ; v }) })
    } else {
      genomicScoresReferenceLabel <- list()
      for (i in referenceLabel){
        genomicScoresReferenceLabel[[i]] <- do.call(rbind, lapply(names(LN[[i]]), function(patient) LgenomicScores[[patient]][LN[[i]][[patient]],]))
      }
      genomicScoresReferenceLabel <- do.call(rbind, genomicScoresReferenceLabel)
      Q01Q99 <- apply(genomicScoresReferenceLabel, 2, stats::quantile, "probs"=c(0+thresholdPercentile,1-thresholdPercentile))
      LgenomicScoresTrimmed <- lapply(LgenomicScores, function(gs){
        apply(gs, 1, function(v) { v[which(v >= Q01Q99[1,] & v <= Q01Q99[2,])] <- 0 ; v }) })
    }
  } else {
    genomicScoresAll <- do.call(rbind, lapply(names(LgenomicScores), function(patient) LgenomicScores[[patient]]))
    Q01Q99 <- apply(genomicScoresAll,2,quantile,"probs"=c(0+thresholdPercentile,1-thresholdPercentile))
    LgenomicScoresTrimmed <- lapply(LgenomicScores,function(gs){
      apply(gs,1,function(v) { v[which(v >= Q01Q99[1,] & v <= Q01Q99[2,])] <- 0 ; v }) })
    }

  for (i in 1:length(seuratList)) {
      genomicAssay <- Seurat::CreateAssayObject(counts = LgenomicScoresTrimmed[i][[1]])
      seuratList[[i]][["genomicScores"]] <- genomicAssay
      seuratList[[i]][["cnv_fraction"]] <- colMeans(abs(LgenomicScoresTrimmed[i][[1]]) > 0)
  }
  return (seuratList)
}
