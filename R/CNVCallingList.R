#' CNVcalling but for a list of seurat objects
#'
#' @param seuratList List of seurat objects
#' @param assay Name of the assay to run the CNV on. Takes the results of `prepareCountsForCNVAnalysis` by default if available
#' @param referenceVar The variable name of the annotations in the Seurat metadata
#' @param referenceLabel  The label given to the observations wanted as reference (can be any type of annotation)
#' @param scaleOnReferenceLabel If you want to scale the results depending on the reference observations
#' @param denoise If the data needs to be denoised (default = `TRUE`)
#' @param thresholdPercentile Which quantiles to take (if 0.01 it will take 0.01-0.99). Background noise appears with higher numbers.
#' @param geneMetadata List of genes from ensembl and their metadata
#' @param windowSize Size of the genomic windows
#' @param windowStep Step between the genomic windows
#' @param saveGenomicWindows If the information of the genomic windows need to be saved in the current directory (default = `FALSE`).
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
                           denoise = TRUE,
                           thresholdPercentile = 0.01,
                           geneMetadata=getGenes(),
                           windowSize=150,
                           windowStep=10,
                           saveGenomicWindows = FALSE,
                           topNGenes=7000
){
  LrawcountsByPatient <- lapply(seuratList, function(x) {
    if (!is.null(assay)) {
      as.matrix(Seurat::GetAssay(x, assay = assay)["counts"])
    } else if ("AggregatedCounts" %in% Seurat::Assays(x)) {
      as.matrix(Seurat::GetAssay(x, assay = "AggregatedCounts")["counts"])
    } else {
      as.matrix(Seurat::GetAssay(x, assay = Seurat::Assays(x)[1])["counts"])
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

  # preparation of gene information
  geneMetadata <- geneMetadata[which(geneMetadata$gene_biotype %in% c("protein_coding","lncRNA") & geneMetadata$chromosome_name %in% c(1:22,"X") & geneMetadata$hgnc_symbol !=""),]
  geneMetadata$chromosome_num <- geneMetadata$chromosome_name
  geneMetadata$chromosome_num[which(geneMetadata$chromosome_num=="X")]<- 23
  geneMetadata$chromosome_num <- as.numeric(geneMetadata$chromosome_num)
  geneMetadata2 <- unique(geneMetadata[,c("hgnc_symbol","chromosome_num","start_position", "chr_arm")])

  # internal functions
  funTrim <- function(normcounts,lo=-3,up=3){
    t(apply(normcounts, 1, function(z) {z[which(z < lo)]<-lo ; z[which(z > up)]<-up ; z}))
  }

  funGenomicScore <- function(normcounts,GW=genomicWindows){
    res <- sapply(GW, function(g) colMeans(normcounts[g,]))
    res
  }

  # Analyses

  commonGenes <- Reduce(intersect, c(lapply(LrawcountsByPatient, rownames), list(geneMetadata2$hgnc_symbol)))
  LrawcountsByPatient <- lapply(LrawcountsByPatient, function(x) x[commonGenes,])

  aveExpr <- compute_average_expression(LN, LrawcountsByPatient)
  aveExpr <- do.call(cbind,aveExpr)
  aveExpr <- rowMeans(aveExpr, na.rm = T)
  if (length(aveExpr) < topNGenes) {topNGenes = length(aveExpr)}
  topExprGenes <- commonGenes[order(aveExpr, decreasing = T)[1:topNGenes]]

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
      top_arm_genes <- remaining_genes[order(aveExpr[commonGenes %in% remaining_genes], decreasing = TRUE)[1:200]]
      genes_by_arm[[arm]] <- unique(c(genes_by_arm[[arm]], top_arm_genes))
    }
  }
  final_selected_genes <- unlist(genes_by_arm)

  LrawcountsByPatient <- lapply(LrawcountsByPatient, function(x) x[final_selected_genes,])

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

  LgenomicScores <- lapply(LnormcountsByPatient,funGenomicScore,"GW"=genomicWindows)

  if (denoise) {
    if(scaleOnReferenceLabel){
      if (length(referenceLabel) == 1){
        genomicScoresReferenceLabel <- do.call(rbind, lapply(names(LN), function(patient) LgenomicScores[[patient]][LN[[patient]],]))
      } else {
        genomicScoresReferenceLabel <- list()
        for (i in referenceLabel){
          genomicScoresReferenceLabel[[i]] <- do.call(rbind, lapply(names(LN[[i]]), function(patient) LgenomicScores[[patient]][LN[[i]][[patient]],]))
        }
        genomicScoresReferenceLabel <- do.call(rbind, genomicScoresReferenceLabel)
      }
      Q01Q99 <- apply(genomicScoresReferenceLabel, 2, stats::quantile, "probs"=c(0+thresholdPercentile,1-thresholdPercentile))
      LgenomicScoresTrimmed <- lapply(LgenomicScores, function(gs){
        apply(gs, 1, function(v) { v[which(v >= Q01Q99[1,] & v <= Q01Q99[2,])] <- 0 ; v }) })

    } else {
      genomicScoresAll <- do.call(rbind, lapply(names(LgenomicScores), function(patient) LgenomicScores[[patient]]))
      Q01Q99 <- apply(genomicScoresAll,2,quantile,"probs"=c(0+thresholdPercentile,1-thresholdPercentile))
      LgenomicScoresTrimmed <- lapply(LgenomicScores,function(gs){
        apply(gs,1,function(v) { v[which(v >= Q01Q99[1,] & v <= Q01Q99[2,])] <- 0 ; v }) })
    }
  } else {
    LgenomicScoresTrimmed <- lapply(LgenomicScores, function(x) t(x))
  }

  if (saveGenomicWindows){
    save(genomicWindows, file = paste0("genomicWindows_size",windowSize,"_step",windowStep,".RData"))
  }

  for (i in 1:length(seuratList)) {
    genomicAssay <- Seurat::CreateAssayObject(counts = LgenomicScoresTrimmed[i][[1]])
    seuratList[[i]][["genomicScores"]] <- genomicAssay
    seuratList[[i]][["cnv_fraction"]] <- colMeans(abs(LgenomicScoresTrimmed[i][[1]]) > 0)
  }

  return (seuratList)
}
