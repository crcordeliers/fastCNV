#' fastCNV calls all of the internal functions needed to compute the putative CNV on a seurat object or a list of seurat objects
#'
#' @param seuratObj Seurat object or list of seurat objects
#' @param sampleName Name of the sample or list of names corresponding to the samples in the `seuratObj`
#' @param referenceVar The variable name of the annotations in the Seurat metadata
#' @param referenceLabel The label given to the observations wanted as reference (can be any type of annotation)
#' @param clustersVar The variable name of the clusters in the Seurat object metadata
#' @param splitPlotOnVar The name given to the meta.data slot to split the observations on during the `plotCNVResults` step, if different from `referenceVar`.
#' @param doPlot If `TRUE` will build a heatmap for each of the samples (default = `TRUE`)
#' @param printPlot If the heatmap should be printed in the console (default = `FALSE`, the plot will only be in a PDF).
#' @param prepareCounts If `FALSE` will not run the `prepareCountsForCNVAnalysis` function. Default `TRUE`.
#' @param pooledReference Default `TRUE`. Will build a pooled reference across all samples if `TRUE`.
#' @param assay Name of the assay to run the CNV on. Takes the results of `prepareCountsForCNVAnalysis` by default if available
#' @param getCNVPerChromosomeArm If `TRUE`, will save into the metadata the CNV per chromosome arm
#' @param getCNVClusters If `TRUE`, will perform clustering on the CNV scores, and save them in the metadata of the seurat object as `cnv_clusters`.
#' @param savePath Default `.`. If `NULL` the heatmap won't be saved as a pdf.
#' @param aggregFactor The number of counts per spot desired (default `30 000`). If less than `1000`, will not run the `prepareCountsForCNVAnalysis` function.
#' @param seuratClusterResolution The resolution wanted for the seurat clusters (default 0.8)
#' @param aggregateByVar If `referenceVar` is given, whether to use it to pool the observations
#' @param reClusterSeurat Whether to re-cluster if the Seurat object given already has a `seurat_clusters` slot in its meta.data
#' @param denoise If the data needs to be denoised (default = `TRUE`)
#' @param referencePalette The color palette that should be used for `referenceVar`
#' @param clusters_palette The color palette that should be used for `clustersVar`
#' @param scaleOnReferenceLabel If you want to scale the results depending on the normal observations
#' @param thresholdPercentile Which quantiles to take (if 0.01 it will take 0.01-0.99). Background noise appears with higher numbers.
#' @param geneMetadata List of genes and their metadata (default genes from ensembl version 113)
#' @param windowSize Size of the genomic windows
#' @param windowStep Step between the genomic windows
#' @param saveGenomicWindows If the information of the genomic windows need to be saved in the current directory (default = `FALSE`).
#' @param topNGenes Number of top expressed genes to keep
#'
#' @return This function returns a list of the seurat objects after all the analysis and it creates the heatmaps of the CNVs of every object in seuratObj and saves them in a pdf in the current working directory
#'
#' @export
#'


fastCNV <- function (seuratObj,
                     sampleName,
                     referenceVar = NULL,
                     referenceLabel = NULL,
                     clustersVar = "cnv_clusters",
                     splitPlotOnVar = clustersVar,
                     doPlot = TRUE,
                     printPlot = FALSE,
                     prepareCounts = TRUE,
                     pooledReference = TRUE,
                     assay = NULL,
                     getCNVPerChromosomeArm = TRUE,
                     getCNVClusters = TRUE,
                     savePath = ".",
                     denoise = TRUE,
                     referencePalette = as.character(paletteer::paletteer_d("pals::glasbey")),
                     clusters_palette = scales::hue_pal(),
                     aggregFactor=15000,
                     seuratClusterResolution = 0.8,
                     aggregateByVar = T,
                     reClusterSeurat = F,
                     scaleOnReferenceLabel = TRUE,
                     thresholdPercentile = 0.01,
                     geneMetadata=getGenes(),
                     windowSize=150,
                     windowStep=10,
                     saveGenomicWindows = FALSE,
                     topNGenes=7000){

  if (getCNVClusters == FALSE) {
    if (!is.null(clustersVar) && clustersVar == "cnv_clusters") {
      if(splitPlotOnVar == clustersVar) {
        splitPlotOnVar = referenceVar
      }
      clustersVar = NULL
    }
  }

  if(!length(seuratObj)==length(sampleName)) stop("error-fastCNV : seuratObj & sampleName should have the same length")

  options(future.globals.maxSize = 8000*1024^2)

  if(length(seuratObj) == 1){
    seuratObj <- list(seuratObj) ; names(seuratObj) <- sampleName
  }
  if (prepareCounts == TRUE & aggregFactor>=1000) {
    message("Aggregating counts matrix.")
    for (i in 1:length(seuratObj)) {
      seuratObj[[i]] <- prepareCountsForCNVAnalysis(seuratObj[[i]], sampleName = sampleName[[i]],
                                                   referenceVar = referenceVar,
                                                   aggregateByVar = aggregateByVar ,
                                                   aggregFactor=aggregFactor,
                                                   seuratClusterResolution = seuratClusterResolution,
                                                   reClusterSeurat = reClusterSeurat  )
      invisible(gc())
    }
  } else {
    for (i in 1:length(seuratObj)) {
      seuratObj[[i]]@project.name = sampleName[[i]]
    }
  }


  if (length(seuratObj) > 1) {
    message("Running CNVAnalysis")
    seuratObj <- CNVanalysis(seuratObj, referenceVar = referenceVar, referenceLabel = referenceLabel,
                              #doRecapPlot = doRecapPlot,
                              pooledReference = pooledReference,
                              scaleOnReferenceLabel = scaleOnReferenceLabel, denoise = denoise,
                              assay = assay, thresholdPercentile = thresholdPercentile, geneMetadata = geneMetadata,
                              windowSize = windowSize, windowStep = windowStep, saveGenomicWindows = saveGenomicWindows, topNGenes = topNGenes)
    invisible(gc())
    message("CNVAnalysis done!")
  } else {
    message("Running CNVAnalysis")
    seuratObj <- CNVanalysis(seuratObj[[1]], referenceVar = referenceVar, referenceLabel = referenceLabel,
                         #doRecapPlot = doRecapPlot,
                         pooledReference = pooledReference, denoise = denoise,
                         scaleOnReferenceLabel = scaleOnReferenceLabel, assay = assay,
                         thresholdPercentile = thresholdPercentile, geneMetadata = geneMetadata,
                         windowSize = windowSize, windowStep = windowStep, saveGenomicWindows = saveGenomicWindows, topNGenes = topNGenes)
    invisible(gc())
    message("CNVAnalysis done!")
  }

  if (getCNVPerChromosomeArm == TRUE) {
    if (length(seuratObj) == 1) {
          seuratObj <- CNVPerChromosomeArm(seuratObj)
    } else {
      for (i in 1:length(seuratObj)) {
        seuratObj[[i]] <- CNVPerChromosomeArm(seuratObj[[i]])
      }
    }
    invisible(gc())
  }

  if (getCNVClusters == TRUE) {
    if (length(seuratObj) == 1) {
      seuratObj <- CNVcluster(seuratObj)
    } else {
      for (i in 1:length(seuratObj)) {
        seuratObj[[i]] <- CNVcluster(seuratObj[[i]])
      }
    }
    invisible(gc())
  }

  if (doPlot == TRUE) {
    message("Plotting CNV results. This step may take some time.")
    if (length(seuratObj) > 1) {
      for (i in 1:length(seuratObj)) {
        if (Seurat::Project(seuratObj[[i]]) == "SeuratProject") {Seurat::Project(seuratObj[[i]]) = paste0("Sample",i)}
        if ("cnv_clusters" %in% names(seuratObj[[i]]@meta.data)) {splitPlotOnVar = "cnv_clusters"}
        plotCNVResults(seuratObj[[i]], referenceVar = referenceVar, splitPlotOnVar = splitPlotOnVar,
                       savePath = savePath, printPlot = printPlot, referencePalette = referencePalette, clusters_palette = clusters_palette)
        invisible(gc())
      }
    } else {
      if ("cnv_clusters" %in% names(seuratObj@meta.data)) {splitPlotOnVar = "cnv_clusters"}
      plotCNVResults(seuratObj = seuratObj, referenceVar = referenceVar, splitPlotOnVar = splitPlotOnVar, clustersVar = clustersVar,
                     savePath = savePath, printPlot = printPlot, referencePalette = referencePalette, clusters_palette = clusters_palette)
      invisible(gc())
    }
  }



  return (seuratObj)
}
