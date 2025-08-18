#' fastCNV calls all of the internal functions needed to compute the putative CNV on a Seurat object or a list of Seurat objects
#'
#' This function orchestrates the CNV analysis on a Seurat object (or multiple objects). It calls internal functions such as
#' `prepareCountsForCNVAnalysis`, `CNVAnalysis`, `CNVPerChromosomeArm`, `CNVcluster`, and `PlotCNVResults` to compute the CNVs,
#' perform clustering, and generate heatmaps. The results are saved in the metadata of the Seurat object(s), with options for
#' generating and saving plots.
#'
#' @param seuratObj Seurat object or list of Seurat objects to perform the CNV analysis on.
#' @param sampleName Name of the sample or a list of names corresponding to the samples in the `seuratObj`.
#' @param referenceVar The variable name of the annotations in the Seurat metadata to be used as reference.
#' @param referenceLabel The label given to the observations you want as reference (can be any type of annotation).
#' @param assay Name of the assay to run the CNV on. Takes the results of `prepareCountsForCNVAnalysis` by default if available.
#' @param prepareCounts If `FALSE`, will not run the `prepareCountsForCNVAnalysis` function (default = `TRUE`).
#' @param aggregFactor The number of counts per spot desired (default = 15 000). If less than 1,000, will not run the `prepareCountsForCNVAnalysis` function.
#' @param seuratClusterResolution The resolution wanted for the Seurat clusters (default = 0.8).
#' @param aggregateByVar If `referenceVar` is given, determines whether to use it to pool the observations (default = `TRUE`).
#' @param reClusterSeurat Whether to re-cluster if the Seurat object given already has a `seurat_clusters` slot in its metadata (default = `FALSE`).
#' @param pooledReference Default is `TRUE`. Will build a pooled reference across all samples if `TRUE`.
#' @param denoise If `TRUE`, the data will be denoised (default = `TRUE`).
#' @param scaleOnReferenceLabel If `TRUE`, scales the results depending on the normal observations (default = `TRUE`).
#' @param thresholdPercentile Which quantiles to take (default 0.01). For example, `0.01` will take quantiles between 0.01-0.99. Background noise appears with higher numbers.
#' @param geneMetadata List of genes and their metadata (default uses genes from Ensembl version 113).
#' @param chrArmsToForce A chromosome arm (e.g., `"8p"`, `"3q"`) or a list of chromosome arms (e.g., `c("3q", "8p", "17p")`) to force into the analysis.
#' @param genesToForce A list of genes to force into the analysis (e.g. `c("FOXP3","MUC16","SAMD15")`).
#' @param regionToForce Chromosome region to force into the analysis (vector containing chr, start, end).
#' @param windowSize Size of the genomic windows for CNV analysis (default = 150).
#' @param windowStep Step between the genomic windows (default = 10).
#' @param saveGenomicWindows If `TRUE`, saves the information of the genomic windows in the current directory (default = `FALSE`).
#' @param topNGenes Number of top expressed genes to keep (default = 7000).
#' @param getCNVPerChromosomeArm If `TRUE`, will save the CNV per chromosome arm into the metadata.
#' @param getCNVClusters If `TRUE`, will perform clustering on the CNV scores and save them in the metadata of the Seurat object as `cnv_clusters`.
#' @param tumorLabel The label within `referenceVar` that specifies the tumor/malignant population (can be any type of annotation).
#' @param k_clusters Optional. Number of clusters to cut the dendrogram into. If `NULL`, the optimal number of clusters is determined automatically using the elbow method.
#' @param h_clusters Optional. The height at which to cut the dendrogram for clustering. If both `k` and `h` are provided, `k` takes precedence.
#' @param plotDendrogram Logical. Whether to plot the dendrogram (default = `FALSE`).
#' @param plotClustersOnDendrogram Logical. Whether to highlight clusters on the dendrogram (default = `FALSE`).
#' @param plotElbowPlot Logical. Whether to plot the elbow plot used for determining the optimal number of clusters (default = `FALSE`).
# @param doRecapPlot Logical. If `TRUE` (default), generates CNV heatmaps grouped by annotation.
#' @param doPlot If `TRUE`, will build a heatmap for each of the samples (default = `TRUE`).
#' @param printPlot If `TRUE`, the heatmap will be printed in the console (default = `FALSE`, the plot will only be saved in a PDF).
#' @param savePath Path to save the heatmap plot. If `NULL`, the plot won't be saved (default = `.`).
#' @param outputType Specifies the file format for saving the plot, either `"png"` or `"pdf"` (default = `"png"`).
#' @param clustersVar The variable name of the clusters in the Seurat metadata (default = `"cnv_clusters"`).
#' @param splitPlotOnVar The name of the metadata column to split the observations during the `plotCNVResults` step, if different from `referenceVar`.
#' @param referencePalette The color palette that should be used for `referenceVar` (default = `"default"`).
#' @param clusters_palette The color palette that should be used for `clustersVar` (default = `"default"`).
#'
#' @return A list of Seurat objects after all the analysis is complete. Heatmaps of the CNVs for every object in `seuratObj` are generated and saved in the specified path (default = current working directory).
#'
#' @importFrom crayon red yellow green black
#'
#' @export

fastCNV <- function (seuratObj,
                     sampleName,
                     referenceVar = NULL,
                     referenceLabel = NULL,
                     assay = NULL,

                     prepareCounts = TRUE,
                     aggregFactor = 15000,
                     seuratClusterResolution = 0.8,
                     aggregateByVar = TRUE,
                     reClusterSeurat = FALSE,

                     pooledReference = TRUE,
                     denoise = TRUE,
                     scaleOnReferenceLabel = TRUE,
                     thresholdPercentile = 0.01,
                     geneMetadata = getGenes(),
                     windowSize = 150,
                     windowStep = 10,
                     saveGenomicWindows = FALSE,
                     topNGenes = 7000,
                     chrArmsToForce = NULL,
                     genesToForce = NULL,
                     regionToForce = NULL,

                     getCNVPerChromosomeArm = TRUE,

                     getCNVClusters = TRUE,
                     tumorLabel = NULL,
                     k_clusters = NULL,
                     h_clusters = NULL,
                     plotDendrogram = FALSE,
                     plotClustersOnDendrogram = FALSE,
                     plotElbowPlot = FALSE,

                     #doRecapPlot = TRUE,
                     doPlot = TRUE,
                     printPlot = FALSE,
                     savePath = ".",
                     outputType = "png",
                     clustersVar = "cnv_clusters",
                     splitPlotOnVar = clustersVar,
                     referencePalette = "default",
                     clusters_palette = "default"
                     ){

  if(!length(seuratObj)==length(sampleName)) stop(crayon::red("seuratObjHD & sampleName should have the same length"))

  options(future.globals.maxSize = 8000*1024^2)

  if(length(seuratObj) == 1){
    seuratObj <- list(seuratObj) ; names(seuratObj) <- sampleName
  }
  if (prepareCounts == TRUE & aggregFactor>=1000) {
    message(crayon::yellow(paste0("[",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"]"," Aggregating counts matrix...")))
    for (i in 1:length(seuratObj)) {
      seuratObj[[i]] <- prepareCountsForCNVAnalysis(seuratObj[[i]], sampleName = sampleName[[i]],
                                                   referenceVar = referenceVar,
                                                   aggregateByVar = aggregateByVar ,
                                                   aggregFactor=aggregFactor,
                                                   seuratClusterResolution = seuratClusterResolution,
                                                   reClusterSeurat = reClusterSeurat  )
      invisible(gc())
    }
    message(crayon::green(paste0("[",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"]"," Done !")))
  } else {
    for (i in 1:length(seuratObj)) {
      seuratObj[[i]]@project.name = sampleName[[i]]
    }
  }


  if (length(seuratObj) > 1) {
    seuratObj <- CNVanalysis(seuratObj, referenceVar = referenceVar, referenceLabel = referenceLabel,
                              #doRecapPlot = doRecapPlot,
                              pooledReference = pooledReference,
                              scaleOnReferenceLabel = scaleOnReferenceLabel, denoise = denoise,
                              assay = assay, thresholdPercentile = thresholdPercentile, geneMetadata = geneMetadata,
                              chrArmsToForce = chrArmsToForce, windowSize = windowSize, windowStep = windowStep,
                              saveGenomicWindows = saveGenomicWindows, topNGenes = topNGenes)
    invisible(gc())
  } else {
    seuratObj <- CNVanalysis(seuratObj[[1]], referenceVar = referenceVar, referenceLabel = referenceLabel,
                         #doRecapPlot = doRecapPlot,
                         pooledReference = pooledReference, denoise = denoise,
                         scaleOnReferenceLabel = scaleOnReferenceLabel, assay = assay,
                         thresholdPercentile = thresholdPercentile, geneMetadata = geneMetadata,
                         chrArmsToForce = chrArmsToForce, windowSize = windowSize, windowStep = windowStep,
                         saveGenomicWindows = saveGenomicWindows, topNGenes = topNGenes)
    invisible(gc())
  }

  if (getCNVPerChromosomeArm == TRUE) {
    message(crayon::yellow(paste0("[",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"]"," Computing CNV per chromosome arm...")))
    if (length(seuratObj) == 1) {
          seuratObj <- CNVPerChromosomeArm(seuratObj)
    } else {
      for (i in 1:length(seuratObj)) {
        seuratObj[[i]] <- CNVPerChromosomeArm(seuratObj[[i]])
      }
    }
    invisible(gc())
    message(crayon::green(paste0("[",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"]"," Done !")))
  }

  if (getCNVClusters == TRUE) {
    message(crayon::yellow(paste0("[",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"]"," Clustering CNVs...")))
    if (length(seuratObj) == 1) {
      seuratObj <- CNVcluster(seuratObj, referenceVar = referenceVar, tumorLabel = tumorLabel, k = k_clusters, h = h_clusters, plotDendrogram = plotDendrogram,
                              plotClustersOnDendrogram = plotClustersOnDendrogram, plotElbowPlot = plotElbowPlot)
    } else {
      for (i in 1:length(seuratObj)) {
        seuratObj[[i]] <- CNVcluster(seuratObj[[i]], referenceVar = referenceVar, tumorLabel = tumorLabel, k = k_clusters, h = h_clusters, plotDendrogram = plotDendrogram,
                                     plotClustersOnDendrogram = plotClustersOnDendrogram, plotElbowPlot = plotElbowPlot)
      }
    }
    invisible(gc())
    message(crayon::green(paste0("[",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"]"," Done !")))
  }

  if (doPlot == TRUE) {
    message(crayon::yellow(paste0("[",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"]"," Plotting CNV heatmap...")))
    if (length(seuratObj) > 1) {
      for (i in 1:length(seuratObj)) {
        if (Seurat::Project(seuratObj[[i]]) == "SeuratProject") {Seurat::Project(seuratObj[[i]]) = paste0("Sample",i)}
        if ("cnv_clusters" %in% names(seuratObj[[i]]@meta.data)) {splitPlotOnVar = "cnv_clusters"}
        plotCNVResults(seuratObj[[i]], referenceVar = referenceVar, splitPlotOnVar = splitPlotOnVar,
                       savePath = savePath, printPlot = printPlot, referencePalette = referencePalette, clusters_palette = clusters_palette, outputType = outputType)
        invisible(gc())
      }
    } else {
      if ("cnv_clusters" %in% names(seuratObj@meta.data)) {splitPlotOnVar = "cnv_clusters"}
      plotCNVResults(seuratObj = seuratObj, referenceVar = referenceVar, splitPlotOnVar = splitPlotOnVar, clustersVar = clustersVar,
                     savePath = savePath, printPlot = printPlot, referencePalette = referencePalette, clusters_palette = clusters_palette, outputType = outputType)
      invisible(gc())
    }
    message(crayon::green(paste0("[",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"]"," Done !")))
  }



  return (seuratObj)
}
