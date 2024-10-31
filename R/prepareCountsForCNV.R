#' A function to aggregate observations with the same cell types to get higher counts per observation.
#' It will help to have higher counts to compute the CNV
#'
#' @param seuratObj Seurat object containing the data
#' @param sampleName Name of the sample
#' @param referenceVar  The variable name of the annotations in the Seurat metadata
#' @param aggregateByVar If a referenceVar is given, whether to use it to aggregate the observations depending on their type. Default `TRUE`.
#' @param aggregFactor The number of counts per observation desired (default 30 000)
#' @param seuratClusterResolution The resolution wanted for the seurat clusters (default 0.8)
#' @param reClusterSeurat If `TRUE`, clustering is re-done for the Seurat object.
#'
#' @return This function returns a Seurat object with the modified count matrix in a new assay called `AgregatedCounts` and the Seurat clusters in the metadata
#'
#' @import Seurat
#' @import hdf5r
#'
#' @export
#'
prepareCountsForCNVAnalysis <- function(seuratObj, sampleName = NULL, referenceVar = NULL,
                                        aggregateByVar = T, aggregFactor=15000, seuratClusterResolution = 0.8,
                                        reClusterSeurat = F ){

    assay <- Seurat::Assays(seuratObj)[1]
      if (is.null(referenceVar) || aggregateByVar == F) {
        if (!"seurat_clusters" %in% colnames(seuratObj[[]]) | reClusterSeurat) {
          print("Running Seurat SCTransform and clustering. This could take some time.")
          if (!is.null(sampleName)){print(paste0("Sample : ", sampleName))}
          seuratObj <- Seurat::SCTransform(seuratObj, assay = assay)
          seuratObj <- Seurat::RunPCA(seuratObj, assay = "SCT")
          seuratObj <- Seurat::FindNeighbors(seuratObj, reduction = "pca", dims = 1:10)
          seuratObj <- Seurat::FindClusters(seuratObj, resolution = seuratClusterResolution)
          if (!is.null(sampleName)){print(paste0("Seurat SCTransform and clustering done for sample ", sampleName))}
          else {print("Seurat SCTransform and clustering done.")}
        }

        countsMat <- as.matrix(Seurat::GetAssay(seuratObj, assay = assay)["counts"])
        LC = split(Seurat::Cells(seuratObj), Seurat::FetchData(seuratObj, vars = "seurat_clusters"))
        LC = lapply(lapply(LC, function(cells) {nc <- Seurat::FetchData(seuratObj, vars = paste0("nCount_",assay))[cells,] ; names(nc) <- cells ; nc}), sort)

        seuratObj <- Seurat::AddMetaData(seuratObj, metadata = 0, col.name = "metaSpots")

        nbMetaSpots <- 1

        for (X in LC) {
          if (length(X) != 0) {
            barcodes = c()
            taille = 0

            for(i in 1:length(X)) {
              taille = taille + as.numeric(X[i])
              barcodes = append(barcodes, names(X[i]))

              if (taille >= aggregFactor) {
                if (length(barcodes)>1) {
                  countsMat[,barcodes] <- rowSums(countsMat[,barcodes])
                }
                seuratObj@meta.data[barcodes,"metaSpots"] <- nbMetaSpots
                nbMetaSpots = nbMetaSpots + 1
                taille = 0
                barcodes = c()
              }
              if (i == length(X)){
                if (length(barcodes)>1) {
                  countsMat[,barcodes] <- rowSums(countsMat[,barcodes])
                }
                seuratObj@meta.data[barcodes,"metaSpots"] <- nbMetaSpots
                nbMetaSpots = nbMetaSpots + 1
              }
            }
          }
        }
        aggregAssay <- Seurat::CreateAssayObject(counts = countsMat)
        seuratObj[["AggregatedCounts"]] <- aggregAssay

      } else {

        if (!"seurat_clusters" %in% colnames(seuratObj[[]]) | reClusterSeurat) {
          print("Running Seurat SCTransform and clustering. This could take some time.")
          if (!is.null(sampleName)){print(paste0("Sample : ", sampleName))}
          seuratObj <- Seurat::SCTransform(seuratObj, assay = assay)
          seuratObj <- Seurat::RunPCA(seuratObj, assay = "SCT")
          seuratObj <- Seurat::FindNeighbors(seuratObj, reduction = "pca", dims = 1:10)
          seuratObj <- Seurat::FindClusters(seuratObj, resolution = seuratClusterResolution)
          if (!is.null(sampleName)){print(paste0("Seurat SCTransform and clustering done for sample ", sampleName,"."))}
          else {print("Seurat SCTransform and clustering done.")}
        }

        countsMat <- as.matrix(Seurat::GetAssay(seuratObj, assay = Seurat::Assays(seuratObj)[1])["counts"])
        spotsCateg <- split(Seurat::Cells(seuratObj),Seurat::FetchData(seuratObj, vars = referenceVar))
        LC <- lapply(spotsCateg, function(x) split(x, Seurat::FetchData(seuratObj, vars = "seurat_clusters", cells = x)))

        remove_empty_lists <- function(list_of_lists) {
          lapply(list_of_lists, function(sub_list) {
            Filter(function(x) length(x) > 0, sub_list)
          })
        }

        LC <- remove_empty_lists(LC)
        LC = lapply(LC, function(x) {lapply(
                lapply(x, function(cells) {
                  nc <- Seurat::FetchData(seuratObj, vars = paste0("nCount_",assay))[cells,]
                  names(nc) <- cells
                  nc}
                ), sort)})

        seuratObj <- Seurat::AddMetaData(seuratObj, metadata = 0, col.name = "metaSpots")

        nbMetaSpots <- 1

        for (i in LC) {
          for (X in i) {
            if (length(X) != 0) {
              barcodes = c()
              taille = 0

              for(i in 1:length(X)) {
                taille = taille + as.numeric(X[i])
                barcodes = append(barcodes, names(X[i]))

                if (taille >= aggregFactor) {
                  if (length(barcodes)>1) {
                    countsMat[,barcodes] <- rowSums(countsMat[,barcodes])
                  }
                  seuratObj@meta.data[barcodes,"metaSpots"] <- nbMetaSpots
                  nbMetaSpots = nbMetaSpots + 1
                  taille = 0
                  barcodes = c()
                }
                if (i == length(X)){
                  if (length(barcodes)>1) {
                    countsMat[,barcodes] <- rowSums(countsMat[,barcodes])
                  }
                  seuratObj@meta.data[barcodes,"metaSpots"] <- nbMetaSpots
                  nbMetaSpots = nbMetaSpots + 1
                }
              }
            }
          }
        }
        aggregAssay <- Seurat::CreateAssayObject(counts = countsMat)
        seuratObj[["AggregatedCounts"]] <- aggregAssay

      }

  if (!is.null(sampleName)) {
    if (Seurat::Project(seuratObj) == "SeuratProject") {
       Seurat::Project(seuratObj) <- sampleName
    }}

  return(seuratObj)
}
