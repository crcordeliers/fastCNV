#' Genes data from ensemble version 111
#'
#' Data downloaded from the ensembl website version 113
#' Data containing ensembl gene id, HUGO nomenclature and entrez gene id,
#' location on chromosome, biotype and length for each gene (about 76 000 genes)
#'
#' @docType data
#'
#' @usage data(geneMetadata)
#'
#' @format An object of class list.
#'
#' @keywords datasets
#'
#' @source https://www.ensembl.org/index.html
#'
#' @examples
#' data(geneMetadata)
#' hgnc <- geneMetadata$hgnc_symbol
#' entrez <- geneMetadata$entrezgene_id
"geneMetadata"
