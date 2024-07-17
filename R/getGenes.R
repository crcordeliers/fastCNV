#' Downloads the gene information from the ensembl site
#'
#' @param filters Filters (one or more) that should be used in the query
#' @param cache Whether to use the data from the cache or to download the latest version
#'
#' @return The function returns a list of all the genes and the information related to these genes
#' @import biomaRt
#' @importFrom utils data
#' @export

getGenes <- function(filters=NULL, cache=TRUE){
  if(cache) {
    message("Using genes data from ensembl version 111.")
    data(genes,envir = environment())
  } else {
    message("Retrieving data from ensembl current version...")
    message("This could take some time.")
    ensembl = biomaRt::useMart("ensembl")
    ensembl = biomaRt::useDataset("hsapiens_gene_ensembl",mart=ensembl)
    # attributes = listAttributes(ensembl)
    columnstokeep <- c('ensembl_gene_id',   'hgnc_symbol',   'entrezgene_id',"chromosome_name" ,'start_position','end_position','gene_biotype')
    #   filters = listFilters(ensembl)
    #   'with_entrezgene' , 'with_hgnc',  'hgnc_symbol' ...
    #
    if(is.null(filters)){
      genes <- biomaRt::getBM(attributes=columnstokeep,values=T,mart = ensembl)
    }else{
      genes <- biomaRt::getBM(attributes=columnstokeep,filters = filters,values=T,mart = ensembl)
    }
    genes$length <- 1+genes$end_position - genes$start_position
    genes$length.kb=genes$length/1000
  }

  return(genes)
}
