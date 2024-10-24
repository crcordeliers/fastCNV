#' CNVclassify
#'
#' @param cnv_vector Vector of CNVs
#' @param peaks_info peaks info
#'
#' @return
#'


classify_cnv <- function(cnv_vector, peaks) {
  cnv_classification <- c()  # Initialize a vector to store the classification

  for (cnv in cnv_vector){
    if (cnv < peaks[1]) {
      classification = "loss"
    } else if (cnv > peaks[3]) {
      classification = "gain"
    } else {
      classification = "no_alteration"
    }
    cnv_classification <- c(cnv_classification, classification)
  }

  return(cnv_classification)
}
