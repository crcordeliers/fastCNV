#' CNVclassify
#'
#' @param cnv_vector Vector of CNVs
#' @param peaks A list containing the thresholds for loss or gain.
#'
#' @return The classification of each spot for a given chromosome arm
#'
#' @keywords internal
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
