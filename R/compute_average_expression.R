#' Compute average expression for patients
#'
#' This function calculates the average gene expression for each patient across
#' different cell types. It first retrieves patient data from LN, then extracts
#' the corresponding count data from LrawcountsByPatient, and calculates the mean expression.
#'
#' @param LN A list where each element represents a cell type with sublists containing patient data.
#' @param LrawcountsByPatient A named list where each element contains count data for a specific patient.
#'
#' @return A named vector containing the average expression for each patient.
#'
compute_average_expression <- function(LN, LrawcountsByPatient) {

  patients = names(LrawcountsByPatient)
  results <- list()

  get_patient_data <- function(patient_name, LN) {
    for (type_cellulaire in names(LN)) {
      if (patient_name %in% names(LN[[type_cellulaire]])) {
        return(LN[[type_cellulaire]][[patient_name]])
      }
    }
    return(NULL)
  }

  for (patient_name in patients) {
    patient_data <- get_patient_data(patient_name, LN)

    if (!is.null(patient_data)) {
      patient_counts <- LrawcountsByPatient[[patient_name]]
      avg_expr <- rowMeans(patient_counts, na.rm = TRUE)
      results[[patient_name]] <- avg_expr
    } else { results[[patient_name]] <- NA }
  }
  return(results)
}

