#' CNVclassify
#'
#' @param cnv_vector Vector of CNVs
#' @param peaks_info peaks info
#'
#' @return
#'


classify_cnv <- function(cnv_vector, peaks_info) {
  cnv_classification <- c()  # Initialize a vector to store the classification

  # Extract the peak-related data from the peaks_info list
  big_peaks <- peaks_info$`x abciss big peaks`
  inter_peaks <- peaks_info$`x abciss inter-peaks`
  peaks_size <- peaks_info$`peaks x size`

  # Check the number of peaks
  num_peaks <- peaks_info$`nb big peaks`

  # Iterate through each CNV value and classify it
  for (cnv_value in cnv_vector) {
    classification <- "no_alteration"  # Default classification

    # If there's only 1 peak:
    if (num_peaks == 1) {
      peak <- big_peaks[1]

      if (peak < -0.1) {
        mean_loss_threshold <- mean(c(peak, peaks_size$`right born`[1]))
        if (cnv_value < mean_loss_threshold) {
          classification <- "loss"
        }
      } else if (peak > 0.1) {
        mean_gain_threshold <- mean(c(peak, peaks_size$`left born`[1]))
        if (cnv_value > mean_gain_threshold) {
          classification <- "gain"
        }
      }

      # If there are 2 peaks:
    } else if (num_peaks == 2) {
      peak1 <- big_peaks[1]
      peak2 <- big_peaks[2]

      # Apply conditions for the first peak (loss)
      if (peak1 < -0.1) {
        mean_loss_threshold <- mean(c(peak1, peaks_size$`right born`[1]))
        if (cnv_value < mean_loss_threshold) {
          classification <- "loss"
        }
      }

      # Apply conditions for the second peak (gain)
      if (peak2 > 0.1) {
        mean_gain_threshold <- mean(c(peak2, peaks_size$`left born`[2]))
        if (cnv_value > mean_gain_threshold) {
          classification <- "gain"
        }
      }

      # If there are 3 peaks:
    } else if (num_peaks == 3) {
      peak1 <- big_peaks[1]
      peak2 <- big_peaks[2]
      peak3 <- big_peaks[3]

      # Condition for first peak (loss)
      if (peak1 < -0.1) {
        mean_loss_threshold <- mean(c(peak1, peaks_size$`right born`[1]))
        if (cnv_value < mean_loss_threshold) {
          classification <- "loss"
        }
      }

      # Condition for second peak (no alteration)
      if (peak2 >= -0.1 && peak2 <= 0.1) {
        mean_no_alteration_left <- mean(c(peak2, peaks_size$`left born`[2]))
        mean_no_alteration_right <- mean(c(peak2, peaks_size$`right born`[2]))
        if (cnv_value >= mean_no_alteration_left && cnv_value <= mean_no_alteration_right) {
          classification <- "no_alteration"
        }
      }

      # Condition for third peak (gain)
      if (peak3 > 0.1) {
        mean_gain_threshold <- mean(c(peak3, peaks_size$`left born`[3]))
        if (cnv_value > mean_gain_threshold) {
          classification <- "gain"
        }
      }
    }

    # Append the classification result
    cnv_classification <- c(cnv_classification, classification)
  }

  return(cnv_classification)
}
