#' Checks if working on ubuntu 22
#'
#' This function checks if the user is working on ubuntu 22
#' Raster_by_magick does not work under ubuntu 22 so the parameter on the heatmap is changed if it returns TRUE
#'
#' @return TRUE or FALSE
#'

is_ubuntu22 <- function() {
  if (Sys.info()["sysname"] != "Linux") return(FALSE)
  info <- try(system("lsb_release -d", intern = TRUE), silent = TRUE)
  if (inherits(info, "try-error")) return(FALSE)
  grepl("ubuntu 22", tolower(info))
}

