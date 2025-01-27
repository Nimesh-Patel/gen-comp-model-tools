#' Validate intervention start and end times
#'
#' Helper function to validate intervention start and end times, used as part
#' of a validation function
#' @param param description
#' @return stops with information or indicates checks passed
#' }
check_intervention_times <- function(intervention_start_time,
                                     intervention_end_time) {
  if (!is.null(
    intervention_start_time
  ) &&
    !(
      is.numeric(intervention_start_time) &&
        length(intervention_start_time) == 1
    )) {
    stop("intervention_start_time must be a single numeric value or NULL.")
  }

  if (!is.null(
    intervention_end_time
  ) &&
    !(
      is.numeric(intervention_end_time) &&
        length(intervention_end_time) == 1
    )) {
    stop("intervention_end_time must be a single numeric value or NULL.")
  }

  print("intervention start and end time checks passed")
}
