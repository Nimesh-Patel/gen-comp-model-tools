#' Aggregate incidence at reporting times
#'
#' @param t_incident numeric vector of times when incidence changes
#' @param t_report numeric vector of reporting times to get aggregate incidence
#' @param cum_incidence numeric vector of cumulative incidence
#' matching t_incident
#' defaults to NA which treats t_incident as times of single incidence cases
#' @param t0 optional numeric for initial reporting time
#' default is 0 meaning the first box contains all incident times between 0 and
#' the first element of t_report
#' @param inc0 optional numeric for cumulative incidence reporting time
#' default is 0 meaning the first box contains all incident times between 0 and
#' the first element of t_report
#' @return numeric vector of new incidence at reporting times
#' @family calculations
#' @export
report_incidence <- function(t_incident,
                             t_report,
                             cum_incidence = NA,
                             t0 = 0,
                             inc0 = 0) {
  if (identical(cum_incidence, NA)) {
    message("No input cumulative incidence, treating input t_incident times as
      single incident cases.")

    cum_incidence <- seq_along(t_incident)
  }
  if (!identical(t_report, sort(t_report))) {
    stop("Report times (t_report) must be ordered. Sort input first.")
  }

  if (!identical(t_incident, sort(t_incident))) {
    stop("Incident times (t_incident) must be ordered.
      Check how it was generated.")
  }

  if (!identical(cum_incidence, sort(cum_incidence))) {
    stop("Input incidence (cum_incidence) must be cumulative i.e., must
    monotonically increase. Check how it is generated.")
  }

  currtbefore <- t0
  incbefore <- inc0
  newinc <- incbefore

  cnt <- 0

  reported_incidence <- numeric(length(t_report))
  for (currtend in t_report) {
    cnt <- cnt + 1
    currlogic <- (t_incident <= currtend)
    cuminc_allbefore <- cum_incidence[currlogic]
    if (length(cuminc_allbefore) > 0) {
      newinc <- cuminc_allbefore[length(cuminc_allbefore)]
    }
    reported_incidence[[cnt]] <- newinc - incbefore
    currtbefore <- currtend
    incbefore <- newinc
  }

  reported_incidence
}
