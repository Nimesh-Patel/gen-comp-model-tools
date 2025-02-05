#' add travel between metapopulations in model
#'
#' @param peterlist list of instructions for piping |>
#' @param rates character/numeric vector movement rate/s
#' @param metapop_from character vector of origin
#' metapopulations
#' default is "" which specifies all metapopulations
#' @param metapop_to character vector of destination
#' metapopulations. default is "" which specifies all
#' metapopulations
#' @param travelstates character vector of basestates that
#' trave between metapopulations
#' default is "" which specifies all basestates
#' @return updated instruction list
#' @export
add_travel <- function(peterlist,
                       rates,
                       metapop_from = "",
                       metapop_to = "",
                       travelstates = "") {
  numtravel <- length(metapop_from)
  if (length(rates) == 1) {
    rates <- rep(rates, numtravel)
  }
  # bad, clean this for edge cases
  if (is.vector(travelstates)) {
    travelstates <- rep(list(travelstates), numtravel)
  }
  if ((is.list(travelstates)) && (length(travelstates) != numtravel)) {
    stop("List of travel states is
    not same length as input metapopulations")
  }
  if (length(metapop_from) != length(metapop_to)) {
    stop("to and from metapopulations are not same length.")
  }

  peterlist$travel <- rbind(
    peterlist$travel,
    tibble::tibble(
      frommetapopulation = metapop_from,
      tometapopulation = metapop_to,
      rate = rates,
      travelpops = travelstates
    )
  )
  return(peterlist)
}
