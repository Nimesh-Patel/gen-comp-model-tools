#' Initializes model instructions with given base states
#'
#' @param base_state character vector of base states.
#' @param environment_name character vector environment names
#' default is ""
#' @param metapop_name character vector of metapopulation names
#' default is ""
#' @return list of model instructions
#' @export
define_states <- function(
    base_state,
    environment_name = "", metapop_name = "") {
  peterlist <- list()
  peterlist$states <- base_state

  peterlist$space <- tibble::tibble(
    metapop_names = metapop_name,
    environment_names = environment_name,
    interactionscale = 1,
    transitionscale = 1,
    basestates = list(base_state)
  )
  peterlist$travel <- tibble::tibble(
    frommetapopulation = character(),
    tometapopulation = character(),
    rate = character(),
    travelpops <- list()
  )

  peterlist$interactions <- tibble::tibble(
    states_in = list(),
    states_out = list(),
    ratefunction = character(),
    normlogic = logical(),
    scalefunction = character(),
    environment = character(),
    metapopulation = character(),
    name = character()
  )
  peterlist$transitions <- tibble::tibble(
    fromstate = character(),
    tobase = character(),
    tostate = character(),
    fromchain = numeric(),
    tochain = numeric(),
    rate = character(),
    percapitastate = character(),
    metapopulation = character(),
    environment = character(),
    name = character(),
    groupname = character(),
    grouptype = character()
  )
  peterlist$groups <- tibble::tibble(
    groupname = character(),
    grouptype = character(),
    basestates = list(),
    interactionscale = character(),
    transitionscale = character()
  )
  # named list to specify, '' for all interactions

  peterlist$combinetypes <- list()

  return(peterlist)
}
