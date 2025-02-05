#' add infection to model
#'
#' Helper function to add specific type of interaction
#' @param peterlist list of instructions for piping |>
#' @param infector base state that is infectious
#' @param infectee base state that is susceptible
#' @param infected base state that infectee
#' converts to following infectious
#' @param rate character/numeric value specifying
#' PER CAPITA rate of interaction
#' @param normlogic logical specifying if rate should be
#' divided by total population in given environment/metapopulation.
#' default is TRUE which scales rate by
#' total population (e.g., beta*S*I/N)
#' @param environment character vector of environments where
#' this interaction occurs.
#' default is "" which specifies all environments
#' @param metapopulation character vector of metapopulations
#' where this interaction occurs.
#' default is "" which specifies all metapopulations
#' @param name character value to allow scaling by other functions.
#' default is "infection"
#' @param groups specifies grouptypes that engage in this interaction.
#' default is "" which specifies all types.
#' @return updated instruction list
#' @export
add_infection <- function(
    peterlist, infector,
    infectee, infected, rate, normlogic = TRUE,
    environment = "", metapopulation = "",
    name = "infection", groups = "") {
  # infector -- state name of infectious class
  # infectee -- state name of susceptible class
  # infected -- state name for what susceptible becomes after infection
  # rate -- per capita rate e.g.,
  # normlogic -- normalize the rate by the total population
  peterlist <- add_interaction(peterlist,
    c(infector, infectee),
    c(infector, infected),
    rate = rate, normlogic = normlogic,
    environment = environment,
    metapopulation = metapopulation,
    name = name, groups = groups
  )
  return(peterlist)
}
