#' add interactions to model
#'
#' general form where users can specify states up
#' and down. Assumes per capita rates of interaction!
#'
#' @param peterlist list of instructions for piping |>
#' @param states_in character vector of base state names
#' that interact
#' @param states_out character vector of base
#' state names after interactions. Assume states_in
#' and states_out are ordered (e.g., S,I -> I,I)
#' @param rate character/numeric value specifying
#' PER CAPITA rate of interaction
#' @param normlogic logical specifying if rate
#' should be divided by total population in given
#' environment/metapopulation
#' @param environment character vector of environments
#' where this interaction occurs.
#' default is "" which specifies all environments
#' @param metapopulation character vector of
#' metapopulations where this interaction occurs
#' default is "" which specifies all metapopulations
#' @param name character value to allow scaling
#' by other functions
#' default is ""
#' @param groups specifies grouptypes the engage
#' in this interaction
#' default is "" which specifies all groups
#' @return updated instruction list
#' @export
add_interaction <- function(
    peterlist, states_in, states_out, rate,
    normlogic = TRUE, environment = "",
    metapopulation = "", name = "", groups = "") {
  if (("" %in% states_in) || ("" %in% states_in)) {
    stop("Empty state listed in interaction.
    Birth/death functionality currently not supported.
    Use add_import or add_export for importation/exportation.")
  }
  if (length(states_in) != length(states_out)) {
    stop("interaction can only lead to changes in state.
    states_in and states_out must be same length")
  }
  tibble2add <- namedlist2tibblerow(list(
    states_in = states_in,
    states_out = states_out,
    rate = rate,
    normlogic = normlogic,
    environment = environment,
    metapopulation = metapopulation,
    name = name,
    groups = groups
  ))
  peterlist$interactions <- rbind(peterlist$interactions, tibble2add)
  return(peterlist)
}
