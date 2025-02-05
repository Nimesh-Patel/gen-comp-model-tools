1 #' Define metapopulation names
#'
#' Removes all prior add_metapopulation instructions and sets
#' them according to input
#'
#' @param peterlist list of instructions for piping |>
#' @param metapop_names character vector specifying
#' metapopulation names
#' default is ""
#' @param environment_names character vector specifying
#' environment names
#' default is ""
#' @param interactionscale numeric vector that scales
#' interactions in
#' each specified metapopulation
#' default is 1
#' @param transitionscale numeric vector that scales
#' transitions in
#' each specified metapopulation
#' default is 1
#' @param basestates character vector specifying
#' which basestates exist
#' in each metapopulations
#' default is "" which specifies all basestates
#' @return model instructions
#' @export
define_metapopulations <- function(peterlist,
                                   metapop_names = "",
                                   environment_names = "",
                                   interactionscale = 1,
                                   transitionscale = 1,
                                   basestates = "") {
  # Removes all prior metapopulations and sets
  # them according to input. There are multiple ways to define metapopulations
  # by defining "metapop_names" we give explit
  # names to the metapopulations e.g., adding 2 metpopulations c("US","UK")
  # environment_names sets environment names in metapopulations.
  # if its a character vector or list the same number
  # of elements as metapop_names then
  # there is 1 environment in each metapopulation.
  # If it is a list of lists, then the names in each
  # element list give multiple environments in the metapopulations
  # interactionscale -- scales interaction rates,
  # can be numeric or character (i.e., different infection rates
  # between locations),cmust be same size as metapop_names,
  # environment_names or a scalar
  # transitionscale -- scales transition rates, can be
  # numeric or character (e.g., parameter name). This allows,
  # different recovery rates between locations, consider the biological
  # reasonableness of using this. Must be same size as metapop_names,
  # environment_names or a scalar
  # basestates -- assigns which populations interact in given
  # environments or metapopulations

  coercedinputs <- make_metapopulations(
    metapop_names = metapop_names,
    environment_names = environment_names,
    interactionscale = interactionscale,
    transitionscale = transitionscale,
    basestates = basestates
  )
  peterlist$space <- tibble::as_tibble(coercedinputs)
  return(peterlist)
}
