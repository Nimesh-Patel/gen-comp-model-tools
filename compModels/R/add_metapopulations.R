#' Adds metapopulations to model instructions
#'
#' Allows nested input of metapopulations and
#' environments. Maybe clean this up.
#'
#' @param peterlist list of instructions for piping |>
#' @param metapop_names character vector of
#' metapopulation names to add.
#' default is "" where compilation automatically
#' generates a placeholder name for
#' each add_metapopulation use
#' @param environment_names character vector of
#' environment names to add
#' default is ""
#' @param interactionscale character/numeric vector
#' specifying multiplicative scaling of interaction
#' rates in each specifed metapopulation
#' default is 1
#' @param transitionscale character/numeric vector
#' specifying multiplicative scaling of transition
#' rates in each specifed metapopulation
#' default is 1
#' @param basestates specifies which basestates can
#' visit metapopulation/environment
#' default is "" which species all basestates
#' @return updated instruction list
#' @export
add_metapopulations <- function(
    peterlist, metapop_names = "",
    environment_names = "", interactionscale = 1,
    transitionscale = 1, basestates = "") {
  # keep old doc string for extra info...
  # Adds metapopulations to existing peterlist.
  # There are multiple ways to define metapopulations...
  # by defining "metapop_names" we give explit
  # names to the metapopulations e.g.,
  # adding 2 metpopulations c("US","UK")
  # environment_names sets environment names
  # in metapopulations. if its a character vector
  # or list the same number of elements as metapop_names then
  # there is 1 environment in each metapopulation.
  # If it is a list of lists, then the names in each
  # element list give multiple environments in the metapopulations
  # interactionscale -- scales interaction rates,
  # can be numeric or character (i.e., different
  # infection rates between locations), must be same
  # size as metapop_names, environment_names or a scalar
  # transitionscale -- scales transition rates, can
  # be numeric or character (e.g., parameter name).
  # This allows, different recovery rates between
  # locations, consider the biological reasonableness
  # of using this. Must be same size as metapop_names,
  # environment_names or a scalar
  # basestates -- assigns which populations
  # interact in given environments/metapopulations
  if ("" %in% peterlist$space$metapop_names) {
    warning("Prior metapopulation not named. Renaming to Metapopulation 1.")
    peterlist$space$metapop_names <- "Metapopulation 1"
  }
  newinputs <-
    tibble::as_tibble(make_metapopulations(
      metapop_names = metapop_names,
      environment_names =
        environment_names,
      interactionscale = interactionscale,
      transitionscale = transitionscale,
      basestates = basestates
    ))
  if ("" %in% newinputs) {
    filteredspace <- peterlist$space |>
      dplyr::filter(grepl("Metapopulation ", metapop_names))
    if (nrow(filteredspace) > 0) {
      # this will break if someone inputs "metapopulation_string",
      # where string is non-numeric. fix later
      maxidx <- max(as.numeric(sub(
        ".*Metapopulation ", "",
        filteredspace |>
          dplyr::pull(metapop_names)
      )))
      newidx <- maxidx + 1
      warning(paste0("Added metapopulation not named.
                     Renaming to metapopulation_", as.character(newidx), "."))
      newinputs$metapop_names <- paste0("Metapopulation ", as.character(newidx))
    } else {
      warning("Added metapopulation not named. Renaming to Metapopulation 1.")
      newinputs$metapop_names <- "Metapopulation 1"
    }
  }
  peterlist$space <- rbind(peterlist$space, newinputs)
  return(peterlist)
}
