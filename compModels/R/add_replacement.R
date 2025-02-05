#' add replacement in model
#'
#' Helper fucntion that is same as transiton,
#' but allows simplifies specifying independent transitions
#' from multiple states to a single state.
#' @param peterlist list of instructions for piping |>
#' @param statesout character vector of state names
#' to transition from
#' @param statesin character vector of basestates to
#' transition to
#' @param meantime character/numeric value specifying
#' average time to go from fromstate to all tostates
#' @param name character value to allow scaling by
#' other functions
#' default is "replacement"
#' @param chainlength total steps from fromstate to
#' tostates (i.e., boxcars or length of pitchfork).
#' Alters distribution of transition times during
#' stochastic implementation
#' default is 1 for exponentially distributed wait times
#' @param chaintimescale character/numeric vector with
#' length chainlength scales rate transitioning between
#' steps along chains. Allows deviating transition time
#' distribution from a gamma distributions
#' default is NA which defaults to the same transition
#' times between all chained states
#' @param percapita specifies rates are given as per capita
#' value so the total change to population should be
#' scaled by the tostate
#' default to TRUE
#' @param metapopulation character vector of
#' metapopulation
#' names this transitions occurs in. default "" specifies
#' all metapopulations
#' default to "" which specifies all metapopulations
#' @param environment character vector of environment names
#' this transitions occurs in
#' default to "" which specifies all environments
#' @param groupname specifies which stratified groups transition
#' default to "" which specifies all groups.
#' @param grouptype specifies which stratified
#' group types transition
#' default to "" which specifies all grouptypes.
#' @return updated instruction list
#' @export
add_replacement <- function(
    peterlist,
    statesout,
    statesin,
    meantime,
    name = "replacement",
    chainlength = 1,
    chaintimescale = NA,
    percapita = TRUE,
    metapopulation = "",
    environment = "",
    groupname = "",
    grouptype = "") {
  if (length(statesin) != 1) {
    if (length(statesout) != length(statesin)) {
      stop("Error: vectors defining states
    replacing and replaced must be same length
    or the replacing vector must be length 1.")
    }
  } else {
    statesin < rep(statesin, length(statesout))
  }

  for (seqidx in seq_along(statesout)) {
    peterlist <- peterlist |>
      add_transition(statesout[[seqidx]],
        statesin[[seqidx]], meantime,
        name = name,
        chainlength = chainlength, chaintimescale = chaintimescale,
        percapita = percapita, metapopulation = metapopulation,
        environment = environment, groupname = groupname,
        grouptype = grouptype
      )
  }
}
