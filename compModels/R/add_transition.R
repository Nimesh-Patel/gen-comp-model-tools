#' add transition in model
#'
#' These assume a "pitchfork" shape to transtions
#' E.g., recovery/vaccination/hospitalization
#' to 1 or multiple states
#' @param peterlist list of instructions for piping |>
#' @param fromstate basestate to transition from
#' @param tostates character vector of basestates to
#' transition to
#' @param meantime character/numeric value specifying
#' average time to go from fromstate to all tostates
#' @param name character value to allow scaling by
#' other functions
#' default to ""
#' @param chainlength total steps from fromstate to
#' tostates (i.e., boxcars or length of pitchfork)
#' Alters distribution of transition times during
#' stochastic implementation
#' default to 1 which is exponentially
#' distrubuted transition times
#' @param chaintimescale character/numeric vector
#' with length chainlength scales rate transitioning
#' between steps along chains. Allows deviating
#' transition time distribution from a gamma distributions
#' default to NA which is same transition times between
#' all chained states (e.g., gamma distributed wait times)
#' @param forkprobability character/numeric vector
#' specifying relative weight to transtion to
#' different tostates when multiple are given
#' default to NA which specifies equal probability
#' transitioning to each final state
#' @param percapita specifies rates are given as per
#' capita value so the total change to population
#' should be scaled by the tostate
#' default is TRUE
#' @param metapopulation character vector of metapopulation
#' names this transitions occurs in.
#' default "" is which specifies all metapopulations
#' @param environment character vector of environment
#' names this transitions occurs in.
#' default is "" which specifies all environments
#' @param groupname specifies which stratified groups
#' transition
#' default is "" which specifies all groups.
#' @param grouptype specifies which stratified group
#' types transition
#' default is "" which specifies all grouptypes.
#' @return updated instruction list
#' @export
add_transition <- function(
    peterlist, fromstate, tostates,
    meantime, name = "", chainlength = 1, chaintimescale = NA,
    forkprobability = NA, percapita = TRUE, metapopulation = "",
    environment = "", groupname = "", grouptype = "") {
  # Check consistency between and account for optional input
  if (length(fromstate) > 1) {
    stop("Transitions can only come from one state, multiple states are input.")
  }
  if (
    (!(NA %in% chaintimescale)) &&
      (chainlength > 1) &&
      (length(chaintimescale) != chainlength)
  ) {
    stop("Length of rate scaling between transitions
         don't match input number of transitions.")
  }
  if (!(NA %in% chaintimescale)) {
    pitchforklength <- length(chaintimescale)
    chaintimescale <- normalize2mean1(chaintimescale)
  } else {
    pitchforklength <- chainlength
  }
  if (!(NA %in% forkprobability)) {
    if (length(forkprobability) != length(tostates)) {
      stop("Length of states to transition into is
      different than the input forkprobability rates.")
    }
  } else {
    # assume equal probability of splitting between states
    if (length(tostates) > 1) {
      warning("Mulitple final states but relative rates (forkprobability)
      are not defined or contain NA. Assuming equal
      probability to transition between them")
    }
    forkprobability <- rep(1, length(tostates))
  }
  forkprobability <- normalize2probability(forkprobability)

  if (percapita) {
    percapitastate <- fromstate
  } else {
    percapitastate <- NA
  }
  baserate <- time2rate(meantime)
  # pitchfork shaft
  if (pitchforklength > 1) {
    baserate <- paste0(
      as.character(pitchforklength),
      "*",
      as.character(baserate)
    )
    if (!(NA %in% chaintimescale)) {
      # rates vary across pitchfork
      baserate <- paste0(
        as.character(chaintimescale),
        "*",
        as.character(baserate)
      )
    } else {
      baserate <- rep(as.character(baserate), pitchforklength)
    }
    toptibble <- tibble::tibble(
      fromstate = fromstate,
      tostate = fromstate,
      fromchain = seq(pitchforklength - 1),
      tochain = seq(pitchforklength - 1) + 1,
      rate = baserate[1:(pitchforklength - 1)],
      percapitastate = percapitastate,
      metapopulation = metapopulation,
      environment = environment,
      name = name,
      groupname = groupname,
      grouptype = grouptype
    )
  } else {
    toptibble <- tibble::tibble(
      fromstate = character(),
      tostate = character(),
      fromchain = numeric(),
      tochain = numeric(),
      rate = character(),
      percapitastate = character(),
      metapopulation = character(),
      environment = character(),
      name = name,
      groupname = groupname,
      grouptype = grouptype
    )
  }
  # pitchfork fork
  bottomtibble <- tibble::tibble(
    fromstate = fromstate,
    tostate = tostates,
    fromchain = pitchforklength,
    tochain = 1, # NA,
    rate = paste0(
      as.character(forkprobability),
      "*",
      baserate[length(baserate)]
    ),
    percapitastate = percapitastate,
    metapopulation = metapopulation,
    environment = environment,
    name = name,
    groupname = groupname,
    grouptype = grouptype
  )
  tibbleout <- rbind(toptibble, bottomtibble)
  peterlist$transitions <- rbind(peterlist$transitions, tibbleout)
  return(peterlist)
}
