#' Stratify population into groups
#'
#' @param peterlist list of instructions for piping |>
#' @param groupnames character vector of groups
#' in a given grouptype (e.g., c("Old","Young"))
#' @param grouptype character vector of grouptype (e.g., c("Age"))
#' default is NA which automatically generates a type during
#' model compilation
#' @param combinewithgrouptype character vector
#' specifying which group types to create combinations
#' of state names. default is "" which stratifies given
#' grouptype with all other grouptypes
#' @param basestates character vector of basestates
#' to stratify. default is "" which specifies all basestates
#' @param interactionscale character/numeric specifies
#' how to scale interaction rates for each group. default is 1
#' @param transitionscale character/numeric specifies
#' how to scale transition rates for each group. default is 1.
#' @return updated model instructions
#' @export
add_group <- function(
    peterlist, groupnames, grouptype = NA,
    combinewithgrouptype = "", basestates = "",
    interactionscale = 1, transitionscale = 1) {
  # maybe add metapopulation? split group in some
  # places and they never travel, dunno if needed
  # groupnames -- names of groups which
  # splut the basestate, must me a character vector
  # grouptype
  # should basestates be optional? when would
  # groups not have all states (if you don't vax certain groups?
  # that exhibit variation in type)

  # this can be simplified by pulling in column names
  if (!is.vector(groupnames)) {
    stop("Input groupnames must be a character or numeric vector")
  }
  if (!is.vector(combinewithgrouptype)) {
    stop("Input combinegrouptypes must be a character or numeric vector")
  }
  if (!is.vector(basestates)) {
    stop("Input basestates must be a character or numeric vector,
     applying to all groupnames for a  single grouptype")
  }
  if (length(grouptype) > 1) {
    stop("More than one grouptype specified, only one type (e.g., Age) may be
         specified even if inputting multiple groupnames (e.g.,c(Young,Old))")
  }

  peterlist$combinetypes[[length(peterlist$combinetypes) + 1]] <-
    c(grouptype, combinewithgrouptype)

  # coerce input to correct shape to apply
  numgroups <- length(groupnames)
  grouptype <- rep(grouptype, numgroups)
  interactionscale <- matchlength(interactionscale, numgroups)
  transitionscale <- matchlength(transitionscale, numgroups)

  bindthese <- list()
  for (seqidx in seq(numgroups)) {
    currnamedlist <- list(
      groupname = groupnames[[seqidx]],
      grouptype = grouptype[[seqidx]],
      basestates = basestates,
      interactionscale = interactionscale[[seqidx]],
      transitionscale = transitionscale[[seqidx]]
    )
    bindthese[[seqidx]] <- namedlist2tibblerow(currnamedlist)
  }
  tblout <- do.call("rbind", bindthese)
  peterlist$groups <- rbind(peterlist$groups, tblout)

  return(peterlist)
}
