#' Compiles model instructions to generate dynamical model
#' (updated state names, petersen matrix, and process rates)
#'
#'
#' @param inputpeter model instructions list.
#' @return list with updated instructions $modelinstructions to
#' aid user in verifying model and $modeloustrcutions list
#' with petersen matrix, process rates, and updated state names
#' @export
#' @importFrom rlang .data
compilemodel <- function(inputpeter) { # nolint: cyclocomp_linter.
  updatedstates <- inputpeter$states
  tblupdatestate <- tibble::tibble(
    basestates = character(),
    interactionscale = character(),
    transitionscale = character()
  )

  # initialize tables for reference
  # fill in default values for missing inputs
  allbasestates <- inputpeter$states

  tblgroups <- inputpeter$groups
  combinelist <- inputpeter$combinetypes
  dummylogic <- FALSE
  if (nrow(tblgroups) == 0) {
    dummylogic <- TRUE
    tblgroups <- tibble::tibble(
      groupname = "dummy",
      grouptype = "dummy",
      basestates = "",
      interactionscale = "",
      transitionscale = ""
    )
    combinelist <- list("dummy")
  }

  tblgroups <- cleanscale(tblgroups) |>
    dplyr::mutate(basestates = replace(
      .data$basestates,
      .data$basestates == "", list(allbasestates)
    ))
  if (is.logical(tblgroups$grouptype)) {
    tblgroups <- tblgroups |>
      dplyr::mutate(grouptype = as.character(.data$grouptype)) |>
      tidyr::replace_na(list(grouptype = "group"))
  }
  tbltravel <- inputpeter$travel
  tblspace <- inputpeter$space |>
    dplyr::mutate(basestates = replace(
      .data$basestates,
      .data$basestates == "", list(allbasestates)
    ))
  allmetapops <- tblspace |>
    dplyr::distinct(.data$metapop_names) |>
    dplyr::pull(.data$metapop_names)
  tblspace <- cleanscale(tblspace) |>
    dplyr::rename(
      spacetransitionscale = .data$transitionscale,
      spaceinteractionscale = .data$interactionscale
    )
  tblspace_unnest <- tblspace |> tidyr::unnest(.data$basestates)
  tbltransition <- inputpeter$transitions
  tbltransition_expand <- tbltransition
  tblinteractions <- inputpeter$interactions |>
    dplyr::rename(grouptype = .data$groups)

  distinctgrouptypes <- tblgroups |>
    dplyr::distinct(.data$grouptype) |>
    dplyr::pull()
  numgrouptype <- length(distinctgrouptypes)
  if (numgrouptype != 0) {
    if (!dummylogic) {
      print("updating names -- groups")
    }
    if (numgrouptype == 1) {
      combinelist <- list(distinctgrouptypes)
      listgroupcombinations <- lapply(
        combinelist,
        groupcombinationsacrosstypes,
        tblgroups
      )
    } else {
      # WHAT IF: no types?
      # replace empty string with all types
      combinelist <- lapply(combinelist, empty2states, distinctgrouptypes)
      # remove duplicate types in types, sort, and keep unique sets of types
      # combinelist is a list of all combination of groups and
      # is used throughout the script to determine consistency of input
      combinelist <- unique(lapply(combinelist, function(x) {
        sort(unique(x))
      }))

      # check prescribed basestates are consistent within combined types
      for (currtypes in combinelist) {
        # replace "" and sort and get distinct elements
        basestates2expandtbl <- tblgroups |>
          dplyr::filter(.data$grouptype %in% currtypes) |>
          dplyr::mutate(basestates = replace(
            .data$basestates,
            .data$basestates == "",
            list(allbasestates)
          )) |>
          dplyr::mutate(basestates = lapply(.data$basestates, sort)) |>
          dplyr::distinct(.data$basestates)

        if (nrow(basestates2expandtbl) > 1) {
          stop(paste0(
            "The {",
            paste(currtypes, collapse = ","),
            "} types do not share the same basestates.
                        Check when adding groups."
          ))
        }
      }
      # make combinations of groups for each grouptype combinations
      listgroupcombinations <- lapply(
        combinelist,
        groupcombinationsacrosstypes,
        tblgroups
      )
    }
    tblupdatestate <-
      dplyr::bind_rows(lapply(
        listgroupcombinations,
        function(x) {
          x |> tidyr::unnest(.data$basestates)
        }
      ))

    currtransitionexpandlist <- list()
    for (listidx in seq_along(combinelist)) {
      currtransitionexpandlist[[listidx]] <- tbltransition_expand |>
        replaceglobal("grouptype", combinelist[[listidx]])
    }
    tbltransition_expand <- dplyr::bind_rows(currtransitionexpandlist)
    tbltransition_expand_split <- tbltransition_expand |>
      dplyr::group_split(.data$grouptype)

    currsplit2typecols <- function(currtbl) {
      tblout <- split2typecols(currtbl, combinelist, listgroupcombinations)
      return(tblout)
    }

    tbltransition_expand <- dplyr::bind_rows(lapply(
      tbltransition_expand_split,
      currsplit2typecols
    )) |>
      dplyr::select(-tidyselect::one_of(
        "groupname",
        "transitionscale",
        "interactionscale",
        "basestates"
      )) |>
      suppressWarnings()
  }
  # update with basestates without groups
  statesintblupdatestate <- tblupdatestate |>
    dplyr::distinct(.data$basestates) |>
    dplyr::pull()
  unaccountedstates <- setdiff(inputpeter$states, statesintblupdatestate)
  # add replace NA term to deal without groups
  tblupdatestate <- dplyr::bind_rows(
    tblupdatestate,
    tibble::tibble(basestates = unaccountedstates)
  ) |>
    tidyr::replace_na(list(interactionscale = "", transitionscale = ""))
  # current columns:
  # basestates,interactionscale,transitionscale and optional grouptypes


  tblspacelogic <- nrow(tblspace) > 1
  if (tblspacelogic) {
    print("names -- space")
  }

  tbltransition_expand <- tbltransition_expand |>
    replaceglobal("metapopulation", allmetapops) |>
    tidyr::unnest(.data$metapopulation)

  tblupdatestate <- tblupdatestate |>
    dplyr::full_join(tblspace_unnest,
      by = "basestates",
      relationship = "many-to-many"
    ) |>
    dplyr::mutate(
      interactionscale = paste0(
        .data$spaceinteractionscale,
        .data$interactionscale
      ),
      transitionscale = paste0(
        .data$spacetransitionscale,
        .data$transitionscale
      )
    ) |>
    dplyr::select(-.data$spacetransitionscale, -.data$spaceinteractionscale)

  # Update names based on chaining/boxcars
  # SI model doesn't have transitions
  if (nrow(tbltransition_expand) > 0) {
    # check whether to update names due to chaining
    # switched tbltransition -> tbltransition_expand
    haschains <- tbltransition_expand |> dplyr::filter(.data$fromchain > 1)
    if (nrow(haschains) > 0) {
      print("names -- chained transition")
      lastchains <- haschains |> dplyr::filter(.data$fromstate != .data$tostate)
      # add functionality for multiple environments below

      countlastchains <- lastchains |>
        dplyr::count(
          .data$fromstate,
          .data$tostate,
          .data$metapopulation,
          dplyr::across(tidyselect::all_of(distinctgrouptypes))
        )

      # check if multiple long transitions go to and from same pair
      # of states in a given metapopulation
      # if so, throw error.
      # later allow combinations and demand processes are
      # given different names to update states
      # e.g.,"name1chain1_name2chain1"
      # "name1chain2_name2chain1"
      # "name1chain1_name2chain2"
      # "name1chain2_name2chain2"
      if (nrow(countlastchains |> dplyr::filter(.data$n > 1)) > 0) {
        stop("Two or more different chainlengths for transitions from
              the same fromstate among the same grouptype (e.g., two
              different gamma distributions are defined for transitioning
              from the same state). This is illogical, review transitions
              and consider changing state names")
      }
      # lastchain has relevant final fromchain
      # fill from chain
      lastchains2expand <- lastchains |>
        dplyr::select(
          -.data$tostate,
          -.data$tochain,
          -.data$rate,
          -.data$percapitastate,
          -.data$name,
          -.data$environment
        ) |>
        dplyr::distinct()
      # clean by making consistent names
      lastchainsexpanded <- lastchains2expand |>
        dplyr::mutate(dplyr::across(
          .data$fromchain,
          ~ map(.x, ~ seq_len(.x))
        )) |>
        tidyr::unnest(.data$fromchain) |>
        dplyr::rename(
          basestates = .data$fromstate,
          metapop_names = .data$metapopulation
        )

      # maybe remove below
      tbltransition_expand_distinct <- tbltransition_expand |>
        dplyr::select(
          -.data$tostate,
          -.data$tochain,
          -.data$rate,
          -.data$percapitastate,
          -.data$name,
          -.data$environment
        ) |>
        dplyr::distinct() |>
        dplyr::rename(
          basestates = .data$fromstate,
          metapop_names = .data$metapopulation
        )
      # maybe remove above
      tblupdatestate <- dplyr::left_join(
        tblupdatestate,
        lastchainsexpanded
      ) |>
        dplyr::rename(chain = .data$fromchain)
      # replace_na is used to expand in next step
    } else {
      tblupdatestate <- tblupdatestate |>
        dplyr::mutate(chain = NA) # switched from 1
    }
  }

  tblupdatestate <- tblupdatestate |>
    dplyr::rename(metapopulation = .data$metapop_names)
  # sort column names for consistency and build new statenames
  # breaks if no groups -- keep old pipes for debugging
  tbl2bind <- tblupdatestate |>
    dplyr::select(
      -.data$interactionscale, -.data$transitionscale,
      -.data$environment_names
    )
  # remove columns with only 1 distinct value -- but keep basestates
  # this automatically removes any dummy variable
  currleft <- tbl2bind |> dplyr::select(.data$basestates)
  currright <- tbl2bind |> dplyr::select(-.data$basestates)
  keeplogic <- currright |> dplyr::summarise_all(dplyr::n_distinct) != 1
  currright <- currright |> dplyr::select_if(function(x) length(unique(x)) > 1)
  tblupdatestate2name <- cbind(currleft, currright)

  newvariablenames <- c()
  if (ncol(tblupdatestate2name) > 1) {
    for (rowidx in seq_len(nrow(tblupdatestate2name))) {
      currrow <- tblupdatestate2name[rowidx, ]
      newvariablenames <- c(
        newvariablenames,
        makenewvariablenamesfromrow(currrow)
      )
    }
  } else {
    newvariablenames <- tblupdatestate$basestates
  }
  tblupdatestate$updatedstate <- newvariablenames


  # UPDATE PETERSEN & PROCESSES -- travel, transitions, interactions
  print("start")
  updatedstatenames <- tblupdatestate |>
    dplyr::distinct(.data$updatedstate) |>
    dplyr::pull(.data$updatedstate)
  currpetercolfun <- function(x, y) {
    states2petercolumn(x, y, updatedstatenames)
  }
  allpetermat <- list()
  allprocessvec <- list()
  # TRAVEL -- SPACE

  if (nrow(tbltravel) > 0) {
    print("Adding travel to dynamics")
    tbltravel <- tbltravel |>
      dplyr::mutate(frommetapopulation = replace(
        .data$frommetapopulation,
        .data$frommetapopulation == "",
        list(allmetapops)
      )) |>
      tidyr::unnest(.data$frommetapopulation) |>
      dplyr::mutate(tometapopulation = replace(
        .data$tometapopulation,
        .data$tometapopulation == "",
        list(allmetapops)
      )) |>
      tidyr::unnest(.data$tometapopulation) |>
      dplyr::filter(.data$frommetapopulation != .data$tometapopulation) |>
      dplyr::mutate(travelpops = replace(
        .data$travelpops,
        .data$travelpops == "",
        list(allbasestates)
      )) |>
      tidyr::unnest(.data$travelpops)

    # environment names
    tbl2jointotbltravel <- tblupdatestate |>
      dplyr::distinct_at(dplyr::vars(
        -.data$environment_names,
        -.data$transitionscale,
        -.data$interactionscale
      )) |>
      dplyr::full_join(tblupdatestate,
        by = colnames(tblupdatestate |>
          dplyr::select(
            -.data$metapopulation,
            -.data$updatedstate,
            -.data$transitionscale,
            -.data$interactionscale,
            -.data$environment_names
          ))
      ) |>
      dplyr::filter(.data$metapopulation.y != .data$metapopulation.x) |>
      dplyr::rename(
        frommetapopulation = .data$metapopulation.x,
        tometapopulation = .data$metapopulation.y,
        fromstate = .data$updatedstate.x,
        tostate = .data$updatedstate.y
      )

    tbltravel2 <- tbltravel |>
      dplyr::left_join(
        tbl2jointotbltravel |>
          dplyr::select(
            .data$frommetapopulation,
            .data$tometapopulation,
            .data$fromstate,
            .data$tostate,
            .data$basestates
          ),
        by = c("frommetapopulation",
          "tometapopulation",
          "travelpops" = "basestates"
        )
      )

    petersenmat_travel <- do.call(
      cbind,
      mapply(currpetercolfun, tbltravel2 |> dplyr::pull(.data$fromstate),
        tbltravel2 |> dplyr::pull(.data$tostate),
        USE.NAMES = FALSE
      )
    )
    # updaterates? TODO
    tbltravel_processstr <- paste0(tbltravel2$fromstate, "*", tbltravel2$rate)

    allpetermat[["travel"]] <- petersenmat_travel
    allprocessvec[["travel"]] <- tbltravel_processstr
  }

  # transition
  if (nrow(tbltransition) > 0) {
    print("Adding transitions to dynamics")
    # check scaling
    # split into transitions from chain 1 <-
    ttrans_11 <- tbltransition_expand |>
      dplyr::filter(.data$fromchain == 1, .data$tochain == 1)

    ttrans_11 <- ttrans_11 |>
      dplyr::left_join(tblupdatestate |> tidyr::replace_na(list(chain = 1)),
        by = c(
          "fromstate" = "basestates", "environment" = "environment_names",
          "metapopulation", distinctgrouptypes
        )
      ) |>
      tidyr::replace_na(list(transitionscale = "")) |>
      dplyr::mutate(rate = paste0(.data$transitionscale, .data$rate)) |>
      dplyr::select(
        -.data$transitionscale, -.data$interactionscale,
        -.data$fromchain
      ) |>
      dplyr::rename(
        fromupdatedstate = .data$updatedstate,
        fromchain = .data$chain
      ) |>
      dplyr::left_join(tblupdatestate |> tidyr::replace_na(list(chain = 1)),
        by = c(
          "tostate" = "basestates",
          "tochain" = "chain",
          "environment" = "environment_names",
          "metapopulation", distinctgrouptypes
        )
      ) |>
      dplyr::select(-.data$transitionscale, -.data$interactionscale) |>
      dplyr::rename(toupdatedstate = .data$updatedstate)


    ttrans_not11 <- tbltransition_expand |>
      dplyr::filter(.data$fromchain != 1 | .data$tochain != 1)

    ttrans_not11 <- ttrans_not11 |>
      dplyr::left_join(tblupdatestate |> tidyr::replace_na(list(chain = 1)),
        by = c(
          "fromstate" = "basestates",
          "fromchain" = "chain", "environment" = "environment_names",
          "metapopulation", distinctgrouptypes
        )
      ) |>
      tidyr::replace_na(list(transitionscale = "")) |>
      dplyr::mutate(rate = paste0(.data$transitionscale, .data$rate)) |>
      dplyr::select(-.data$transitionscale, -.data$interactionscale) |>
      dplyr::rename(fromupdatedstate = .data$updatedstate) |>
      dplyr::left_join(tblupdatestate |> tidyr::replace_na(list(chain = 1)),
        by = c(
          "tostate" = "basestates",
          "tochain" = "chain",
          "environment" = "environment_names",
          "metapopulation", distinctgrouptypes
        )
      ) |>
      dplyr::select(-.data$transitionscale, -.data$interactionscale) |>
      dplyr::rename(toupdatedstate = .data$updatedstate)

    tbltransition_expand <- rbind(ttrans_11, ttrans_not11)

    petersenmat_transition <- do.call(
      cbind,
      mapply(currpetercolfun,
        tbltransition_expand |> dplyr::pull(.data$fromupdatedstate),
        tbltransition_expand |> dplyr::pull(.data$toupdatedstate),
        USE.NAMES = FALSE
      )
    )

    tbltransition_processstr <- tbltransition_expand$rate
    # change for just a logical
    notnalogic <- !is.na(tbltransition_expand$percapitastate)
    tbltransition_processstr[notnalogic] <-
      paste0(
        tbltransition_expand$fromupdatedstate[notnalogic],
        "*",
        tbltransition_processstr[notnalogic]
      )

    # update for transitions
    allpetermat[["transition"]] <- petersenmat_transition
    allprocessvec[["transition"]] <- tbltransition_processstr
  }

  # expand by groups if necessary
  if (nrow(tblinteractions) > 0) {
    print("Adding interactions to dynamics")
    # check that interactions are interpretable
    # same number of states_in and states_out
    # for now remove interactions with different length instead
    # of padding stop(Don't have functionality
    # for birth death / immigration)
    tblinteractions_expand <- tblinteractions
    if (!identical(
      sapply(tblinteractions$states_in, length),
      sapply(tblinteractions$states_in, length)
    )) {
      stop("Interaction states_in and states_out have different
      lengths. Only change of states are allowed, not birth/death")
    }
    emptyinlogic <- "" %in% unlist(tblinteractions_expand$states_in)
    emptyoutlogic <- "" %in% unlist(tblinteractions_expand$states_out)
    if (emptyinlogic || emptyoutlogic) {
      stop("Specified states_in or states_out can not be empty string")
    }
    tblinteractions_expand$changestateidx <- mapply(function(x, y) {
      which(x != y)
    }, tblinteractions_expand$states_in, tblinteractions_expand$states_out)
    # split by metapopulation
    tblinteractions_expand <- tblinteractions_expand |>
      replaceglobal("metapopulation", allmetapops) |>
      tidyr::unnest(.data$metapopulation)
    # then maybe environments

    # scale rates
    tblinteractions_expand <- tblinteractions_expand |>
      replaceglobal("grouptype", distinctgrouptypes)
    rowpetermat <- list()
    rowprocessvec <- list()
    for (rowidx in seq_len(nrow(tblinteractions_expand))) {
      currrow <- tblinteractions_expand[rowidx, ]
      currmetapop <- currrow$metapopulation[[1]]
      currenv <- currrow$environment[[1]]
      currgrouptypes <- currrow$grouptype[[1]]
      currrowrate <- currrow$rate[[1]]
      # is it necessary to check NA?
      select4logic <- tblupdatestate |>
        dplyr::select(tidyselect::any_of(currgrouptypes))
      templogic <- !is.na(select4logic)
      currtblupdatestate <- tblupdatestate |>
        dplyr::filter(
          Matrix::rowSums(templogic) > 0,
          .data$metapopulation == currmetapop,
          .data$environment_names == currenv
        )
      # get number of rows for each
      currstatesin <- currrow$states_in[[1]]
      currstatesintbls <- lapply(currstatesin, function(x) {
        currtblupdatestate |> dplyr::filter(.data$basestates == x)
      })
      currstatesout <- currrow$states_out[[1]]
      currstatesouttbls <- lapply(currstatesout, function(x) {
        currtblupdatestate |> dplyr::filter(.data$basestates == x)
      })
      numstatesin <- sapply(currstatesintbls, nrow)
      currgeoscale <- length(numstatesin)
      combo_idxs <- expand.grid(lapply(numstatesin, function(x) {
        seq(x)
      }))

      # process rate
      currprocessscale <- rep(
        paste0("1)^(1/", as.character(currgeoscale), ")"), nrow(combo_idxs)
      )

      for (stateinidx in seq_along(currstatesintbls)) {
        if (stateinidx == 1) {
          currcomboidx <- combo_idxs[, stateinidx]
          currpercapscale <-
            currstatesintbls[[stateinidx]]$updatedstate[currcomboidx]
          currprocessscale <- paste0(
            currstatesintbls[[stateinidx]]$interactionscale[currcomboidx],
            currprocessscale
          )
        } else {
          currcomboidx <- combo_idxs[, stateinidx]
          currpercapscale <- paste0(
            currstatesintbls[[stateinidx]]$updatedstate[currcomboidx],
            "*", currpercapscale
          )
          currprocessscale <- paste0(
            currstatesintbls[[stateinidx]]$interactionscale[currcomboidx],
            currprocessscale
          )
        }
      }
      currprocessscale <- paste0("(", currprocessscale)
      currrowprocessvec <- paste0(
        paste0(currprocessscale, "*", currpercapscale),
        "*", currrowrate
      )

      if (currrow$normlogic[[1]]) {
        normstr <- paste(currtblupdatestate$updatedstate, collapse = "+")
        currrowprocessvec <- paste0(currrowprocessvec, "/(", normstr, ")")
      }
      rowprocessvec[[rowidx]] <- currrowprocessvec
      # petermat
      # initialize
      currrowpetermat <- Matrix::Matrix(0,
        nrow = length(updatedstatenames),
        ncol = length(currrowprocessvec),
        sparse = TRUE
      )
      allchangeidx <- currrow$changestateidx[[1]]
      for (currchangeidx in allchangeidx) {
        nowcurrstatesintbl <- currstatesintbls[[currchangeidx]] |>
          dplyr::select(
            -.data$transitionscale, -.data$interactionscale,
            -.data$basestates
          )
        nowcurrstatesouttbl <- currstatesouttbls[[currchangeidx]] |>
          dplyr::filter((.data$chain == 1) | is.na(.data$chain)) |>
          dplyr::select(
            -.data$transitionscale,
            -.data$interactionscale, -.data$chain, -.data$basestates
          )
        col2join <-
          colnames(nowcurrstatesouttbl |> dplyr::select(-.data$updatedstate))
        tbl2join <-
          dplyr::left_join(nowcurrstatesintbl,
            nowcurrstatesouttbl,
            by = col2join
          ) |>
          dplyr::rename(
            statedown = .data$updatedstate.x,
            stateup = .data$updatedstate.y
          ) |>
          dplyr::select(.data$statedown, .data$stateup)
        tmpcmbidx <- combo_idxs[, currchangeidx]
        currchangeprocesstbl <-
          tibble::tibble(
            processnumber = seq(1, length(currrowprocessvec)),
            statedown = nowcurrstatesintbl$updatedstate[tmpcmbidx]
          ) |>
          dplyr::left_join(tbl2join) |>
          dplyr::mutate(rowdown = sapply(.data$statedown, function(x) {
            which(updatedstatenames == x)
          }), rowup = sapply(.data$stateup, function(x) {
            which(updatedstatenames == x)
          }))
        for (j in seq_len(nrow(currchangeprocesstbl))) {
          currrowdown <- currchangeprocesstbl$rowdown[[j]]
          currrowup <- currchangeprocesstbl$rowup[[j]]
          currprocessnum <- currchangeprocesstbl$processnumber[[j]]
          currrowpetermat[currrowdown, currprocessnum] <-
            currrowpetermat[currrowdown, currprocessnum] - 1
          currrowpetermat[currrowup, currprocessnum] <-
            currrowpetermat[currrowup, currprocessnum] + 1
        }
      }
      rowpetermat[[rowidx]] <- currrowpetermat
    }
    allpetermat[["interaction"]] <- do.call(cbind, rowpetermat)
    allprocessvec[["interaction"]] <- do.call(c, rowprocessvec)
  }

  # output these
  # updatedstate -- cheatsheet for all states
  # updatedprocess -- cheatsheet for process, include basestate
  # updatedpetermatrix -- matrix ordered by updatedstates
  if (dummylogic) {
    tblupdatestate <- tblupdatestate |> dplyr::select(-.data$dummy)
    tblinteractions_expand <- tblinteractions_expand |>
      dplyr::select(-.data$grouptype)
    tbltransition_expand <- tbltransition_expand |> dplyr::select(-.data$dummy)
  }
  peterlist <- list()
  peterlist$updatedstates <- updatedstatenames
  peterlist$petermatrix <- do.call(cbind, allpetermat)
  peterlist$processrates <- do.call(c, allprocessvec)
  modelinstructions <- list()
  modelinstructions$tblupdatedstates <- tblupdatestate
  modelinstructions$tbltransitions <- tbltransition_expand
  modelinstructions$tblinteractions <- tblinteractions_expand
  modelinstructions$tblspace <- tblspace
  outputlist <- list(
    modelinstructions = modelinstructions,
    modeloutstructions = peterlist
  )
  return(outputlist)
}
