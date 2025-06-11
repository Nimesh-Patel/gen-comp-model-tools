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
    scaleinteractions = character(),
    scaletransitions = character()
  )

  # PREP TABLES FOR COMPILATION
  # e.g., fill in default values for missing inputs, fill in global values, etc.
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
      scaleinteractions = "",
      scaletransitions = "",
      scalemigrations = "",
      scaleprocessbyname = list(list()),
      scaleprocessbygroup = list(list())
    )
    combinelist <- list("dummy")
  }

  # get list of groups by type to replace ""
  temptbl <- tblgroups |>
    dplyr::select("groupname", "grouptype") |>
    dplyr::group_by(.data$grouptype) |>
    dplyr::summarize(tempcol = list(.data$groupname)) |>
    dplyr::ungroup()
  allgroupslist <- temptbl[["tempcol"]]
  names(allgroupslist) <- temptbl$grouptype

  tblgroups <- cleanscale(tblgroups) |>
    dplyr::mutate(basestates = replace(
      .data$basestates,
      .data$basestates == "",
      list(allbasestates)
    ))

  tbltravel <- inputpeter$travel
  tblprocess_travel <- tbltravel

  tblspace <- inputpeter$space |>
    dplyr::mutate(
      basestates = replace(
        .data$basestates,
        .data$basestates == "",
        list(allbasestates)
      ),
      groups = replace(
        .data$groups,
        .data$groups == "",
        list(allgroupslist)
      )
    )

  allmetapops <- tblspace |>
    dplyr::distinct(.data$metapopulation) |>
    dplyr::pull("metapopulation")

  tblspace <- cleanscale(tblspace) |>
    dplyr::rename(
      spacescaletransitions = "scaletransitions",
      spacescaleinteractions = "scaleinteractions"
    )
  tblspace_unnest <- tblspace |> tidyr::unnest("basestates")
  tbltransition <- inputpeter$transitions
  tbltransition_expand <- tbltransition |>
    dplyr::mutate(
      metapopulation = replace(
        .data$metapopulation,
        .data$metapopulation == "",
        list(allmetapops)
      ),
      fromgroups = replace(
        .data$fromgroups,
        .data$fromgroups == "",
        list(allgroupslist)
      ),
      togroups = replace(
        .data$togroups,
        .data$togroups == "",
        list(allgroupslist)
      )
    ) |>
    tidyr::unnest("metapopulation")

  tblinteractions <- inputpeter$interactions

  distinctgrouptypes <- tblgroups |>
    dplyr::distinct(.data$grouptype) |>
    dplyr::pull()
  numgrouptype <- length(distinctgrouptypes)


  # COMPILE STATES
  print("Compiling names")

  # UPDATE NAMES -- GROUPS
  if (!dummylogic) {
    print("updating names -- groups")
  }
  if (numgrouptype == 1) {
    combinelist <- list(distinctgrouptypes)
    listgroupcombinations <- lapply(
      combinelist,
      groupcombinations,
      tblgroups
    )
  } else {
    # replace empty string with all types
    combinelist <- lapply(combinelist, empty2states, distinctgrouptypes)
    # remove duplicate types in types, sort, and keep unique sets of types
    # combinelist is a list of all combination of groups and
    # is used throughout the script to determine consistency of input
    combinelist <- unique(lapply(
      combinelist,
      function(x) {
        sort(unique(x))
      }
    ))
    allcombtypes <- unlist(combinelist)
    missingcombtypes <- setdiff(distinctgrouptypes, allcombtypes)
    combinelist <- c(combinelist, lapply(missingcombtypes, c))
    # check prescribed basestates are consistent within combined types
    for (currtypes in combinelist) {
      # replace "" and sort and get distinct elements
      basestates2expandtbl <- tblgroups |>
        dplyr::filter(.data$grouptype %in% currtypes) |>
        dplyr::mutate(basestates = lapply(.data$basestates, sort)) |>
        dplyr::distinct(.data$basestates)

      if (nrow(basestates2expandtbl) > 1) {
        stop(paste0(
          "The {",
          paste(currtypes, collapse = ","),
          "} types do not share the same basestates. Check when
                    adding groups."
        ))
      }
    }
    # make combinations of groups for each grouptype combinations
    listgroupcombinations <- lapply(combinelist, groupcombinations, tblgroups)
  }
  tblupdatestate <- dplyr::bind_rows(lapply(
    listgroupcombinations,
    function(x) {
      x |>
        tidyr::unnest("basestates")
    }
  ))

  # update with basestates without groups
  statesintblupdatestate <- tblupdatestate |>
    dplyr::distinct(.data$basestates) |>
    dplyr::pull()
  unaccountedstates <- setdiff(inputpeter$states, statesintblupdatestate)
  tblupdatestate <- dplyr::bind_rows(
    tblupdatestate,
    tibble::tibble(
      basestates =
        unaccountedstates
    )
  )


  # UPDATE NAMES -- METAPOPULATIONS
  tblspacelogic <- nrow(tblspace) > 1
  if (tblspacelogic) {
    print("updating names -- metapopulations")
  }
  tblspace2expand <-
    tblspace |> dplyr::select("metapopulation", "basestates", "groups")
  tblspace2expand_base <- tblspace2expand |>
    tidyr::unnest("basestates") |>
    dplyr::select("metapopulation", "basestates")
  tblspace2expand_group <- lapply(
    tblspace2expand$groups,
    grabfromtbl,
    tblupdatestate
  )
  for (metapopidx in seq_len(nrow(tblspace2expand))) {
    tblspace2expand_group[[metapopidx]] <-
      tblspace2expand_group[[metapopidx]] |>
      dplyr::mutate(
        metapopulation = tblspace2expand$metapopulation[[metapopidx]]
      )
  }
  tblspace2expand_group <- dplyr::bind_rows(tblspace2expand_group)
  tblupdatestate <-
    dplyr::full_join(tblspace2expand_group, tblspace2expand_base)

  # UPDATE NAMES -- CHAINED TRANSITIONS
  tbltransition_last <- tbltransition_expand |>
    dplyr::filter(.data$fromchain > 1, .data$tochain == 1) |>
    dplyr::distinct(.data$fromstate,
      .data$compileid,
      .data$metapopulation,
      .keep_all = TRUE
    )
  # update if there are chained transitions
  # this should be same for changing groups
  if (nrow(tbltransition_last) > 0) {
    print("updating names -- chained transitions")
    # get compileid step comparable to type group
    tblchains <- tbltransition_last |>
      dplyr::distinct(
        .data$fromchain,
        .data$compileid,
        .data$metapopulation
      ) |>
      dplyr::mutate(
        step = purrr::map(.data$fromchain, ~ seq_len(.x)),
        compileid = paste0("chainid", as.character(.data$compileid))
      ) |>
      tidyr::unnest("step") |>
      dplyr::mutate(step = as.character(.data$step)) |>
      dplyr::select(-"fromchain")

    tbltransition_grab <- lapply(
      tbltransition_last$fromgroups,
      grabfromtbl, tblupdatestate
    )

    # loop because dplyr seems to not work as expected
    for (loopidx in seq_along(tbltransition_grab)) {
      tbltransition_grab[[loopidx]] <- tbltransition_grab[[loopidx]] |>
        dplyr::filter(
          .data$basestates == tbltransition_last$fromstate[[loopidx]],
          .data$metapopulation == tbltransition_last$metapopulation[[loopidx]]
        ) |>
        dplyr::mutate(compileid = tbltransition_last$compileid[[loopidx]])
    }
    tbltransition_chaingroup <- dplyr::bind_rows(tbltransition_grab)

    bynames <- setdiff(
      names(tbltransition_chaingroup),
      c("compileid", "dummykillname")
    )

    # remove by left_joining in case multiple transitions from same state
    tbltransition_chaingroup <-
      dplyr::bind_rows(tbltransition_chaingroup, tblupdatestate) |>
      dplyr::left_join(
        tbltransition_chaingroup |>
          dplyr::select(-"compileid") |>
          dplyr::mutate(dummykillname = TRUE),
        by = bynames
      ) |>
      dplyr::filter(!((is.na(.data$compileid)) & (!is.na(.data$dummykillname)))) |> # nolint: line_length_linter.
      dplyr::select(-"dummykillname")

    tblchains_na <- tbltransition_chaingroup |>
      dplyr::filter(is.na(.data$compileid)) |>
      dplyr::select(-"compileid")
    tblchains_notna <- tbltransition_chaingroup |>
      dplyr::filter(!is.na(.data$compileid))


    tblchains_notnagroup <- tblchains_notna |>
      dplyr::mutate(compileid = paste0(
        "chainid",
        as.character(.data$compileid)
      )) |>
      dplyr::group_by(dplyr::across(-"compileid")) |>
      dplyr::summarise(compileid = list(.data$compileid)) |>
      dplyr::ungroup()

    # get chain combinations, in a loop for ease
    listtosquash <- list()
    for (rowidx in seq_len(nrow(tblchains_notnagroup))) {
      currrow <- tblchains_notnagroup[rowidx, ]
      currtblchains <- tblchains |>
        dplyr::filter(.data$metapopulation == currrow$metapopulation[[1]])
      add2row <- grouptypecombinations(
        currrow$compileid[[1]],
        currtblchains,
        "compileid",
        "step"
      )
      currrow_rep <- dplyr::bind_rows(replicate(nrow(add2row),
        currrow |>
          dplyr::select(-"compileid"),
        simplify = FALSE
      ))
      listtosquash[[rowidx]] <- dplyr::bind_cols(currrow_rep, add2row)
    }
    tblchains_notstack <- dplyr::bind_rows(listtosquash)

    tblupdatestate <- dplyr::bind_rows(tblchains_notstack, tblchains_na)
  }


  # remove columns with only 1 distinct value -- but keep basestates
  # this automatically removes any dummy variable
  tbl2bind <- tblupdatestate
  currleft <- tbl2bind |> dplyr::select("basestates")
  currright <- tbl2bind |> dplyr::select(-"basestates")
  currright <- currright |> dplyr::select_if(function(x) length(unique(x)) > 1)
  tblupdatestate2name <- dplyr::bind_cols(currleft, currright)
  tblstates <- tblupdatestate2name
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
  tblstates$updatedstate <- newvariablenames

  # combine by metapopulation to get total
  tblntotal <- tblupdatestate |>
    dplyr::group_by(.data$metapopulation) |>
    dplyr::summarize(totalpop = paste0(.data$updatedstate, collapse = "+")) |>
    dplyr::mutate(
      totalpop = paste0("(", .data$totalpop, ")"),
      varname = paste0("Ntotal_", make.names(.data$metapopulation))
    )
  if (nrow(tblntotal) == 1) {
    tblntotal$varname <- "Ntotal"
  }

  # SCALE PROCESSES
  tblprocesslist <- list(
    transition = tbltransition,
    migration = tbltravel,
    interaction = tblinteractions
  )
  tblprocessnames <-
    dplyr::bind_rows(lapply(
      names(tblprocesslist),
      function(x) {
        tblprocesslist[[x]] |>
          dplyr::distinct(
            .data$processname,
            .data$processgroup
          ) |>
          dplyr::mutate(processlabel = x)
      }
    ))

  # remove non-existent processes
  tblprocessnames <- tblprocessnames |>
    dplyr::filter(
      !((is.na(tblprocessnames$processname)) &
        (is.na(tblprocessnames$processgroup)) # nolint: indentation_linter.
      )
    )


  cntprocessname <- tblprocessnames |>
    dplyr::distinct() |>
    dplyr::count(.data$processname) |>
    dplyr::filter(.data$n > 1)
  if (nrow(cntprocessname) > 0) {
    stop(paste0(
      "The named processes ",
      paste(cntprocessname$processname, collapse = ", "),
      " occur for different types or labels. Keep them unique."
    ))
  }

  # check if process names are used across processgroup/processlabel
  # dangerous choice and should be discouraged

  # COMPILE PROCESSES
  print("Compiling processes")

  # Travel
  if (nrow(tbltravel) > 0) {
    print("Adding travel")
    tbltravel <- tbltravel |>
      dplyr::mutate(frommetapopulation = replace(
        .data$frommetapopulation,
        .data$frommetapopulation == "",
        list(allmetapops)
      )) |>
      tidyr::unnest("frommetapopulation") |>
      dplyr::mutate(tometapopulation = replace(
        .data$tometapopulation,
        .data$tometapopulation == "",
        list(allmetapops)
      )) |>
      tidyr::unnest("tometapopulation") |>
      dplyr::filter(.data$frommetapopulation != .data$tometapopulation) |>
      dplyr::mutate(travelbasestates = replace(
        .data$travelbasestates,
        .data$travelbasestates == "",
        list(allbasestates)
      )) |>
      dplyr::mutate(travelgroups = replace(
        .data$travelgroups,
        .data$travelgroups == "",
        list(allgroupslist)
      ))
    # grab groups
    tbltravel_grab <- lapply(
      tbltravel$travelgroups,
      grabfromtbl,
      tblupdatestate
    )

    # loop because dplyr seems to not work as expected
    for (loopidx in seq_along(tbltravel_grab)) {
      currtbltravel <- tbltravel[loopidx, ]
      tbltravel_grab_pre <- tbltravel_grab[[loopidx]] |>
        dplyr::filter(.data$basestates %in% tbltravel$travelbasestates[[loopidx]]) # nolint: line_length_linter.

      tbltravel_grab_from <-
        tbltravel_grab_pre |>
        dplyr::filter(.data$metapopulation == currtbltravel$frommetapopulation) |> # nolint: line_length_linter.
        dplyr::rename(
          frommetapopulation = "metapopulation",
          fromupdatedstate = "updatedstate"
        )
      tbltravel_grab_to <-
        tbltravel_grab_pre |>
        dplyr::filter(.data$metapopulation == currtbltravel$tometapopulation) |>
        dplyr::rename(
          tometapopulation = "metapopulation",
          toupdatedstate = "updatedstate"
        )
      bynames <- setdiff(
        colnames(tbltravel_grab_to),
        c("tometapopulation", "toupdatedstate")
      )
      tbltravel_grab_join <- dplyr::full_join(tbltravel_grab_to,
        tbltravel_grab_from,
        by = bynames
      ) |>
        dplyr::left_join(
          currtbltravel |>
            dplyr::select(
              "frommetapopulation",
              "tometapopulation",
              "rate",
              "processname",
              "processgroup"
            ),
          by = c("frommetapopulation", "tometapopulation")
        )

      tbltravel_grab[[loopidx]] <- tbltravel_grab_join
    }

    tblprocess_travel <- dplyr::bind_rows(tbltravel_grab)

    tblprocess_travel <-
      scaleprocesstblbygroups(tblprocess_travel, tblgroups, "scalemigrations")
    # Don't scale by metapopulation because that is ambiguous for travel
  }

  # TRANSITIONS
  # do in loop because easy
  if (nrow(tbltransition_expand) > 0) {
    print("Adding transitions")
    tbltransition_expand_fromgroup <- lapply(
      tbltransition_expand$fromgroups,
      grabfromtbl,
      tblupdatestate
    )
    tbltransition_expand_togroup <- lapply(
      tbltransition_expand$togroups,
      grabfromtbl,
      tblupdatestate
    )
    tbltransition_join <- list()
    for (loopidx in seq_len(nrow(tbltransition_expand))) {
      currtbltransition_expand <- tbltransition_expand[loopidx, ]
      currfromstate <- currtbltransition_expand$fromstate
      currtostate <- currtbltransition_expand$tostate
      currchaincol <-
        paste0("chainid", as.character(currtbltransition_expand$compileid))

      curr_fromgroup <- tbltransition_expand_fromgroup[[loopidx]] |>
        dplyr::filter(.data$basestates == currfromstate)
      curr_togroup <- tbltransition_expand_togroup[[loopidx]] |>
        dplyr::filter(.data$basestates == currtostate)
      # filter fromstates based on chain if chained
      if (currchaincol %in% colnames(tblupdatestate)) {
        curr_fromgroup <- curr_fromgroup |>
          dplyr::filter(
            .data[[currchaincol]] == currtbltransition_expand$fromchain
          )
      }
      # differentiate if end or shaft of fork
      if (currtbltransition_expand$tochain == 1) {
        # ignore chaining and send to first
        curr_togroup <- curr_togroup |>
          dplyr::select(-dplyr::contains("chainid")) |>
          dplyr::distinct(dplyr::across(-"updatedstate"), .keep_all = TRUE)
        shared_names <- intersect(
          colnames(curr_fromgroup),
          colnames(curr_togroup)
        )

        matching_columns <-
          sapply(
            shared_names,
            function(col) {
              identical(
                curr_fromgroup[[col]],
                curr_togroup[[col]]
              )
            }
          )
        bynames <- shared_names[matching_columns]
        # }
      } else {
        curr_togroup <-
          curr_togroup |>
          dplyr::filter(
            .data[[currchaincol]] == currtbltransition_expand$tochain
          )
        bynames <- setdiff(
          colnames(curr_fromgroup),
          c(currchaincol, "updatedstate")
        )
      }
      tbltransition_join[[loopidx]] <-
        dplyr::left_join(curr_fromgroup,
          curr_togroup,
          by = bynames
        ) |>
        dplyr::rename(
          updatedstatefrom = "updatedstate.x",
          updatedstateto = "updatedstate.y"
        ) |>
        dplyr::select("updatedstatefrom", "updatedstateto") |>
        dplyr::filter(!is.na(.data$updatedstateto)) |>
        dplyr::bind_cols(currtbltransition_expand |>
          dplyr::select(
            "percapitastate",
            "rate",
            "processname",
            "processgroup",
            "compileid"
          ))
    }
    tblprocess_transition <- dplyr::bind_rows(tbltransition_join)

    tblprocess_transition <-
      dplyr::left_join(tblprocess_transition,
        tblupdatestate,
        by = c("updatedstatefrom" = "updatedstate")
      )

    # scale by GROUPS
    tblprocess_transition <-
      scaleprocesstblbygroups(
        tblprocess_transition,
        tblgroups,
        "scaletransitions"
      )

    # scale by SPACE
    tblprocess_transition <-
      scaleprocesstblbyspace(
        tblprocess_transition,
        tblspace,
        "spacescaletransitions"
      )
  }

  # Interactions
  tblinteractions <-
    inputpeter$interactions |> dplyr::rename(groupnames = .data$groupname)

  if (nrow(tblinteractions) > 0) {
    print("Adding interactions to dynamics")

    if (!identical(
      sapply(tblinteractions$states_in, length),
      sapply(tblinteractions$states_in, length)
    )) {
      stop("Interaction states_in and states_out have different
      lengths. Only change of states are allowed, not birth/death")
    }
    emptyinlogic <- "" %in% unlist(tblinteractions$states_in)
    emptyoutlogic <- "" %in% unlist(tblinteractions$states_out)
    if (emptyinlogic || emptyoutlogic) {
      stop("Specified states_in or states_out can not be empty string")
    }

    tblinteractions$changestateidx <- mapply(function(x, y) {
      which(x != y)
    }, tblinteractions$states_in, tblinteractions$states_out)
    # split by metapopulation
    tblinteractions <- tblinteractions |>
      replaceglobal("metapopulation", allmetapops) |>
      tidyr::unnest("metapopulation") |>
      dplyr::mutate(groupnames = replace(
        .data$groupnames,
        .data$groupnames == "",
        list(allgroupslist)
      ))

    # grab groups
    tblinteractions_grab <-
      lapply(tblinteractions$groupnames, grabfromtbl, tblupdatestate)
    list2bind <- list()
    for (loopidx in seq_along(tblinteractions_grab)) {
      currtblinteractions <- tblinteractions[loopidx, ]
      currmetapopulation <- currtblinteractions$metapopulation
      currtblinteractions_grab <- tblinteractions_grab[[loopidx]]
      statesin2combo <- lapply(currtblinteractions$states_in[[1]], function(x) {
        currtblinteractions_grab |>
          dplyr::filter(
            .data$basestates == x,
            .data$metapopulation == currmetapopulation
          ) |>
          dplyr::pull("updatedstate")
      })
      combostatesin <- expand.grid(statesin2combo)
      # loop through changestateidx, make 2nd table
      combostatesin_after <- combostatesin
      for (changeidx in currtblinteractions$changestateidx) {
        currcol2change <- combostatesin_after[, changeidx]
        newstate <- currtblinteractions$states_out[[1]][changeidx]
        currleft <- dplyr::tibble(updatedstate = currcol2change) |>
          dplyr::left_join(tblupdatestate, by = "updatedstate") |>
          dplyr::select(-dplyr::contains("chainid")) |>
          dplyr::select(-"basestates")
        currright <- tblupdatestate |>
          dplyr::filter(
            .data$basestates == newstate,
            .data$metapopulation == currmetapopulation
          ) |>
          dplyr::select(-dplyr::contains("chainid")) |>
          dplyr::distinct(dplyr::across(-"updatedstate"), .keep_all = TRUE) |>
          dplyr::select(-"basestates")
        bynames <- setdiff(colnames(currleft), "updatedstate")
        currjoin <- dplyr::left_join(currleft, currright, by = bynames)
        combostatesin_after[, changeidx] <- currjoin$updatedstate.y
      }

      # on og table set rate, then scale by processes/etc
      colnames2loop <- colnames(combostatesin)
      combostatesin_test <-
        combostatesin |>
        dplyr::mutate(
          rate = 1,
          processname = currtblinteractions$processname,
          processgroup = currtblinteractions$processgroup
        )

      # scale by each group
      for (currcolname in colnames2loop) {
        currtbl2fix <-
          combostatesin_test[
            ,
            c(
              currcolname,
              "rate",
              "processname",
              "processgroup"
            )
          ] |>
          dplyr::rename(updatedstate = tidyselect::all_of(currcolname))
        currtbl2fix <-
          dplyr::left_join(currtbl2fix, tblupdatestate, by = "updatedstate")
        currtbl2fix <-
          scaleprocesstblbygroups(currtbl2fix, tblgroups, "scaleinteractions")
        combostatesin_test$rate <- currtbl2fix$rate
      }
      # take harmonic mean
      combostatesin_test$rate <- paste0(
        "((",
        combostatesin_test$rate,
        ")^(1/",
        as.character(length(colnames2loop)),
        "))*",
        currtblinteractions$rate
      )
      # combine, re-list the var cols
      states_in <- combostatesin |>
        dplyr::rowwise() |>
        dplyr::mutate(states_in = list(dplyr::c_across(dplyr::everything()))) |>
        dplyr::select("states_in")
      states_out <- combostatesin_after |>
        dplyr::rowwise() |>
        dplyr::mutate(
          states_out =
            list(dplyr::c_across(dplyr::everything()))
        ) |>
        dplyr::select("states_out")

      tblnowexpand <-
        dplyr::bind_cols(
          states_in,
          states_out,
          currtblinteractions |>
            dplyr::select(
              -"rate",
              -"states_in",
              -"states_out",
              -"groupnames"
            )
        )
      tblnowexpand[, "rate"] <- combostatesin_test$rate
      # scale by metapopulation
      tblnowexpand <-
        scaleprocesstblbyspace(tblnowexpand, tblspace, "spacescaleinteractions")

      list2bind[[loopidx]] <- tblnowexpand
    }
  }
  tblprocess_interactions <- dplyr::bind_rows(list2bind)

  print("Generating output")
  tblprocesses <- tibble::tibble(
    states_up = list(),
    states_down = list(),
    rate = character(),
    processname = character(),
    processgroup = character(),
    processlabel = character(),
    metapopulation_up = character(),
    metapopulation_down = character(),
    percapitastates = list(),
    percapitarate = character(),
    metapopulation_percapita = character()
  )
  # clean interactions
  if (nrow(tblprocess_interactions) > 0) {
    tibble2bind <- list()
    for (loopidx in seq_len(nrow(tblprocess_interactions))) {
      curr_inter <- tblprocess_interactions[loopidx, ]
      currstates_up <-
        as.vector(curr_inter$states_out[[1]][curr_inter$changestateidx])
      currstates_down <-
        as.vector(curr_inter$states_in[[1]][curr_inter$changestateidx])
      currmetapop <- curr_inter$metapopulation

      setdiff_states <- setdiff(
        as.vector(curr_inter$states_in[[1]]),
        currstates_down
      )

      currpercapitastate <-
        list(as.vector(curr_inter$states_in[[1]]))
      currpercapitarate <- curr_inter$rate
      currpercapitarate_clean <- currpercapitarate
      if (curr_inter$normlogic) {
        currtblntotal <- tblntotal |>
          dplyr::filter(.data$metapopulation == currmetapop)
        currpercapitarate <- paste0(
          currpercapitarate,
          "/",
          currtblntotal$totalpop
        )
        currpercapitarate_clean <- paste0(
          currpercapitarate_clean,
          "/",
          currtblntotal$varname
        )
      }
      currrate <- paste0(
        paste(currpercapitastate[[1]], collapse = "*"),
        "*(",
        currpercapitarate,
        ")"
      )

      currpercapitarate_clean <- paste0(
        paste(setdiff_states, collapse = "*"),
        "*(",
        currpercapitarate_clean,
        ")"
      )
      tibble2bind[[loopidx]] <- tibble::tibble(
        states_up = list(currstates_up),
        states_down = list(currstates_down),
        rate = currrate,
        processname = curr_inter$processname,
        processgroup = curr_inter$processgroup,
        processlabel = "interaction",
        metapopulation_up = currmetapop,
        metapopulation_down = currmetapop,
        percapitastates =
          list(as.vector(curr_inter$states_in[[1]])),
        percapitarate = currpercapitarate_clean,
        metapopulation_percapita = currmetapop
      )
    }
    tibbleinteraction <- dplyr::bind_rows(tibble2bind)
    tblprocesses <- dplyr::bind_rows(tblprocesses, tibbleinteraction)
  }

  if (nrow(tbltransition) > 0) {
    tibble2bind <- list()
    for (loopidx in seq_len(nrow(tblprocess_transition))) {
      currtblprocess_transition <- tblprocess_transition[loopidx, ]
      currstates_up <- currtblprocess_transition$updatedstateto
      currstates_down <- currtblprocess_transition$updatedstatefrom
      currmetapop <- currtblprocess_transition$metapopulation

      currpercapitastate <- currtblprocess_transition$updatedstatefrom
      currpercapitarate <- currtblprocess_transition$rate
      currrate <- paste0("(", currpercapitarate, ")*", currpercapitastate)

      tibble2bind[[loopidx]] <- tibble::tibble(
        states_up = list(currstates_up),
        states_down = list(currstates_down),
        rate = currrate,
        processname = currtblprocess_transition$processname,
        processgroup = currtblprocess_transition$processgroup,
        processlabel = "transition",
        metapopulation_up = currmetapop,
        metapopulation_down = currmetapop,
        percapitastates = list(currpercapitastate),
        percapitarate = currpercapitarate,
        metapopulation_percapita = currmetapop
      )
    }
    tibbletransition <- dplyr::bind_rows(tibble2bind)
    tblprocesses <- dplyr::bind_rows(tblprocesses, tibbletransition)
  }

  if (nrow(tbltravel) > 0) {
    tibble2bind <- list()
    for (loopidx in seq_len(nrow(tblprocess_travel))) {
      currtblprocess_travel <- tblprocess_travel[loopidx, ]
      currstates_up <- currtblprocess_travel$toupdatedstate
      currstates_down <- currtblprocess_travel$fromupdatedstate
      currpercapitarate <- currtblprocess_travel$rate
      currrate <- paste0("(", currpercapitarate, ")*", currstates_down)

      tibble2bind[[loopidx]] <- tibble::tibble(
        states_up = list(currstates_up),
        states_down = list(currstates_down),
        rate = currrate,
        processname = currtblprocess_travel$processname,
        processgroup = currtblprocess_travel$processgroup,
        processlabel = "migration",
        metapopulation_up = currtblprocess_travel$tometapopulation,
        metapopulation_down = currtblprocess_travel$frommetapopulation,
        percapitastates = list(currstates_down),
        percapitarate = currpercapitarate,
        metapopulation_percapita = currtblprocess_travel$frommetapopulation
      )
    }
    tibbletravel <- dplyr::bind_rows(tibble2bind)
    tblprocesses <- dplyr::bind_rows(tblprocesses, tibbletravel)
  }


  # scale_process
  # set_process

  # outstructions, instructions
  currpetercolfun <- function(statesdown, statesup) {
    states2petercolumn(statesdown, statesup, tblupdatestate$updatedstate)
  }
  petersenmat <- do.call(
    cbind,
    mapply(currpetercolfun,
      tblprocesses |> dplyr::pull(.data$states_down),
      tblprocesses |> dplyr::pull(.data$states_up),
      USE.NAMES = FALSE
    )
  )

  # change list to vector if all length 1
  tblprocesses_clean <- tblprocesses
  uplogic <- FALSE %in% (sapply(tblprocesses$states_up, length) == 1)
  downlogic <- FALSE %in% (sapply(tblprocesses$states_down, length) == 1)
  if (!uplogic && !downlogic) {
    tblprocesses_clean$states_up <- unlist(tblprocesses$states_up)
    tblprocesses_clean$states_down <- unlist(tblprocesses$states_down)
  }
  tblupdatestate_clean <-
    dplyr::bind_cols(
      tblupdatestate |>
        dplyr::select(
          "updatedstate",
          "basestates",
          "metapopulation"
        ),
      tblupdatestate |>
        dplyr::select(
          -"updatedstate",
          -"basestates",
          -"metapopulation"
        )
    )

  # keep columns with multiple distinct values
  tblupdatestate_clean <-
    tblupdatestate_clean |> dplyr::select_if(function(x) length(unique(x)) > 1)

  peterlist <- list()
  peterlist$updatedstates <- tblupdatestate$updatedstate
  peterlist$petermatrix <- petersenmat
  peterlist$processrates <- tblprocesses$rate
  modelinstructions <- list()
  modelinstructions$tblupdatedstates <- tblupdatestate_clean
  modelinstructions$tblprocesses <- tblprocesses_clean
  modelinstructions$tbltransitions <- tblprocess_transition
  modelinstructions$tblinteractions <- tblprocess_interactions
  modelinstructions$tbltravel <- tblprocess_travel
  modelinstructions$tblspace <- tblspace
  modelinstructions$tblntotal <- tblntotal
  outputlist <- list(
    modelinstructions = modelinstructions,
    modeloutstructions = peterlist
  )
  return(outputlist)
}
