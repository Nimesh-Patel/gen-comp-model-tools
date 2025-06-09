#' Move popsize for between states matching input features
#'
#' Multi-argument function that flexibly moves populations between states. Users
#' can specify features that change before and after the move, along with
#' features that do not change before and after the move, but specify
#' which populations to move between. Multiple states can move into a single
#' state, allowing a more-to-fewer mapping across the whole popsize table.
#' Checks ensure proper input argument combinations.
#'
#' @param tblpopsize compiled updatedstate table with popsize value
#' @param movepopsize value to make all popsizes that match input features
#' @param basestates_before character vector of basestate to move movepopsize
#' from. Many values may be provided.
#' default is NA which ignores this input
#' @param basestates_after character vector of basestate to move movepopsize
#' into. At most one value may be provided.
#' default is NA which ignores this input
#' @param basestates vector of basestates specifying moving populations.
#' default is NA which designates all basestates, unless user provides _before
#' and _after input arguments
#' @param groupnames_before character vector/list of groupnames to move
#' movepopsize from. Many values may be specified.
#' default is NA which ignores this input
#' @param groupnames_after character vector/list of groupnames to move
#' movepopsize into. Single value must be specified (altogether,
#' or for each name in list)
#' default is NA which ignores this input
#' @param groupnames list of groupnames specifying moving populations.
#' default is NA which designates all basestates, unless user provides _before
#' and _after input arguments
#' @param metapopulation_before character vector of metapopulation to move
#' movepopsize from. Many values may be provided.
#' default is NA which ignores this input
#' @param metapopulation_after character vector of metapopulations to move
#' movepopsize into. At most one value may be provided.
#' default is NA which ignores this input
#' @param metapopulation vector of metapopulations specifying moving
#' populations. default is NA which designates all metapopulation, unless user
#' provides _before and _after input arguments
#' @param chains_before vector/named list of chains steps to move
#' movepopsize from. Many values may be specified.
#' default is NA which ignores this input
#' @param chains_after vector/named list of chains steps to move movepopsize
#' into. Vector input must be length 1. Named values must be length 1 for
#' named list input, but the list itself can have length >1. Names correspond to
#' chainids.
#' default is NA which ignores this input
#' @param chains vector/named list of chains steps to fix when moving
#' populations.
#' default is NA which assumes moving from first step in all chains to
#' first step in other chains as stratified by other input. Users can specify
#' chains _before and _after using the respective input arguments
#' @param movefromallchains flag to specify moving from all steps in before
#' chains. The default is set to false so movements flow from first step in
#' before chains.
#' @return tibble with columns of updated state names
#' (updatedstates) and current conditions (popsize)
#' @export
#' @importFrom rlang .data
movepopsize_byfeature <- function(tblpopsize, # nolint: cyclocomp_linter.
                                  movepopsize,
                                  basestates_before = NA,
                                  basestates_after = NA,
                                  basestates = NA,
                                  groupnames_before = NA,
                                  groupnames_after = NA,
                                  groupnames = NA,
                                  metapopulation_before = NA,
                                  metapopulation_after = NA,
                                  metapopulation = NA,
                                  chains_before = NA,
                                  chains_after = NA,
                                  chains = NA,
                                  movefromallchains = FALSE) {
  check_inputpopsize(tblpopsize)
  # compare popsize type before and after applying
  intlogic_before <- intlogic(tblpopsize[["popsize"]])

  # check acceptable input and id columns to ignore when joining.
  # These columns change values before and after
  col2ignore <- c()

  logic_bstate <- !identical(basestates, NA)
  logic_bstate_before <- !identical(basestates_before, NA)
  logic_bstate_after <- !identical(basestates_after, NA)
  if (logic_bstate_before || logic_bstate_after) {
    if (logic_bstate) {
      stop("Can not specify both basestates and either
              basestates_before or basestates_after.")
    }
    if (!(logic_bstate_before && logic_bstate_after)) {
      stop("Both basestates_before and basestates_after must either be
        specified or NA")
    }
    col2ignore <- c(col2ignore, "basestates")
  }

  if (logic_bstate) {
    basestates_before <- basestates
    basestates_after <- basestates
  }

  logic_meta <- !identical(metapopulation, NA)
  logic_meta_before <- !identical(metapopulation_before, NA)
  logic_meta_after <- !identical(metapopulation_after, NA)
  if (logic_meta_before || logic_meta_after) {
    if (logic_meta) {
      stop("Can not specify both metapopulation and either
              metapopulation_before or metapopulation_after.")
    }
    if (!(logic_meta_before && logic_meta_after)) {
      stop("Both metapopulation_before and metapopulation_after must either be
        specified or NA")
    }
    col2ignore <- c(col2ignore, "metapopulation")
  }

  if (logic_meta) {
    metapopulation_before <- metapopulation
    metapopulation_after <- metapopulation
  }

  # filter popsize table based on metapopulation and basestate
  # states can have
  tblpopsize_before <-
    filterpopsize_byfeature(tblpopsize,
      basestates = basestates_before,
      metapopulation = metapopulation_before
    )
  tblpopsize_after <-
    filterpopsize_byfeature(tblpopsize,
      basestates = basestates_after,
      metapopulation = metapopulation_after
    )

  # groupnames, coerce input into named list if not already
  # groupnames can both change or be a conditioning feature
  # grouptypes can not be same when conditioning and when changing
  # if they change, we force changes within the same type e.g., old to young
  logic_gnames <- !identical(groupnames, NA)
  logic_gnames_before <- !identical(groupnames_before, NA)
  logic_gnames_after <- !identical(groupnames_after, NA)

  if (logic_gnames_before || logic_gnames_after || logic_gnames) {
    # check proper input
    groupnames_before <-
      fixinput_groupnames(groupnames_before, tblpopsize_before)
    groupnames_after <-
      fixinput_groupnames(groupnames_after, tblpopsize_after)
    groupnames_condition_before <-
      fixinput_groupnames(groupnames, tblpopsize_before)
    groupnames_condition_after <-
      fixinput_groupnames(groupnames, tblpopsize_after)
    if (logic_gnames_before || logic_gnames_after) {
      if (!(logic_gnames_before && logic_gnames_after)) {
        stop("Both groupnames_before and groupnames_after must either be
          specified or NA")
      }
      if (!setequal(names(groupnames_before), names(groupnames_after))) {
        stop("groupnames_before and groupnames_after differ in type names. Can
        only move populations between groupvalues in the same type.")
      }
      templogic <-
        names(groupnames_condition_before) %in% names(groupnames_before)
      if (TRUE %in% templogic) {
        stop("groupnames and groupnames_before share the same grouptype. Can not
        condition and vary within the same grouptype.")
      }
      templogic <-
        names(groupnames_condition_after) %in% names(groupnames_after)
      if (TRUE %in% templogic) {
        stop("groupnames and groupnames_after share the same grouptype. Can not
        condition and vary within the same grouptype.")
      }
      col2ignore <- c(col2ignore, names(groupnames_before))
    }

    # both condition and vary by groupname
    bothgnames_before <- unlist(
      list(
        groupnames_before,
        groupnames_condition_before
      ),
      recursive = FALSE
    )
    bothgnames_after <- unlist(
      list(
        groupnames_after,
        groupnames_condition_after
      ),
      recursive = FALSE
    )
    bothgnames_before <- bothgnames_before[!is.na(bothgnames_before)]
    bothgnames_after <- bothgnames_after[!is.na(bothgnames_after)]

    tblpopsize_before <-
      filterpopsize_byfeature(tblpopsize_before,
        groupnames = bothgnames_before
      )
    tblpopsize_after <-
      filterpopsize_byfeature(tblpopsize_after,
        groupnames = bothgnames_after
      )
  }

  # chains, coerce input into named list if not already
  # chains can both change or be a conditioning feature, but not when
  # sharing the same chainid
  # chains can be defined only before or after or both
  # if they change, we force changes within the same type e.g., old to young
  logic_chains <- !identical(chains, NA)
  logic_chains_before <- !identical(chains_before, NA)
  logic_chains_after <- !identical(chains_after, NA)
  # check proper input
  chains_before <- fixinput_chains(chains_before, tblpopsize_before)
  chains_after <- fixinput_chains(chains_after, tblpopsize_after)
  chains_condition_before <-
    fixinput_chains(chains, tblpopsize_before)
  chains_condition_after <-
    fixinput_chains(chains, tblpopsize_after)
  if (logic_chains_before || logic_chains_after || logic_chains) {
    # check consistency for stable chain values
    if (!identical(chains_condition_before, chains_condition_after)) {
      stop("chains argument specifies chain step values to remain constant
      when moving populations. The value does not exist when conditioning on
      the other _before and _after input arguments. If chain exists only
      for before or after features specify using chains_before or chains_after
      respectively.")
    }

    if (logic_chains_before && logic_chains_after) {
      if (!setequal(names(chains_before), names(chains_after))) {
        stop("chains_before and chains _after differ in type names. Can
        only move populations between groupvalues in the same type.")
      }
      templogic <- names(chains_condition_before) %in% names(chains_before)
      if (TRUE %in% templogic) {
        stop("chains and chains_before share the same grouptype. Can not
        condition and vary within the same grouptype.")
      }
      if (TRUE %in% (names(chains_condition_after) %in% names(chains_after))) {
        stop("chains and chains_after share the same grouptype. Can not
        condition and vary within the same grouptype.")
      }
      col2ignore <- c(col2ignore, names(chains_before))
    }

    # both condition and vary by groupname
    bothchains_before <- unlist(
      list(
        chains_before,
        chains_condition_before
      ),
      recursive = FALSE
    )
    bothchains_after <- unlist(
      list(
        chains_after,
        chains_condition_after
      ),
      recursive = FALSE
    )
    bothchains_before <- bothchains_before[!is.na(bothchains_before)]
    bothchains_after <- bothchains_after[!is.na(bothchains_after)]

    if (length(bothchains_before) > 0) {
      tblpopsize_before <-
        filterpopsize_byfeature(tblpopsize_before, chains = bothchains_before)
    }

    if (length(bothchains_after) > 0) {
      tblpopsize_after <-
        filterpopsize_byfeature(tblpopsize_after, chains = bothchains_after)
    }
  }
  # no before or after specified
  if (length(col2ignore) == 0) {
    stop("No before and after feature specfied. Check input")
  }

  # Flow into first step for chains not in before states
  # Flow from first step for chains not in after states
  tblbefore_nona <-
    tblpopsize_before |> dplyr::select(dplyr::where(~ !all(is.na(.))))
  tblafter_nona <-
    tblpopsize_after |> dplyr::select(dplyr::where(~ !all(is.na(.))))
  colnames_before <- colnames(tblbefore_nona)
  colnames_after <- colnames(tblafter_nona)
  chain_before <- colnames_before[grepl("chainid", colnames_before)]
  chain_after <- colnames_after[grepl("chainid", colnames_after)]
  chain_after <- setdiff(chain_after, names(chains_after))
  chain_before <- setdiff(chain_before, names(chains_before))
  newchainafter <- setdiff(chain_after, chain_before)
  if (length(newchainafter) > 0) {
    keeplogic <- rep(TRUE, nrow(tblafter_nona))
    for (currchainafter in newchainafter) {
      keeplogic <- keeplogic & (tblafter_nona[[currchainafter]] == 1)
    }
    # filter
    tblafter_nona <- tblafter_nona[keeplogic, ]
  }
  newchainbefore <- setdiff(chain_before, chain_after)
  if (!movefromallchains) {
    if (length(newchainbefore) > 0) {
      keeplogic <- rep(TRUE, nrow(tblbefore_nona))
      for (currchainbefore in newchainbefore) {
        # with chains==1 and NA chains, important for many to one mapping
        keeplogic <-
          keeplogic & ((tblbefore_nona[[currchainbefore]] == 1) | (is.na(tblbefore_nona[[currchainbefore]]))) # nolint: line_length_linter.
      }
      # filter
      tblbefore_nona <- tblbefore_nona[keeplogic, ]
    }
  }
  # join temporary tables to get updatedstate mapping. Many-to-one logic allowed
  # only join if ambiguous relationship between filtered tables
  if (nrow(tblafter_nona) == 1) {
    tbljoin <-
      tibble::tibble(
        updatedstate_before = tblbefore_nona[["updatedstate"]],
        updatedstate_after = tblafter_nona[["updatedstate"]]
      )
  } else {
    # remove columns that are different, these are determined by before_after
    tblbefore_nona <-
      tblbefore_nona |>
      dplyr::select(-dplyr::one_of(col2ignore)) |>
      dplyr::select(-"popsize") |>
      dplyr::rename(updatedstate_before = "updatedstate") |>
      suppressWarnings()
    tblafter_nona <-
      tblafter_nona |>
      dplyr::select(-newchainafter) |>
      dplyr::select(-dplyr::one_of(col2ignore)) |>
      dplyr::select(-"popsize") |>
      dplyr::rename(updatedstate_after = "updatedstate") |>
      suppressWarnings()

    # double check for many-to-one
    tbljoin <- dplyr::left_join(tblbefore_nona, tblafter_nona)
  }
  # loop to apply movement
  for (rowidx in seq_len(nrow(tbljoin))) {
    namebefore <- tbljoin[["updatedstate_before"]][[rowidx]]
    beforelogic <- tblpopsize[["updatedstate"]] == namebefore
    nameafter <- tbljoin[["updatedstate_after"]][[rowidx]]
    afterlogic <- tblpopsize[["updatedstate"]] == nameafter
    tblpopsize[["popsize"]][beforelogic] <-
      tblpopsize[["popsize"]][beforelogic] - movepopsize
    tblpopsize[["popsize"]][afterlogic] <-
      tblpopsize[["popsize"]][afterlogic] + movepopsize
  }

  if (min(tblpopsize[["popsize"]]) < 0) {
    stop("Negative popsize after moving populations. Check input.")
  }
  tblpopsize <- maintainpopsizetype(tblpopsize, "popsize")
  intlogic_after <- intlogic(tblpopsize[["popsize"]])
  if (intlogic_before != intlogic_after) {
    message("popsize changed between non-integer and integer. Ensure this is
    intentional.")
  }

  tblpopsize
}
