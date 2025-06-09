test_that("movepopsize_byfeature sets initial conditions", {
  basestates <- c("S", "I", "R")
  type1vec <- c("group1", "group2")
  type2vec <- c("group3", "group4")
  metapops <- c("UK", "USA")
  chainlength <- 2
  modelinstructions <- define_states(basestates) |>
    add_transition("I", "R", "tau", chainlength = chainlength) |>
    add_infection("I", "S", "I", "beta") |>
    add_group(type1vec, grouptype = "type1") |>
    add_group(type2vec, grouptype = "type2") |>
    combine_groups() |>
    define_metapopulations(metapops)

  compiledmodel <- compilemodel(modelinstructions)

  tblpopsize <- define_popsize(compiledmodel) |>
    setpopsize_byfeature(1, basestates = "S")

  check1 <- tblpopsize |>
    movepopsize_byfeature(.2,
      basestates_before = "S",
      basestates_after = "I"
    )
  # flow into first chain

  check1same <- tblpopsize |>
    movepopsize_byfeature(.2,
      basestates_before = "S",
      basestates_after = "I",
      chains_after = 1
    )

  expect_equal(check1, check1same)

  expect_error(tblpopsize |> movepopsize_byfeature(.2, basestates_before = "S", basestates_after = "I", chains = 1)) # nolint: line_length_linter.


  checkgroup <- tblpopsize |>
    movepopsize_byfeature(.2,
      basestates_before = "S",
      basestates_after = "I",
      groupnames = "group1"
    )

  checkgroupsame <- tblpopsize |>
    movepopsize_byfeature(.2,
      basestates_before = "S",
      basestates_after = "I",
      groupnames = list(type1 = "group1")
    )

  # moving between types
  expect_error(tblpopsize |> movepopsize_byfeature(.2, basestates_before = "S", basestates_after = "I", groupnames_before = "group1", groupnames_after = "group4")) # nolint: line_length_linter.
  checkn <- tblpopsize |> movepopsize_byfeature(.2, basestates_before = "S", basestates_after = "I", groupnames = "group1", groupnames_before = "group3", groupnames_after = "group4") # nolint: line_length_linter.

  tblpopsize |>
    dplyr::filter(type1 == "group2", basestates == "S") |>
    dplyr::distinct(.data$popsize) |>
    nrow()
  expect_equal(
    tblpopsize |>
      dplyr::filter(type1 == "group2", basestates == "S") |>
      dplyr::distinct(.data$popsize) |>
      dplyr::pull("popsize"),
    1
  )

  expect_equal(
    checkn |>
      dplyr::filter(
        type2 == "group3",
        basestates == "S",
        type1 == "group1"
      ) |>
      dplyr::distinct(.data$popsize) |>
      dplyr::pull(),
    .8
  )

  expect_equal(
    checkn |>
      dplyr::filter(
        type2 == "group4",
        basestates == "I",
        type1 == "group1",
        chainid1 == 1
      ) |>
      dplyr::distinct(.data$popsize) |>
      dplyr::pull(),
    .2
  )

  expect_equal(
    checkn |>
      dplyr::filter(
        type2 == "group4",
        basestates == "I",
        type1 == "group2",
        chainid1 == 1
      ) |>
      dplyr::distinct(.data$popsize) |>
      dplyr::pull(),
    0
  )

  expect_equal(
    checkn |>
      dplyr::filter(type1 == "group2", basestates == "S") |>
      dplyr::distinct(.data$popsize) |>
      dplyr::pull("popsize"),
    1
  )

  expect_equal(sum(checkn$popsize), sum(tblpopsize$popsize))
  # check chaining
  basestates <- c("S", "E", "I", "R")
  chainlengthe <- 2
  chainlengthi <- 2
  modelinstructions <- define_states(basestates) |>
    add_transition("I", "R", "tau", chainlength = chainlengthi) |>
    add_transition("E", "I", "taue", chainlength = chainlengthe) |>
    add_infection("I", "S", "I", "beta")

  compiledmodel <- compilemodel(modelinstructions)

  tblpopsize <- define_popsize(compiledmodel) |>
    setpopsize_byfeature(1, basestates = "S")

  expect_error(tblpopsize |> movepopsize_byfeature(.2, basestates_before = "E", basestates_after = "I", chains = 2)) # nolint: line_length_linter.

  checkchain <-
    tblpopsize |>
    movepopsize_byfeature(.4,
      basestates_before = "S",
      basestates_after = "E",
      chains_after = 2
    ) |>
    movepopsize_byfeature(.1,
      basestates_before = "E",
      basestates_after = "I",
      chains_before = 2
    )

  tblpopsize |>
    movepopsize_byfeature(.4,
      basestates_before = "S",
      basestates_after = "E",
      chains_after = 1
    ) |>
    movepopsize_byfeature(.1,
      basestates_before = "E",
      basestates_after = "I"
    )

  checkmultiple <- tblpopsize |>
    movepopsize_byfeature(.4,
      basestates_before = "S",
      basestates_after = "E"
    ) |>
    movepopsize_byfeature(.1,
      basestates_before = c("S", "E"),
      basestates_after = "R"
    )

  temp_prep <-
    tblpopsize |>
    movepopsize_byfeature(.4, basestates_before = "S", basestates_after = "E")
  tem <- c("S", "E")
  expect_error(temp_prep |> movepopsize_byfeature(.1, basestates_before = tem, basestates_after = "R", movefromallchains = TRUE)) # nolint: line_length_linter.

  expect_equal(
    checkchain |> dplyr::filter(basestates == "S") |> dplyr::pull(popsize),
    .6
  )
  expect_equal(
    checkchain |>
      dplyr::filter(basestates == "E", chainid2 == 2) |>
      dplyr::pull(popsize),
    .3
  )
  expect_equal(
    checkchain |>
      dplyr::filter(basestates == "I", chainid1 == 1) |>
      dplyr::pull(popsize),
    .1
  )
})
