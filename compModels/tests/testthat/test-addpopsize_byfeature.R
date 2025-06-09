test_that("addpopsize_byfeature sets initial conditions", {
  basestates <- c("S", "I", "R")
  groups <- c("group1", "group2")
  metapops <- c("UK", "USA")
  modelinstructions <- define_states(basestates) |>
    add_transition("I", "R", "tau") |>
    add_infection("I", "S", "I", "beta") |>
    add_group(groups, grouptype = "type1") |>
    define_metapopulations(metapops)
  compiledmodel <- compilemodel(modelinstructions)
  tblpopsize <- define_popsize(compiledmodel)

  currpop <- .2
  check1 <- tblpopsize |> addpopsize_byfeature(currpop,
    basestates = "S"
  )

  check2 <- tblpopsize |> addpopsize_byfeature(currpop,
    basestates = "S",
    groupnames = "group2"
  )

  check3 <- tblpopsize |> addpopsize_byfeature(currpop,
    basestates = "S",
    groupnames = "group2",
    metapopulation = "UK"
  )

  expect_equal(
    sum(check1[["popsize"]]),
    length(groups) * length(metapops) * currpop
  )

  expect_equal(
    sum(check2[["popsize"]]),
    length(metapops) * currpop
  )

  expect_equal(
    sum(check3[["popsize"]]),
    currpop
  )
})
