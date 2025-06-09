test_that("removepopsize_byfeature sets initial conditions", {
  basestates <- c("S", "I", "R")
  groups <- c("group1", "group2")
  metapops <- c("UK", "USA")
  modelinstructions <- define_states(basestates) |>
    add_transition("I", "R", "tau") |>
    add_infection("I", "S", "I", "beta") |>
    add_group(groups, grouptype = "type1") |>
    define_metapopulations(metapops)
  compiledmodel <- compilemodel(modelinstructions)

  currnames <- compiledmodel$modeloutstructions$updatedstates
  maxpop <- 100
  totx0 <- maxpop * length(currnames)
  x0vec <- rep(maxpop, length(currnames))
  names(x0vec) <- currnames
  tblpopsize <- define_popsize(compiledmodel, inputpops = x0vec)

  currpop <- 25
  check1 <- tblpopsize |> removepopsize_byfeature(currpop,
    basestates = "S"
  )

  check2 <- tblpopsize |> removepopsize_byfeature(currpop,
    basestates = "S",
    groupnames = "group2"
  )

  check3 <- tblpopsize |> removepopsize_byfeature(currpop,
    basestates = "S",
    groupnames = "group2",
    metapopulation = "UK"
  )

  expect_equal(
    sum(check1[["popsize"]]),
    totx0 - length(groups) * length(metapops) * currpop
  )

  expect_equal(
    sum(check2[["popsize"]]),
    totx0 - length(metapops) * currpop
  )

  expect_equal(
    sum(check3[["popsize"]]),
    totx0 - currpop
  )
})
