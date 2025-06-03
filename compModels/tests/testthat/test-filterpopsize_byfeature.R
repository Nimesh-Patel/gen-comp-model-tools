test_that("filterpopsize_byfeature filters tables", {
  basestates <- c("S", "I", "R")
  modelinstructions <- define_states(basestates) |>
    add_transition("I", "R", "tau", chainlength = 2) |>
    add_infection("I", "S", "I", "beta") |>
    add_group(c("group1", "group2"), grouptype = "type1") |>
    define_metapopulations(c("UK", "USA"))
  compiledmodel <- compilemodel(modelinstructions)
  tblpopsize <- compiledmodel[["modelinstructions"]][["tblupdatedstates"]]

  single <- filterpopsize_byfeature(tblpopsize,
    basestates = "I",
    groupnames = "group2",
    metapopulation = "UK",
    chains = 2
  )

  expect_equal(nrow(single), 1)

  eq1_1 <- filterpopsize_byfeature(tblpopsize,
    basestates = NA,
    groupnames = "group1",
    metapopulation = NA,
    chains = NA
  )
  eq1_2 <- filterpopsize_byfeature(tblpopsize,
    basestates = NA,
    groupnames = list(type1 = "group1"),
    metapopulation = NA,
    chains = NA
  )

  expect_equal(eq1_1, eq1_2)

  eq2_1 <- filterpopsize_byfeature(tblpopsize,
    basestates = NA,
    groupnames = NA,
    metapopulation = NA,
    chains = 1
  )
  eq2_2 <- filterpopsize_byfeature(tblpopsize,
    basestates = NA,
    groupnames = NA,
    metapopulation = NA,
    chains = list(chainid1 = 1)
  )

  expect_equal(eq2_1, eq2_2)

  eq3_1 <- filterpopsize_byfeature(tblpopsize,
    basestates = "S",
    groupnames = NA,
    metapopulation = NA,
    chains = 1
  )

  eq3_2 <- filterpopsize_byfeature(tblpopsize,
    basestates = "S",
    groupnames = NA,
    metapopulation = NA,
    chains = NA
  )
  expect_equal(eq3_1, eq3_2)
})
