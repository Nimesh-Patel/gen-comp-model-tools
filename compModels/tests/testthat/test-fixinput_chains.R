test_that("fixinput_chains converts vectors", {
  basestates <- c("S", "I", "R")
  modelinstructions <- define_states(basestates) |>
    add_transition("I", "R", "tau", chainlength = 2) |>
    add_infection("I", "S", "I", "beta")
  compiledmodel <- compilemodel(modelinstructions)

  tblwgroupcols <- compiledmodel[["modelinstructions"]][["tblupdatedstates"]]

  fixinput_chains(chains = 1, tblwgroupcols)

  expect_equal(
    fixinput_chains(chains = NA, tblwgroupcols),
    NA
  )

  expect_equal(
    fixinput_chains(chains = 1, tblwgroupcols),
    fixinput_chains(chains = "1", tblwgroupcols)
  )

  expect_equal(
    names(fixinput_chains(chains = 1, tblwgroupcols)),
    "chainid1"
  )

  expect_error(fixinput_chains("not here", tblwgroupcols))
})
