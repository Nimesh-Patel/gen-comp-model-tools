test_that("fixinput_groupnames converts vectors", {
  basestates <- c("S", "I", "R")
  modelinstructions <- define_states(basestates) |>
    add_transition("I", "R", "tau") |>
    add_infection("I", "S", "I", "beta") |>
    add_group(groupnames = c("old", "young"), grouptype = "Age")

  compiledmodel <- compilemodel(modelinstructions)

  tblwgroupcols <- compiledmodel[["modelinstructions"]][["tblupdatedstates"]]

  fixinput_groupnames("old", tblwgroupcols)

  expect_equal(
    fixinput_groupnames(NA, tblwgroupcols),
    NA
  )

  expect_equal(
    fixinput_groupnames(list(Age = "old"), tblwgroupcols),
    fixinput_groupnames("old", tblwgroupcols)
  )

  expect_equal(
    names(fixinput_groupnames("old", tblwgroupcols)),
    "Age"
  )

  expect_error(fixinput_groupnames("not here", tblwgroupcols))
})
