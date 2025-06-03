test_that("trackincidence_byname generates incidence", {
  basestates <- c("S", "I", "R")
  modelinstructions <- define_states(basestates) |>
    add_transition("I", "R", "tau") |>
    add_infection("I", "S", "I", "beta")
  compiledmodel <- compilemodel(modelinstructions)
  compiledmodel_winc <- compiledmodel |>
    trackincidence_byname("S", trackname = "incS") |>
    trackincidence_byname("I", trackname = "incI") |>
    trackincidence_byname("R", trackname = "incR")

  # same processes
  expect_equal(
    compiledmodel[["modeloutstructions"]][["processrates"]],
    compiledmodel_winc[["modeloutstructions"]][["processrates"]]
  )

  newupdatedstates <-
    compiledmodel_winc[["modeloutstructions"]][["updatedstates"]]
  oldupdatedstates <-
    compiledmodel[["modeloutstructions"]][["updatedstates"]]
  expect_equal(
    length(newupdatedstates),
    length(oldupdatedstates) + length(oldupdatedstates)
  )
})
