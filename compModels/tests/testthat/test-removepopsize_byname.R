test_that("removepopsize_byname removes appropriately", {
  basestates <- c("S", "I", "R")
  modelinstructions <- define_states(basestates) |>
    add_transition("I", "R", "tau") |>
    add_infection("I", "S", "I", "beta") |>
    add_group(c("group1", "group2"), grouptype = "type1")
  compiledmodel <- compilemodel(modelinstructions)
  tblpopsize <- define_popsize(compiledmodel)
  setvec <- namevec <- c("S_type1group1" = 1, "I_type1group1" = 1.5)
  namevec <- c("S_type1group1" = .5, "I_type1group1" = 1)
  tblpopsize <-
    tblpopsize |>
    setpopsize_byname(setvec) |>
    removepopsize_byname(namevec)

  expect_equal(
    tblpopsize |>
      dplyr::filter(updatedstate == "S_type1group1") |>
      dplyr::pull("popsize"),
    tblpopsize |>
      dplyr::filter(updatedstate == "I_type1group1") |>
      dplyr::pull("popsize")
  )
})
