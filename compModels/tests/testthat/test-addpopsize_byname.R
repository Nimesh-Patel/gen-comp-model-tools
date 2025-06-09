test_that("addpopsize_byname adds appropriately", {
  basestates <- c("S", "I", "R")
  modelinstructions <- define_states(basestates) |>
    add_transition("I", "R", "tau") |>
    add_infection("I", "S", "I", "beta") |>
    add_group(c("group1", "group2"), grouptype = "type1")
  compiledmodel <- compilemodel(modelinstructions)
  tblpopsize <- define_popsize(compiledmodel)
  namevec <- c("S_type1group1" = .99, "I_type1group1" = .01)
  tblpopsize <- tblpopsize |> addpopsize_byname(namevec)

  expect_equal(sum(tblpopsize[["popsize"]]), sum(namevec))
})
