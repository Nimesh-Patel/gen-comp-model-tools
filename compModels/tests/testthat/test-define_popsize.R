test_that("define_popsize initiates population sizes", {
  basestates <- c("S", "I", "R")
  modelinstructions <- define_states(basestates) |>
    add_transition("I", "R", "tau") |>
    add_infection("I", "S", "I", "beta")
  compiledmodel <- compilemodel(modelinstructions)

  tblpopsize <- define_popsize(compiledmodel)
  expect_equal(sum(tblpopsize$popsize), 0)

  namedinput_int <- c("S" = 100, "I" = 400)
  tblpopsize_int <- define_popsize(compiledmodel, inputpops = namedinput_int)
  expect_equal(is.integer(tblpopsize_int$popsize), TRUE)
  expect_equal(sum(tblpopsize_int$popsize), sum(namedinput_int))

  namedinput_flt <- c("S" = .99, "I" = .01)
  tblpopsize_flt <- define_popsize(compiledmodel, inputpops = namedinput_flt)
  expect_equal(is.integer(tblpopsize_flt$popsize), FALSE)
  expect_equal(sum(tblpopsize_flt$popsize), sum(namedinput_flt))
})
