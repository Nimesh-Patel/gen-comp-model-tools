test_that("movepopsize_byname moves appropriately", {
  basestates <- c("S", "I", "R")
  modelinstructions <- define_states(basestates) |>
    add_transition("I", "R", "tau") |>
    add_infection("I", "S", "I", "beta") |>
    add_group(c("group1", "group2"), grouptype = "type1")
  compiledmodel <- compilemodel(modelinstructions)

  beforename <- "S_type1group1"
  aftername <- "I_type1group1"
  x0 <- 1
  names(x0) <- beforename
  x0move <- .2
  tblpopsize <- define_popsize(compiledmodel, inputpops = x0)
  tblpopsizemove <-
    tblpopsize |>
    movepopsize_byname(x0move, beforename, aftername)

  tblpopsizesame <- tblpopsize |>
    movepopsize_byname(x0move, beforename, aftername) |>
    movepopsize_byname(x0move, aftername, beforename)

  expect_equal(sum(tblpopsize[["popsize"]]), sum(tblpopsizemove[["popsize"]]))
  expect_equal(tblpopsize, tblpopsizesame)
})
