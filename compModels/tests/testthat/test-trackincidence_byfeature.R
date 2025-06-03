test_that("trackincidence_byfeature generates incidence", {
  basestates <- c("S", "I", "R")
  modelinstructions <- define_states(basestates) |>
    add_transition("I", "R", "tau") |>
    add_infection("I", "S", "I", "beta")
  compiledmodel <- compilemodel(modelinstructions)
  compiledmodel_byname <- compiledmodel |>
    trackincidence_byname("S", trackname = "incS") |>
    trackincidence_byname("I", trackname = "incI") |>
    trackincidence_byname("R", trackname = "incR")

  compiledmodel_winc <- compiledmodel |>
    trackincidence_byfeature(basestates = "S", trackname = "incS") |>
    trackincidence_byfeature(basestates = "I", trackname = "incI") |>
    trackincidence_byfeature(basestates = "R", trackname = "incR")

  expect_equal(compiledmodel_byname, compiledmodel_winc)

  # same processes
  modelinstructions_complex <- modelinstructions |>
    define_metapopulations(c("UK", "USA"))

  compiledmodel_complex <- compilemodel(modelinstructions_complex)

  winc_all <- compiledmodel_complex |>
    trackincidence_byfeature(basestates = "I", trackname = "incI")

  winc_single <- compiledmodel_complex |>
    trackincidence_byfeature(
      basestates = "I",
      metapopulation = "USA",
      trackname = "incI"
    )

  expect_equal(
    sum(winc_all[["modeloutstructions"]][["petermatrix"]][7, ]),
    2 * sum(winc_single[["modeloutstructions"]][["petermatrix"]][7, ])
  )
})
