test_that("sci test for SIR extinction dynamics", {
  sir <- define_states(c("S", "I", "R")) |>
    add_transition("I", "R", "tau") |>
    add_infection("I", "S", "I", "beta")
  sirc <- compilemodel(sir)

  nsims <- 1000
  n <- 1000
  pars <- c(beta = 2, tau = 1, N = n)
  x0 <- c(S = n - 1, I = 1, R = 0)
  maxt <- 7.5
  sciout <- wrap_adaptivetau(
    x0,
    sirc,
    rate_func = NULL,
    pars,
    maxt,
    nsims,
    "adaptivetau"
  )

  sciout_t <- unlist(sapply(sciout, function(x) {
    x |>
      dplyr::filter(I == 0) |>
      dplyr::distinct(I, .keep_all = TRUE) |>
      dplyr::pull(time)
  }))

  tibbleplot <- tibble::tibble(
    tdeath = sort(sciout_t),
    alpha = seq_along(sciout_t) / nsims
  )
  # bind columns for comparison
  tibbleplot <-
    dplyr::bind_rows(tibble::tibble(tdeath = 0, alpha = 0), tibbleplot)

  prob <- c(alpha = 0)
  tdependent_extinction <- function(time, prob, pars) {
    with(as.list(c(prob, pars)), {
      dalpha <- (beta * alpha - (1 / tau)) * (alpha - 1)
      list(c(dalpha))
    })
  }

  tode <- seq(0.01, maxt, .01)
  emp_alpha <- unlist(sapply(
    tode,
    function(x) {
      tibbleplot |>
        dplyr::filter(tdeath < x) |>
        tail(n = 1) |>
        dplyr::pull(alpha)
    }
  ))

  out <-
    deSolve::ode(
      y = prob,
      times = tode,
      func = tdependent_extinction,
      parms = pars
    )

  kstest_score <- max(abs(emp_alpha - out[, "alpha"]))
  expect_equal(kstest_score, 0, tolerance = .05)
})
