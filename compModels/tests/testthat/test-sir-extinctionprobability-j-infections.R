test_that("calculating sir extinction probability with exactly j infections", {
  sir <- define_states(c("s", "i", "r")) |>
    add_transition("i", "r", "tau") |>
    add_infection("i", "s", "i", "beta")
  sirc <- compilemodel(sir)

  nsims <- 1000
  n <- 1000
  pars <- c(beta = 2, tau = 1, n = n)
  x0 <- c(s = n - 1, i = 1, r = 0)
  maxt <- 25
  print("running stochastic SIR simulation")
  set.seed(389)
  sciout <- wrap_adaptivetau(
    x0,
    sirc,
    rate_func = NULL,
    pars,
    maxt,
    nsims,
    "adaptivetau"
  )

  print("comparing final size probability & probability of j infections")
  exact_finalsize_prob <- function(final_infected, pars) {
    with(as.list(c(pars)), {
      combination_term <- choose(2 * final_infected - 2, final_infected - 1)
      omega_j <- (1 / final_infected) *
        ((pars[[1]]^(final_infected - 1)) * (
          (1 / pars[[2]])^final_infected
        )) /
        ((pars[[1]] + (1 / pars[[2]]))^((2 * final_infected) - 1)) *
        combination_term
      omega_j
    })
  }
  sir_pars <- c(beta = 2, tau = 1)
  exact_finalsize_prob(1, sir_pars)

  last_r <- sort(unlist(lapply(sciout, function(x) tail(x$r, n = 1))))
  table(last_r) / 1000
  probs <- as.numeric(table(last_r) / 1000)
  prob_tbl <- data.frame(seq = unique(last_r), prob = probs)
  blank_tbl <- data.frame(seq = 1:1000, prob = 0)

  final_tbl <- blank_tbl |>
    dplyr::left_join(prob_tbl, by = "seq", suffix = c("", ".new")) |>
    dplyr::mutate(prob = ifelse(!is.na(prob.new), prob.new, prob)) |>
    dplyr::select(seq, prob)

  get_prob_for_j <- function(seq_value, tbl) {
    result <- tbl$prob[tbl$seq == seq_value]
    if (length(result) == 0) {
      warning("seq_value not found in the table.")
      return(NA)
    }
    result
  }

  comparison <- abs(
    exact_finalsize_prob(1, sir_pars) - get_prob_for_j(1, final_tbl)
  )
  expect_true(comparison <= 0.03)
})
