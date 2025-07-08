# Draft model as a function to replicate the modeling work done during the
# Chicago Measles 2024 response. Creators: Catherine Herzog & Bradford Taylor.

# https://www.sciencedirect.com/science/article/pii/S2468266723001305?via%3Dihub#sec1 # nolint
# https://www.cdc.gov/mmwr/volumes/73/wr/mm7319a2.htm
# https://github.com/CDCgov/measles-model-chicago-2024/tree/main

devtools::load_all("compModels")

# Helper functions for dates from Chicago Measles Response 2024
date_to_day <- function(dt, start_date) {
  stopifnot(lubridate::is.Date(dt))
  (dt - start_date) / lubridate::ddays(1) + 1
}
replace_null <- function(lst) {
  purrr::map(lst, function(elt) {
    if (is.null(elt)) {
      NA
    } else {
      elt
    }
  }) |>
    unlist()
}
read_params <- function(path) {
  x <- yaml::read_yaml(path)
  # int_day can have NULLs; replace with NAs before parsing dates
  x$params$int_day <- replace_null(x$params$int_day)
  x$history$start_date <- lubridate::ymd(x$history$start_date)
  str_to_day <- function(s) {
    date_to_day(lubridate::ymd(s), start_date = x$history$start_date)
  }
  x$history <- purrr::modify_at(
    x$history,
    c("import_days", "vaccine_days"),
    str_to_day
  )
  stopifnot(length(x$history$vaccine_days) == length(x$history$vaccine_doses))
  x$params <- purrr::modify_at(x$params, "int_day", str_to_day)
  x
}

# Replicating Chicago Measles Response 2024: same populations & parameters
# Overview:

# Simulation start date 2024-02-01, sim length 365 days, 10,000 sims/scenario
# Parameters
# taue (latent period) = 8 days (mean)
# tau (infectious period) = 5 days (mean),
# r_0 = 25  (recall: beta/gamma = r0; also: beta*tau = r0 = 25)

# Population sizes
# https://github.com/CDCgov/measles-model-chicago-2024/blob/c857905313192741ebff2c61d9a46ca57b808a42/seirMeasles/inst/params.yaml#L6 #nolint
# 1877 total people

# Susceptible population sizes
# https://github.com/CDCgov/measles-model-chicago-2024/blob/c857905313192741ebff2c61d9a46ca57b808a42/seirMeasles/inst/params.yaml#L13 #nolint
# https://github.com/CDCgov/measles-model-chicago-2024/blob/c857905313192741ebff2c61d9a46ca57b808a42/seirMeasles/R/simulate.R#L171 #nolint

# Introduce 1 infected into age group 3 on 2024-02-22 (22 timesteps)
# https://github.com/CDCgov/measles-model-chicago-2024/blob/c857905313192741ebff2c61d9a46ca57b808a42/seirMeasles/R/simulate.R#L26 #nolint
# Interventions:
# 1. Vaccination
#    277, 276, & 329 doses on 2024-03-08, 2024-03-09, 2024-03-10
#    (15, 16, 17 days later)
#    Vaccine efficacy: theta_infant (group 2) and theta_adult (groups 3 & 4)
#    Existing documented_vax_coverage: 0.49
# 2. Active case finding: started 2024-03-08 with efficacy 0.25

# Not included here:
# Tracking compartments: cases, pre-reported cases, reported cases
# Case ascertainment_delay: 2.5 days (mean)
# https://github.com/CDCgov/measles-model-chicago-2024/blob/c857905313192741ebff2c61d9a46ca57b808a42/seirMeasles/R/rates.R#L65C7-L65C55 #nolint
# https://github.com/CDCgov/measles-model-chicago-2024/blob/c857905313192741ebff2c61d9a46ca57b808a42/seirMeasles/R/simulate.R#L345 #nolint

run_measles_chicago_2024 <- function(config_path) {
  base_states <- c("S", "E", "I", "R", "Sv")
  measlesmodel <- define_states(base_states) |>
    add_infection("I", "S", "E", "beta") |>
    add_transition("I", "R", meantime = "tau*(1-int_eff)") |>
    add_transition("E", "I", meantime = "taue", chainlength = 2) |>
    add_infection("I", "Sv", "E", "beta") |>
    add_group(groupnames = as.character(1:5))

  suppressMessages(invisible(capture.output({
    measlescompiled <- compilemodel(measlesmodel) |>
      trackincidence_byfeature(basestates = "I", trackname = "incI_cum_latent")
  })))

  measles_rates <- measlescompiled$modeloutstructions$processrates
  measles_peter <- measlescompiled$modeloutstructions$petermatrix
  measles_states <- measlescompiled$modeloutstructions$updatedstates

  config <- read_params(config_path)

  pop_sizes <- setNames(
    c(
      config$params$pop_1, config$params$pop_2, config$params$pop_3,
      config$params$pop_4, config$params$pop_5
    ),
    measles_states[grep("S_", measles_states)]
  )
  pop_sizes_s <- floor(pop_sizes *
    c(
      config$params$prop_s_1, config$params$prop_s_2,
      config$params$prop_s_3, config$params$prop_s_4,
      config$params$prop_s_5
    ))
  pop_sizes_r <- setNames(
    (pop_sizes - pop_sizes_s),
    measles_states[grep("R_", measles_states)]
  )
  initial_pop <- c(pop_sizes_s, pop_sizes_r)
  tblpopsize <- suppressMessages(define_popsize(measlescompiled,
    inputpops = initial_pop
  ))
  vecpopsize <- setNames(tblpopsize$popsize, tblpopsize$updatedstate)

  # First simulation: start to infectious import
  # int_eff = 0 (case finding intervention has not yet started)
  pre_import <- wrap_adaptivetau(
    vecpopsize,
    measlescompiled,
    rate_func = NULL,
    c(
      beta = config$params$r_0 / config$params$infectious_period,
      tau = config$params$infectious_period, int_eff = 0,
      taue = config$params$latent_period
    ),
    config$history$import_days,
    1,
    "adaptivetau"
  )[[1]]

  lr_preimport <- tail(pre_import, n = 1) |>
    dplyr::select(-1) |> # remove time column
    unlist()
  lr_preimport

  # Importing 1 infection (note: this makes total population size + 1)
  inf_import <- 1 # 1 infected in group 3
  tblpopsize <- define_popsize(measlescompiled, inputpops = lr_preimport) |>
    addpopsize_byfeature(inf_import, basestates = "I", groupnames = "3")
  vecpopsize <- setNames(tblpopsize$popsize, tblpopsize$updatedstate)
  # vecpopsize[grep("I_", names(vecpopsize))] # checking infection #nolint

  # Second simulation: infectious import to first vaccination day
  # int_eff = 0 (case finding intervention has not yet started)
  import_prevacc <- wrap_adaptivetau(
    vecpopsize,
    measlescompiled,
    rate_func = NULL,
    c(
      beta = config$params$r_0 / config$params$infectious_period,
      tau = config$params$infectious_period, int_eff = 0,
      taue = config$params$latent_period
    ),
    (
      (
        config$history$vaccine_days[1] - config$history$import_days[1]
      ) +
        config$history$immunity_onset_delay
    ), # adding on immunity onset delay
    1,
    "adaptivetau"
  )[[1]]

  lr_postimport_prevacc1 <- tail(import_prevacc, n = 1) |>
    dplyr::select(-1) |>
    unlist()
  lr_postimport_prevacc1


  # Vaccination Day 1
  # Case Ascertainment intervention starts (int_eff = 0.25)
  # https://github.com/CDCgov/measles-model-chicago-2024/blob/c857905313192741ebff2c61d9a46ca57b808a42/seirMeasles/R/simulate.R#L82 #nolint
  doses_day1 <- config$history$vaccine_doses[1]
  theta_infant <- config$params$ve_infant
  theta_adult <- config$params$ve
  prop_s <- sum(pop_sizes_s[2:4]) / sum(pop_sizes[2:4])
  vacc_coverage <- config$params$documented_vax_coverage
  prop_vax_to_se <- prop_s / (1.0 - vacc_coverage)

  # Get S and E populations and distribute vaccines to specific compartments
  se <- lr_postimport_prevacc1[grep("^(S_|E_)", names(lr_postimport_prevacc1))]
  se_doses_total <- doses_day1 * prop_vax_to_se
  se_doses <- round(pmin(se_doses_total * se / sum(se), se))

  tblpopsize <- define_popsize(measlescompiled,
    inputpops = lr_postimport_prevacc1
  ) |>
    # assume all E vaccinations are failures so not coding a move to R
    movepopsize_byname(round(se_doses["S_dummygroup1X2"] * theta_infant),
      namebefore = "S_dummygroup1X2",
      nameafter = "R_dummygroup1X2"
    ) |>
    movepopsize_byname(round(se_doses["S_dummygroup1X3"] * theta_adult),
      namebefore = "S_dummygroup1X3",
      nameafter = "R_dummygroup1X3"
    ) |>
    movepopsize_byname(round(se_doses["S_dummygroup1X4"] * theta_adult),
      namebefore = "S_dummygroup1X4",
      nameafter = "R_dummygroup1X4"
    ) |>
    movepopsize_byname(round(se_doses["S_dummygroup1X2"] * (1 - theta_infant)),
      namebefore = "S_dummygroup1X2",
      nameafter = "Sv_dummygroup1X2"
    ) |>
    movepopsize_byname(round(se_doses["S_dummygroup1X3"] * (1 - theta_adult)),
      namebefore = "S_dummygroup1X3",
      nameafter = "Sv_dummygroup1X3"
    ) |>
    movepopsize_byname(round(se_doses["S_dummygroup1X4"] * (1 - theta_adult)),
      namebefore = "S_dummygroup1X4",
      nameafter = "Sv_dummygroup1X4"
    )

  vecpopsize <- setNames(tblpopsize$popsize, tblpopsize$updatedstate)

  # Third simulation: first vaccination day to second vaccination day
  # immunity onset delay already added in previous sim
  postimport_vacc1 <- wrap_adaptivetau(
    vecpopsize,
    measlescompiled,
    rate_func = NULL,
    c(
      beta = config$params$r_0 / config$params$infectious_period,
      tau = config$params$infectious_period, int_eff = config$params$int_eff,
      taue = config$params$latent_period
    ),
    (config$history$vaccine_days[2] - config$history$vaccine_days[1]),
    1,
    "adaptivetau"
  )[[1]]

  lr_postimport_postvacc1 <- tail(postimport_vacc1, n = 1) |>
    dplyr::select(-1) |>
    unlist()
  lr_postimport_postvacc1

  # Vaccination Day 2
  # Case Ascertainment intervention continues (int_eff = 0.25)
  doses_day2 <- config$history$vaccine_doses[2]
  se <- lr_postimport_postvacc1[grep(
    "^(S_|E_)",
    names(lr_postimport_postvacc1)
  )]
  se_doses_total <- doses_day2 * prop_vax_to_se
  se_doses <- round(pmin(se_doses_total * se / sum(se), se))

  tblpopsize <- define_popsize(measlescompiled,
    inputpops = lr_postimport_postvacc1
  ) |>
    # assume all E vaccinations are failures so not coding a move to R
    movepopsize_byname(round(se_doses["S_dummygroup1X2"] * theta_infant),
      namebefore = "S_dummygroup1X2",
      nameafter = "R_dummygroup1X2"
    ) |>
    movepopsize_byname(round(se_doses["S_dummygroup1X3"] * theta_adult),
      namebefore = "S_dummygroup1X3",
      nameafter = "R_dummygroup1X3"
    ) |>
    movepopsize_byname(round(se_doses["S_dummygroup1X4"] * theta_adult),
      namebefore = "S_dummygroup1X4",
      nameafter = "R_dummygroup1X4"
    ) |>
    movepopsize_byname(round(se_doses["S_dummygroup1X2"] * (1 - theta_infant)),
      namebefore = "S_dummygroup1X2",
      nameafter = "Sv_dummygroup1X2"
    ) |>
    movepopsize_byname(round(se_doses["S_dummygroup1X3"] * (1 - theta_adult)),
      namebefore = "S_dummygroup1X3",
      nameafter = "Sv_dummygroup1X3"
    ) |>
    movepopsize_byname(round(se_doses["S_dummygroup1X4"] * (1 - theta_adult)),
      namebefore = "S_dummygroup1X4",
      nameafter = "Sv_dummygroup1X4"
    )
  # lr_postimport_postvacc1[grep("^(S_|Sv_|R)", names(lr_postimport_postvacc1))] #nolint
  vecpopsize <- setNames(tblpopsize$popsize, tblpopsize$updatedstate)

  # Fourth simulation: second vaccination day to third vaccination day
  # immunity onset delay already added in previous sim
  postimport_vacc2 <- wrap_adaptivetau(
    vecpopsize,
    measlescompiled,
    rate_func = NULL,
    c(
      beta = config$params$r_0 / config$params$infectious_period,
      tau = config$params$infectious_period, int_eff = config$params$int_eff,
      taue = config$params$latent_period
    ),
    (config$history$vaccine_days[3] - config$history$vaccine_days[2]),
    1,
    "adaptivetau"
  )[[1]]

  lr_postimport_postvacc2 <- tail(postimport_vacc2, n = 1) |>
    dplyr::select(-1) |>
    unlist()
  lr_postimport_postvacc2

  # Vaccination Day 3
  # Case Ascertainment intervention continues (int_eff = 0.25)
  doses_day3 <- config$history$vaccine_doses[3]
  se <- lr_postimport_postvacc2[grep(
    "^(S_|E_)",
    names(lr_postimport_postvacc2)
  )]
  se_doses_total <- doses_day3 * prop_vax_to_se
  se_doses <- round(pmin(se_doses_total * se / sum(se), se))

  tblpopsize <- define_popsize(measlescompiled,
    inputpops = lr_postimport_postvacc2
  ) |>
    # assume all E vaccinations are failures so not coding a move to R
    movepopsize_byname(round(se_doses["S_dummygroup1X2"] * theta_infant),
      namebefore = "S_dummygroup1X2",
      nameafter = "R_dummygroup1X2"
    ) |>
    movepopsize_byname(round(se_doses["S_dummygroup1X3"] * theta_adult),
      namebefore = "S_dummygroup1X3",
      nameafter = "R_dummygroup1X3"
    ) |>
    movepopsize_byname(round(se_doses["S_dummygroup1X4"] * theta_adult),
      namebefore = "S_dummygroup1X4",
      nameafter = "R_dummygroup1X4"
    ) |>
    movepopsize_byname(round(se_doses["S_dummygroup1X2"] * (1 - theta_infant)),
      namebefore = "S_dummygroup1X2",
      nameafter = "Sv_dummygroup1X2"
    ) |>
    movepopsize_byname(round(se_doses["S_dummygroup1X3"] * (1 - theta_adult)),
      namebefore = "S_dummygroup1X3",
      nameafter = "Sv_dummygroup1X3"
    ) |>
    movepopsize_byname(round(se_doses["S_dummygroup1X4"] * (1 - theta_adult)),
      namebefore = "S_dummygroup1X4",
      nameafter = "Sv_dummygroup1X4"
    )

  vecpopsize <- setNames(tblpopsize$popsize, tblpopsize$updatedstate)

  # Fifth simulation: third vaccination day to end of simulation (365 days)
  # remaining days, after immunity onset delay taken into account
  postimport_vacc3 <- wrap_adaptivetau(
    vecpopsize,
    measlescompiled,
    rate_func = NULL,
    c(
      beta = config$params$r_0 / config$params$infectious_period,
      tau = config$params$infectious_period, int_eff = config$params$int_eff,
      taue = config$params$latent_period
    ),
    (
      365 - config$history$vaccine_days[3] -
        config$history$immunity_onset_delay
    ),
    1,
    "adaptivetau"
  )[[1]]

  pre_import_mod <- pre_import
  pre_import_mod <- pre_import_mod[-nrow(pre_import_mod), ]

  import_prevacc_mod <- import_prevacc
  import_prevacc_mod$time <- import_prevacc_mod$time +
    config$history$import_days
  import_prevacc_mod <- import_prevacc_mod[-nrow(import_prevacc_mod), ]

  postimport_vacc1_mod <- postimport_vacc1
  postimport_vacc1_mod$time <- postimport_vacc1_mod$time +
    (config$history$import_days +
      (
        config$history$import_days + config$history$vaccine_days[1]
      ) +
      config$history$immunity_onset_delay)
  postimport_vacc1_mod <- postimport_vacc1_mod[-nrow(postimport_vacc1_mod), ]

  postimport_vacc2_mod <- postimport_vacc2
  postimport_vacc2_mod$time <- postimport_vacc2_mod$time +
    (config$history$import_days +
      (
        config$history$import_days + config$history$vaccine_days[1]
      ) +
      config$history$immunity_onset_delay +
      (config$history$vaccine_days[2] - config$history$vaccine_days[1]))
  postimport_vacc2_mod <- postimport_vacc2_mod[-nrow(postimport_vacc2_mod), ]

  postimport_vacc3_mod <- postimport_vacc3
  postimport_vacc3_mod$time <- postimport_vacc3_mod$time +
    (config$history$import_days +
      (
        config$history$import_days + config$history$vaccine_days[1]
      ) +
      config$history$immunity_onset_delay +
      (config$history$vaccine_days[2] - config$history$vaccine_days[1]) +
      (config$history$vaccine_days[3] - config$history$vaccine_days[2])
    )
  # don't remove last row on last sim

  all_dyn <- rbind(
    pre_import, import_prevacc_mod, postimport_vacc1_mod,
    postimport_vacc2_mod, postimport_vacc3_mod
  )

  # Incorporating case ascertainment delay
  # Creating obs_df to hold delayed times when incidence changes
  # NOTE: this will have to be plotted separately, as plot_stoch_model()
  # expects a list of data frames and handles simulation ids.

  # nolint start object_name_linter
  if (sum(all_dyn$incI_cum_latent) > 0) {
    incI_latent <- cuminc2inc(
      all_dyn$time,
      all_dyn$incI_cum_latent
    )
    all_dyn$incI_latent <- incI_latent
    time_mod <-
      delaytime_exp(
        all_dyn[all_dyn$"incI_latent" > 0, ]$time,
        config$params$ascertainment_delay
      )
    obs_df1 <- data.frame(
      time_report = time_mod,
      incI_latent_d = all_dyn[which(all_dyn$incI_latent > 0), ]$incI_cum_latent
    )
    sub <- all_dyn |> dplyr::select(time, incI_cum_latent, incI_latent)
    obs_df <- dplyr::full_join(sub, obs_df1, by = c("time" = "time_report")) |>
      dplyr::arrange(time)
  } else {
    print("No incidence in this simulation, returning empty obs_df")
    sub <- all_dyn |> dplyr::select(time, incI_cum_latent)
    obs_df1 <- data.frame(
      time_report = numeric(0),
      incI_latent_d = numeric(0)
    )
    obs_df <- dplyr::full_join(sub, obs_df1, by = c("time" = "time_report"))
  }
  # nolint end object_name_linter

  list(all_dyn, obs_df)
}


measles_dyn <- list(run_measles_chicago_2024("scripts/measles_params.yaml"))
measles_states <- names(measles_dyn[[1]][[1]][-1])
dynamics <- measles_dyn[[1]][1]

# Looking at latent and delayed cumulative incidence
inctbl <- measles_dyn[[1]][[2]]
ggplot(inctbl, aes(x = time)) +
  geom_line(aes(y = incI_cum_latent, color = "Latent")) +
  geom_line(aes(y = incI_latent_d, color = "Observed")) +
  theme_classic() +
  labs(
    x = "Time (days)",
    y = "Cumulative Incidence",
    title = "Chicago Measles 2024 Simulation: Incidence Dynamics"
  ) +
  scale_color_manual(
    name = "",
    values = c("Latent" = "blue", "Observed" = "red")
  )

# Running multiple simulations
config_path <- "scripts/measles_params.yaml"
sims <- list()
sims <- purrr::map(1:50, ~ run_measles_chicago_2024(config_path))
dynamics <- purrr::map(sims, 1)

# Run simulations in parallel
sims <- list()
future::plan(future::multisession, workers = parallel::detectCores() - 1)
sims <- furrr::future_map(1:100, ~ run_measles_chicago_2024(config_path),
  .options =
    furrr::furrr_options(
      seed = TRUE,
      packages = "compModels"
    )
)
dynamics <- purrr::map(sims, 1)


# Plot compartments
plot_stoch_model(dynamics,
  compartments = measles_states[grep("^(S_)", measles_states)],
  colors =
    colorRampPalette(RColorBrewer::brewer.pal(5, "Dark2"))(5)
)

plot_stoch_model(dynamics,
  compartments = measles_states[grep("^(Sv_)", measles_states)],
  colors =
    colorRampPalette(RColorBrewer::brewer.pal(5, "Dark2"))(5)
)

plot_stoch_model(dynamics,
  compartments = measles_states[grep("^(E_)", measles_states)],
  colors =
    colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(10)
)

plot_stoch_model(dynamics,
  compartments = measles_states[grep("^(I_)", measles_states)],
  colors =
    colorRampPalette(RColorBrewer::brewer.pal(5, "Dark2"))(5)
)

plot_stoch_model(dynamics,
  compartments = measles_states[grep("^(R_)", measles_states)],
  colors =
    colorRampPalette(RColorBrewer::brewer.pal(5, "Dark2"))(5)
)

plot_stoch_model(dynamics,
  compartments = measles_states[grep("incI_cum_latent", measles_states)],
  colors =
    colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(1)
)
