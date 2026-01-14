# nolint start
SIR <- define_states(c("S", "I", "R")) |>
  add_transition("I", "R", "tau") |>
  add_infection("I", "S", "I", "beta", normlogic = FALSE)
SIRc <- compilemodel(SIR)
inputpeter <- SIR

SI2R <- define_states(c("S", "I", "R")) |>
  add_transition("I", "R", "tau", chainlength = 2) |>
  add_infection("I", "S", "I", "beta")
SI2Rc <- compilemodel(SI2R)
inputpeter <- SI2R
SI2Rc

SIRgroups <- define_states(c("S", "I", "R")) |>
  add_infection("I", "S", "I", "beta", processname = "infection") |>
  add_transition("I", "R", "tau", processname = "recovery", chainlength = 3) |>
  add_transition("I", "R", "tau", processname = "recovery2", chainlength = 5) |>
  add_group(
    groupname = c("old", "young"), grouptype = "Age",
    scaleprocessbyname = list(infection = c(1, 2))
  )

SIRgroups2 <- SIRgroups |>
  add_group(
    groupname = c("Patient", "HCW"), grouptype = "Status",
    scaletransitions = c(.5, 1)
  ) |>
  combine_groups(c("Age", "Status"))

SIRgroupsc <- compilemodel(SIRgroups)
inputpeter <- SIRgroups2
SIRgroups2c <- compilemodel(SIRgroups2)

SIRgroups2c$modeloutstructions$updatedstates

SIRgroups2c$modelinstructions$tblprocesses

SIRgroups_trans <-
  SIRgroups |>
  add_grouptransition("young",
    c("old", "old2"),
    "tau_old",
    chainlength = 2,
    forkprobability = c(.1, .9)
  )
SIRgroups_transc <- compilemodel(SIRgroups_trans)

SIRgroups2_trans <- SIRgroups2 |>
  add_grouptransition("young", "old", "tau_old") |>
  add_grouptransition("Patient", "HCW", "tau_Patient")
SIRgroups2_transc <- compilemodel(SIRgroups2_trans)

inputpeter <- SIRgroups_trans

compilemodel(SIRgroups2)

SIRmeta <- define_states(c("S", "I", "R")) |>
  add_infection("I", "S", "I", "beta") |>
  add_transition("I", "R", "tau") |>
  define_metapopulations(c("School", "Work", "Church")) |>
  add_travel("mu")

SIRmetac <- compilemodel(SIRmeta)
inputpeter <- SIRmeta

SIRmetac$modelinstructions$tblprocesses |> dplyr::filter(processlabel == "migration")

SIRgroupsc <- compilemodel(SIRgroups)

SI2RV <- define_states(c("S", "I", "R", "V")) |>
  add_transition("I", "R", "tau", chainlength = 2) |>
  add_infection("I", "S", "I", "beta") |>
  add_transition("S", "V", "tauv", processname = "vax") |>
  add_transition("I", "V", "tauv", processname = "vax")
SI2RVc <- compilemodel(SI2RV)
inputpeter <- SI2RV

basestates <- c("S", "E", "I", "R")
metapopulation <- c("USA", "UK")

peterlist <- define_states(basestates)

test <-
  make_metapopulations(metapopulation,
    scaleprocessbyname = list(
      process1 = c(1, 2),
      process2 = c("3", "4")
    )
  )

wat <- define_states(basestates) |>
  add_metapopulations() |>
  add_travel("mu")
wat <- define_states(basestates) |>
  define_metapopulations(metapopulation,
    scaleprocessbyname = list(travel = 2)
  ) |>
  add_travel("mu", "UK", "USA", processname = "travel") |>
  add_group(c("Old", "Young"),
    grouptype = "Age",
    scalemigrations = c(1, 2),
    scaleprocessbyname = list(recovery = c(3.5, 2.5)),
    scaleprocessbygroup = list(progression = c(3, 4))
  ) |>
  add_group(c("Patient", "HCW"),
    grouptype = "Job",
    scaleprocessbyname = list(travel = c(3, 4))
  ) |>
  combine_groups(c("Job", "Age")) |>
  add_transition("I", c("R", "S"),
    chainlength = 3,
    meantime = "tau",
    groupname = list(
      Age = c("Old", "Young"),
      Job = ("Patient")
    ),
    processname = "recovery",
    processgroup = "progression"
  ) |>
  add_transition("E", "I",
    chainlength = 1, meantime = "taue",
    processgroup = "progression"
  ) |>
  add_infection("I", "S", "E", "beta")

inputpeter <- wat
compilewat <- compilemodel(wat)

readr::write_csv(tblupdatestate, "~/Documents/tblstate_complex.csv")

sir <- define_states(c("S", "I", "R")) |>
  add_infection("I", "S", "I", "beta") |>
  add_transition("I", "R", "tau")
inputpeter <- sir
readr::write_csv(tblstates, "~/Documents/tblstate_sir.csv")
si2r <- define_states(c("S", "I", "R")) |>
  add_infection("I", "S", "I", "beta") |>
  add_transition("I", "R", "tau", chainlength = 2)
inputpeter <- si2r
readr::write_csv(tblstates, "~/Documents/tblstate_si2r.csv")

sirmeta <- define_states(c("S", "I", "R")) |>
  add_infection("I", "S", "I", "beta") |>
  add_transition("I", "R", "tau") |>
  define_metapopulations(c("UK", "USA"))
inputpeter <- sirmeta
readr::write_csv(tblstates, "~/Documents/tblstate_sirmeta.csv")

sirgroup <- define_states(c("S", "I", "R")) |>
  add_infection("I", "S", "I", "beta") |>
  add_transition("I", "R", "tau") |>
  add_group(c("Old", "Young"))
inputpeter <- sirgroup
readr::write_csv(tblstates, "~/Documents/tblstate_sirgroup.csv")
tibble::tibble(
  metapopulation = metapopulation,
  interactionscale = 1,
  transitionscale = 1,
  basestates = list(base_state)
)


base_states <- c("S", "E", "I", "R", "Sv")
measlesmodel <- define_states(base_states) |>
  add_infection("I", "S", "E", "beta") |>
  add_transition("I", "R", "tau") |>
  add_transition("E", "I", "taue", chainlength = 2) |>
  add_transition("S", c("Sv", "R"), "1/(theta)",
    forkprobability = c("p", "1-p"),
    groupname = 3:4, processname = "vax"
  ) |>
  add_transition("S", c("Sv", "R"), "1/(theta)",
    forkprobability = c("p_infant", "1-p_infant"),
    groupname = 2, processname = "vax"
  ) |>
  add_infection("I", "Sv", "E", "beta") |>
  add_group(1:5, grouptype = "Age")

measlescompiled <- compilemodel(measlesmodel)
tblupdatedstates <- measlescompiled$modelinstructions$tblupdatedstates

temp <- tblupdatedstates |> dplyr::mutate(x0 = 0)
Slogic <- tblupdatedstates$basestates == "S"
temp$x0[Slogic] <- 100
Ilogic <- tblupdatedstates$basestates == "I"
temp$x0[Ilogic] <- 2

x0 <- temp |> dplyr::pull(x0)

measles_pars <-
  c(beta = 5, tau = 5, taue = 8, theta = .01, p = 1 - 0.925, p_infant = 1 - .84)


tblupdatedstates_firstchain <- tblupdatedstates |> dplyr::filter()


tf <- 100
run_adaptivetau(x0, measlescompiled,
  rate_func = NULL,
  measles_pars, tf,
  method = "exact"
)


dyn <- wrap_adaptivetau(
  c("S" = 999, "I" = 1, "R" = 0, "E" = 0, "Sv" = 0),
  measlescompiled,
  rate_func = NULL, # defaults to compModels::generalized_rates()
  c(
    beta = 5, tau = 5, taue = 8,
    theta = .01, p = 1 - 0.925, p_infant = 1 - .84
  ),
  25,
  1,
  "adaptivetau"
)


test <- modeloutput2odeinput(measlescompiled)
dyn[nrow(dyn), ]

tstitch <- c(1, 10, 20)



dhqp <- define_states(c("S", "E", "A", "I", "R", "Tr", "H")) |>
  add_group(c("Shortstay", "Longstay", "HCP"), grouptype = "Resident") |>
  add_group(c("non", "pro"), grouptype = "Vaxstatus") |>
  combine_groups(c("Resident", "Vaxstatus"))
add_replacement(c("S", "E", "A", "I", "R", "Tr", "H"), "S",
  rate = "mu_SR", groupname = list(), processname = "migrate_in"
)


# SCI-TEST XXX

SIR <- define_states(c("S", "I", "R")) |>
  add_transition("I", "R", "tau") |>
  add_infection("I", "S", "I", "beta")
SIRc <- compilemodel(SIR)

Nsims <- 1000
N <- 1000
pars <- c(beta = 2, tau = 1, N = N)
X0 <- c(S = N - 1, I = 1, R = 0)
maxt <- 25
sciout <- wrap_adaptivetau(
  X0,
  SIRc,
  rate_func = NULL,
  pars,
  maxt,
  Nsims,
  "adaptivetau"
)

sciout_t <- sapply(sciout, function(x) {
  x |>
    dplyr::filter(I == 0) |>
    dplyr::distinct(I, .keep_all = TRUE) |>
    dplyr::pull(time)
})

prob <- c(alpha = 0)
tdependent_extinction <- function(time, prob, pars) {
  with(as.list(c(prob, pars)), {
    dalpha <- (beta * alpha - (1 / tau)) * (alpha - 1)
    list(c(dalpha))
  })
}
tode <- seq(0, maxt, .01)
library(deSolve)
out <-
  deSolve::ode(
    y = prob,
    times = tode,
    func = tdependent_extinction,
    parms = pars
  )




tibbleplot <- tibble::tibble(
  tdeath = sort(unlist(sciout_t)),
  alpha = (1:length(unlist(sciout_t))) / Nsims
)

ggplot(out, aes(x = time, y = alpha), color = blue) +
  geom_line() +
  geom_line(data = tibbleplot, aes(x = tdeath, y = alpha), color = "red")


sciout_t_filter <- sapply(sciout, function(x) {
  x |>
    dplyr::filter(I == 0) |>
    dplyr::distinct(I, .keep_all = TRUE) |>
    dplyr::pull(time)
})

bigthresh <- 10
sciout_threshlogic <- sapply(sciout, function(x) {
  max(x |> dplyr::pull(I)) < bigthresh
})

sciout_filterplot <- unlist(sciout_t_filter[sciout_threshlogic])

tibbleplot_filter <-
  tibble::tibble(
    tdeath = sort(sciout_filterplot),
    alpha = (1:length(sciout_filterplot)) / Nsims
  )

ggplot(out, aes(x = time, y = alpha), color = blue) +
  geom_line() +
  geom_line(
    data = tibbleplot_filter,
    aes(x = tdeath, y = alpha), color = "red"
  ) +
  xlim(0, max(sciout_filterplot))


# exmple popsize script
basestates <- c("S", "I", "R")
modelinstructions <- define_states(basestates) |>
  add_transition("I", "R", "tau", chainlength = 2) |>
  add_infection("I", "S", "I", "beta") |>
  add_group(c("group1", "group2"), grouptype = "type1")
compiledmodel <- compilemodel(modelinstructions)

tblpopsize <- define_popsize(compiledmodel) |>
  setpopsize_byfeature(.99, basestates = c("S"), groupnames = list(type1 = c("group1", "group2")), chains = 1)

# debug
basestates <- c("S", "I", "R")

modelinstructions <- define_states(basestates) |>
  add_transition("I", "R", "tau", chainlength = 2) |>
  add_infection("I", "S", "I", "beta") |>
  add_group(c("group1", "group2"), grouptype = "type1") |>
  add_group(c("group3", "group4"), grouptype = "type2") |>
  combine_groups()
compiledmodel <- compilemodel(modelinstructions)

modelinstructions <- define_states(basestates) |>
  add_transition("I", "R", "tau", chainlength = 1) |>
  add_infection("I", "S", "I", "beta", groupnames = list(type1 = "group1")) |>
  add_group(c("group1", "group2"), grouptype = "type1") |>
  add_group(c("group3", "group4"), grouptype = "type2") |>
  combine_groups() |>
  add_infection("I", "S", "I", "beta", groupnames = "group3", crossgroupnames = "group4")
compiledmodel <- compilemodel(modelinstructions)

modelinstructions <- define_states(basestates) |>
  add_transition("I", "R", "tau", chainlength = 1) |>
  add_infection("I", "S", "I", "beta", groupnames = list("group1", "group2")) |>
  add_group(c("group1", "group2"), grouptype = "type1")
compiledmodel <- compilemodel(modelinstructions)

inputpeter <- modelinstructions
compiledmodel$modelinstructions$tblprocesses$rate
compiledmodel$modelinstructions$tblprocesses$percapitarate




basestates <- c("S", "I", "R")

modelinstructions <- define_states(basestates) |>
  add_transition("I", "R", "tau") |>
  add_infection("I", "S", "I", "beta") |>
  add_group(c("group1", "group2"), grouptype = "type1")
compiledmodel <- compilemodel(modelinstructions)

modelinstructions <- define_states(basestates) |>
  add_transition("I", "R", "tau") |>
  add_infection("I", "S", "I", "beta", groupnames = "group1") |>
  add_group(c("group1", "group2"), grouptype = "type1")
compiledmodel <- compilemodel(modelinstructions)

modelinstructions <- define_states(basestates) |>
  add_transition("I", "R", "tau") |>
  add_infection("I", "S", "I", "beta", groupnames = "group1", crossgroupnames = "group2") |>
  add_group(c("group1", "group2"), grouptype = "type1")
compiledmodel <- compilemodel(modelinstructions)

modelinstructions <- define_states(basestates) |>
  add_transition("I", "R", "tau") |>
  add_infection("I", "S", "I", "beta", groupnames = "group1", crossgroupnames = "group2", symmetric = FALSE) |>
  add_group(c("group1", "group2"), grouptype = "type1")
compiledmodel <- compilemodel(modelinstructions)

compiledmodel$modelinstructions$tblprocesses
compiledmodel$modeloutstructions
inputpeter <- modelinstructions
# nolint end
