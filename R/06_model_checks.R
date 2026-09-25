# ==============================================================================
# 06_model_checks.R
#
# Willem M. (Wim) Otte, w.m.otte@umcutrecht.nl
#
# Executable internal-verification checks for the revised cohort model.
# Writes a human-readable record for the supplement and fails on any violation.
# ==============================================================================

if (!exists("CFG")) {
  cand <- c("R/00_config.R", "00_config.R",
            file.path(Sys.getenv("CAE_ROOT", "."), "R/00_config.R"))
  source(cand[file.exists(cand)][1])
}
source(file.path(CFG$root, "R", "02_markov_model.R"))
source(file.path(CFG$root, "R", "lib_cea.R"))
suppressPackageStartupMessages(library(tidyverse))

mort <- read_tsv(file.path(CFG$dir_derived, "mortality_schedule.tsv"),
                 show_col_types = FALSE)
sim <- readRDS(file.path(CFG$dir_derived, "sim_traces.rds"))

checks <- list()
record <- function(name, pass, detail) {
  checks[[length(checks) + 1L]] <<- tibble(check = name, pass = pass, detail = detail)
}

# Transition matrices remain stochastic at every modelled age and country.
matrix_ok <- TRUE
for (ctry in unique(mort$country)) {
  rates <- mort %>% filter(country == ctry) %>% arrange(cycle) %>% pull(bg_rate)
  for (rate in rates) {
    M <- build_transition_matrix(list(
      bg_rate = rate, smr = CFG$smr,
      treat_uptake = CFG$background_treatment_uptake,
      spont_rem = CFG$spont_remission_prob,
      control_prob = CFG$control_prob[[if (ctry == "UK") "HIC" else "LMIC"]],
      relapse = CFG$relapse_base))
    matrix_ok <- matrix_ok && all(M >= 0) && max(abs(rowSums(M) - 1)) < 1e-12
  }
}
record("Transition matrices are stochastic", matrix_ok,
       "All entries non-negative and every row sums to 1 within 1e-12")

# Initial-state labels equal modelled occupancy exactly.
entry <- sim$traces %>% filter(cycle == 0) %>%
  select(country, severity, gap, iter, state, prob) %>%
  filter(state %in% c("Untreated", "Treated")) %>%
  mutate(target = if_else(state == "Untreated", gap, 1 - gap))
entry_ok <- max(abs(entry$prob - entry$target)) < 1e-12
record("Initial treatment gap matches occupancy", entry_ok,
       "Untreated=gap and Treated=1-gap in every country, severity, and draw")

# Death is absorbing, so cumulative occupancy cannot fall.
death_ok <- sim$traces %>% filter(state == "Death") %>%
  arrange(country, severity, gap, iter, cycle) %>%
  group_by(country, severity, gap, iter) %>%
  summarise(ok = all(diff(prob) >= -1e-12), .groups = "drop") %>%
  summarise(all_ok = all(ok)) %>% pull(all_ok)
record("Death occupancy is monotonic", death_ok,
       "No cohort has a decline in cumulative Death occupancy")

# Triangular sampling is continuous and never produces clamping masses.
ps <- sim$params
# Bounds are read from CFG, never restated here. A literal copy would silently
# agree with a sampler that had drifted away from the declared configuration.
bounds <- tribble(
  ~parameter, ~lo, ~hi,
  "smr_U",      CFG$smr_sens$Untreated[1], CFG$smr_sens$Untreated[2],
  "smr_T",      CFG$smr_sens$Treated[1],   CFG$smr_sens$Treated[2],
  "smr_SF",     CFG$smr_sens$SeizureFree[1], CFG$smr_sens$SeizureFree[2],
  "spont_rem",  CFG$spont_remission_sens[1], CFG$spont_remission_sens[2],
  "dw_sev_UT",  dw_bounds("Untreated_severe")[1], dw_bounds("Untreated_severe")[2],
  "dw_less_UT", dw_bounds("Untreated_less_severe")[1],
                dw_bounds("Untreated_less_severe")[2],
  "dw_control", dw_bounds("Controlled")[1], dw_bounds("Controlled")[2],
  "cost_mult",  CFG$cost_mult_sens[1], CFG$cost_mult_sens[2],
  "case_find_cost", CFG$cost_case_finding_sens[1], CFG$cost_case_finding_sens[2],
  "untreated_cost_frac", CFG$cost_untreated_frac_psa[1],
                         CFG$cost_untreated_frac_psa[2])
boundary_ok <- pmap_lgl(bounds, function(parameter, lo, hi) {
  x <- ps[[parameter]]
  all(x > lo, x < hi) && !any(x == lo | x == hi)
})
record("PSA has no clamped boundary masses", all(boundary_ok),
       paste(bounds$parameter, collapse = ", "))

# Shared quantities have a single sampled column by construction.
shared_ok <- all(c("smr_T", "smr_SF", "dw_sev_UT", "dw_less_UT") %in% names(ps)) &&
  !any(c("smr_TC", "smr_C", "smr_R", "dw_T", "dw_C", "dw_R") %in% names(ps)) &&
  CFG$smr$Controlled == CFG$smr$Remission
record("Shared clinical parameters use shared draws", shared_ok,
       paste("One SMR draw for Controlled and Remission; one DW draw for",
             "less-severe Untreated and Treated; one DW draw for Controlled",
             "and Remission"))

# The mortality schedule follows the index cohort and rises over the horizon.
age_ok <- mort %>% group_by(country) %>%
  summarise(ok = all(diff(age_mid) == CFG$cycle_years) &&
                 first(age_mid) == CFG$index_age + CFG$cycle_years / 2 &&
                 all(diff(bg_rate) >= 0), .groups = "drop") %>%
  summarise(all_ok = all(ok)) %>% pull(all_ok)
record("Age-varying mortality schedule is ordered", age_ok,
       sprintf("Index age %d; rates evaluated at cycle midpoints", CFG$index_age))

# Deterministic toy cases with closed-form answers, to calibrate the machinery
# independently of the fitted parameters.
toy <- function(bg_rate, smr_u, spont, control, gap) {
  run_cohort(list(
    bg_rate = bg_rate, bg_rates = rep(bg_rate, CFG$n_cycles),
    smr = list(Untreated = smr_u, Treated = 1, Controlled = 1, Remission = 1),
    initial_gap = gap, treat_uptake = 0,
    spont_rem = spont, control_prob = control, relapse = 0))
}
occ_at <- function(tr, k, st) {
  tr$prob[tr$cycle == k & tr$state == st]
}

# 1. No mortality and no transitions: the cohort never leaves its initial vector.
frozen <- toy(bg_rate = 0, smr_u = 1, spont = 0, control = 0, gap = 0.7)
frozen_ok <- abs(occ_at(frozen, CFG$n_cycles, "Untreated") - 0.7) < 1e-12 &&
  abs(occ_at(frozen, CFG$n_cycles, "Treated") - 0.3) < 1e-12 &&
  abs(occ_at(frozen, CFG$n_cycles, "Death")) < 1e-12
record("Toy case: inert cohort is stationary", frozen_ok,
       "Zero mortality, remission, and control leaves x_k = x_0 for every cycle")

# 2. Everyone untreated, only mortality acts. Survival is exp(-r*T) exactly.
rate <- 0.01
decay <- toy(bg_rate = rate, smr_u = 1, spont = 0, control = 0, gap = 1)
expected_death <- 1 - exp(-rate * CFG$horizon_y)
decay_ok <- abs(occ_at(decay, CFG$n_cycles, "Death") - expected_death) < 1e-12
record("Toy case: pure exponential survival", decay_ok,
       sprintf("Death at 10 y = 1-exp(-%.2f x %d) = %.6f",
               rate, CFG$horizon_y, expected_death))

# Reward calculations are checked against hand-solvable cases. These use the
# same functions as 04_dalys_costs_cea.R and 10_structural_scenarios.R, so the
# checks cover the accounting code rather than only the state trace.
toy_t <- seq(0, 10, by = 2)
linear_area <- integrate_curve(toy_t, toy_t)
quadratic_hybrid_area <- integrate_curve(toy_t, toy_t^2)
integration_ok <- abs(linear_area - 50) < 1e-12 &&
                  abs(quadratic_hybrid_area - (512 / 3 + 164)) < 1e-12
record("Reward case: within-cycle integration has known area", integration_ok,
       "Linear area is 50; Simpson 0-8 plus trapezoid 8-10 gives 334.666667 for t^2")

frozen_trace <- data.frame(
  Untreated = rep(0.7, 6), Treated = rep(0.3, 6),
  Controlled = rep(0, 6), Remission = rep(0, 6), Death = rep(0, 6))
toy_dw <- c(Untreated = 0.4, Treated = 0.2, Controlled = 0,
            Remission = 0, Death = 0)
toy_yld <- discounted_yld(toy_t, frozen_trace, toy_dw, discount = 0)
yld_ok <- abs(toy_yld - 3.4) < 1e-12
record("Reward case: frozen-cohort YLD has closed form", yld_ok,
       "10 years x (0.7 x 0.4 + 0.3 x 0.2) = 3.4 YLD")

mort_trace <- frozen_trace
mort_trace$Death <- seq(0, 0.3, length.out = 6)
mort_trace$Untreated <- 1 - mort_trace$Death
mort_trace$Treated <- 0
toy_yll <- discounted_yll(toy_t, mort_trace, rep(20, 5), discount = 0)
yll_ok <- abs(toy_yll - 6) < 1e-12
record("Reward case: cycle deaths give known YLL", yll_ok,
       "Cumulative deaths 0.3 x 20 residual years = 6 YLL at zero discount")

toy_reward <- value_trace(
  toy_t, frozen_trace, toy_dw, rep(20, 5), annual_cost = 100,
  untreated_cost_fraction = 0.1, comparator_gap = 0.8,
  intervention_gap = 0.2, initiation_unit_cost = 50, discount = 0)
cost_ok <- abs(toy_reward[["care_cost"]] - 370) < 1e-12 &&
           abs(toy_reward[["detect_cost"]] - 30) < 1e-12 &&
           abs(toy_reward[["cost"]] - 400) < 1e-12
record("Reward case: care and initiation costs have closed form", cost_ok,
       "Care is 10 x 100 x (0.3 + 0.1 x 0.7) = 370; initiation is 30")

toy_cea <- incremental_cea(
  comparator_daly = 10, comparator_cost = 400,
  intervention_daly = 7, intervention_cost = 700, threshold = 200)
cea_ok <- abs(toy_cea$daly_averted - 3) < 1e-12 &&
          abs(toy_cea$incr_cost - 300) < 1e-12 &&
          abs(toy_cea$icer - 100) < 1e-12 && abs(toy_cea$inb - 300) < 1e-12
record("Reward case: incremental CEA and INB close exactly", cea_ok,
       "Delta DALY=3, delta cost=300, ICER=100, and INB=200 x 3 - 300 = 300")

# The disability-weight structure used by the accounting code must equal the one
# declared in CFG$dw. A revert once landed in 04 only, leaving the deterministic
# sensitivity analysis running a different model than the headline while every
# gate stayed green. dw_vec() is now the single definition; this check pins it to
# the declared configuration so a silent edit to either side fails the build.
point <- list(dw_sev_UT  = CFG$dw$Untreated_severe,
              dw_less_UT = CFG$dw$Untreated_less_severe,
              dw_control = CFG$dw$Controlled)
expect_sev  <- c(Untreated = CFG$dw$Untreated_severe, Treated = CFG$dw$Treated,
                 Controlled = CFG$dw$Controlled, Remission = CFG$dw$Remission,
                 Death = 0)
expect_less <- c(Untreated = CFG$dw$Untreated_less_severe, Treated = CFG$dw$Treated,
                 Controlled = CFG$dw$Controlled, Remission = CFG$dw$Remission,
                 Death = 0)
dw_ok <- isTRUE(all.equal(dw_vec(point, "severe"), expect_sev)) &&
         isTRUE(all.equal(dw_vec(point, "less_severe"), expect_less))
record("Disability weights match the declared configuration", dw_ok,
       "dw_vec() at point estimates equals CFG$dw for both severity cohorts")

# The Treated/Controlled SMR is one construct sharing one PSA draw. The scenario
# in which Treated keeps its entry disability weight must therefore leave the
# other states untouched, and must actually change the severe Treated weight.
retain <- dw_vec(point, "severe", treated_keeps_entry_dw = TRUE)
retain_ok <- retain[["Treated"]] == CFG$dw$Untreated_severe &&
             retain[["Untreated"]] == CFG$dw$Untreated_severe &&
             retain[["Controlled"]] == CFG$dw$Controlled &&
             dw_vec(point, "less_severe", treated_keeps_entry_dw = TRUE)[["Treated"]] ==
               CFG$dw$Untreated_less_severe
record("DW-RETAIN structural scenario is well formed", retain_ok,
       "Treated inherits the entry weight; Controlled and Remission are unchanged")

# The conventional YLL annuity of e_x is an upper bound on expected discounted
# life-years, because the annuity is concave (Jensen). Verify the direction holds
# at every modelled death age, so the convention is bounded and documented.
lt_ng <- read_tsv(file.path(CFG$dir_derived, "life_table__Nigeria.tsv"),
                  show_col_types = FALSE)
ages  <- CFG$index_age + (1:CFG$n_cycles) * CFG$cycle_years - CFG$cycle_years / 2
le_ng <- mort %>% filter(country == "Nigeria") %>% arrange(cycle) %>% pull(residual_le)
jensen_ok <- all(mapply(function(a, L) {
  discounted_annuity(L, CFG$discount) >= discounted_le_lifetable(lt_ng, a, CFG$discount) - 1e-9
}, ages, le_ng))
record("YLL annuity bounds the life-table survival integral", jensen_ok,
       "a(e_x) >= integral exp(-ru)S(u)du at every cycle midpoint age (Jensen)")

out <- bind_rows(checks)
write_tsv(out, file.path(CFG$dir_tab, "model_checks.tsv"))
if (!all(out$pass)) {
  print(out)
  stop("One or more model checks failed")
}
message(sprintf("[06_checks] %d/%d internal checks passed", sum(out$pass), nrow(out)))
