# ==============================================================================
# 03_run_simulations.R
#
# Willem M. (Wim) Otte, w.m.otte@umcutrecht.nl
#
# Monte-Carlo probabilistic sensitivity analysis (PSA) over the full scenario
# grid: country (Nigeria/Ecuador/UK) x severity (severe/less-severe) x
# treatment-gap level. Each scenario is run for CFG$n_psa iterations. Published
# point estimates and plausible limits define triangular distributions, avoiding
# artificial point masses at bounds from clamped normal draws.
#
# Output (data/derived/):
#   sim_traces.rds   : list with $traces (cohort occupancy, all iters) and
#                      $params (sampled parameters per iter) and $grid
# ==============================================================================

if (!exists("CFG")) {
  cand <- c("R/00_config.R", "00_config.R",
            file.path(Sys.getenv("CAE_ROOT", "."), "R/00_config.R"))
  source(cand[file.exists(cand)][1])
}
source(file.path(CFG$root, "R", "02_markov_model.R"))
source(file.path(CFG$root, "R", "lib_cea.R"))
suppressPackageStartupMessages(library(tidyverse))

set.seed(CFG$seed)

mort <- read_tsv(file.path(CFG$dir_derived, "mortality_schedule.tsv"),
                 show_col_types = FALSE)

# ---- Scenario grid -----------------------------------------------------------
grid <- CFG$countries %>%
  select(country, iso3, setting) %>%
  crossing(severity = c("severe", "less_severe")) %>%
  rowwise() %>%
  mutate(gap = list(CFG$gap_sweep[[setting]])) %>%
  unnest(gap) %>%
  ungroup()

message(sprintf("[03_run] scenario grid: %d rows x %d PSA iters = %d cohort runs",
                nrow(grid), CFG$n_psa, nrow(grid) * CFG$n_psa))

# ---- Pre-sample parameters per (setting, iter) so paired scenarios share draws
# (common random numbers -> stable incremental differences for CEA). The sampler
# lives in lib_cea.R so the co-author sensitivity script draws from the same one.
samp <- sample_psa_params(CFG$n_psa)

# ---- Run all scenarios x iterations ------------------------------------------
run_scenario_row <- function(country, iso3, setting, severity, gap) {
  ctrl_col <- if (setting == "LMIC") "control_cum_LMIC" else "control_cum_HIC"
  ctry <- country
  bg_rates <- mort %>% filter(.data$country == .env$ctry) %>% arrange(cycle) %>% pull(bg_rate)
  pmap_dfr(samp, function(...) {
    s <- list(...)
    p <- list(
      bg_rates     = bg_rates,
      bg_rate      = bg_rates[1],
      smr          = list(Untreated = s$smr_U, Treated = s$smr_T,
                          Controlled = s$smr_SF, Remission = s$smr_SF),
      initial_gap  = gap,
      treat_uptake = CFG$background_treatment_uptake,
      spont_rem    = s$spont_rem,
      control_prob = 1 - (1 - s[[ctrl_col]])^(1 / CFG$n_cycles)
    )
    run_cohort(p) %>% mutate(iter = s$iter)
  })
}

traces <- grid %>%
  mutate(.row = row_number()) %>%
  group_split(.row) %>%
  map_dfr(function(g) {
    run_scenario_row(g$country, g$iso3, g$setting, g$severity, g$gap) %>%
      mutate(country = g$country, setting = g$setting,
             severity = g$severity, gap = g$gap)
  })

out <- list(traces = traces, params = samp, grid = grid)
saveRDS(out, file.path(CFG$dir_derived, "sim_traces.rds"))

# Spot-check initial untreated occupancy; it must equal the scenario gap exactly.
chk <- traces %>%
  filter(cycle == 0, state == "Untreated") %>%
  group_by(country, setting, gap) %>%
  summarise(initial_untreated = mean(prob), .groups = "drop")
stopifnot(all(abs(chk$initial_untreated - chk$gap) < 1e-12))
message("[03_run] initial untreated share equals scenario gap:")
print(chk, n = 30)
message(sprintf("[03_run] saved %s (%d trace rows)",
                "sim_traces.rds", nrow(traces)))
