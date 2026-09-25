# ==============================================================================
# 04_dalys_costs_cea.R
#
# Willem M. (Wim) Otte, w.m.otte@umcutrecht.nl
#
# DALYs (YLD + YLL), direct costs, DALYs averted by bridging the treatment gap,
# and incremental cost-effectiveness ratios (ICERs) vs WHO 1x/3x GDP thresholds.
#
# Outputs (output/tables/):
#   per_person_dalys.tsv  : per-person discounted DALY/cost per scenario (mean+CI)
#   cea_results.tsv       : DALYs averted, incr cost, ICER vs gap_base, threshold
#   mortality_table.tsv   : realized per-cycle & annual death prob per state/ctry
#   population_burden.tsv  : population-level DALYs averted
# ==============================================================================

if (!exists("CFG")) {
  cand <- c("R/00_config.R", "00_config.R",
            file.path(Sys.getenv("CAE_ROOT", "."), "R/00_config.R"))
  source(cand[file.exists(cand)][1])
}
suppressPackageStartupMessages(library(tidyverse))
# Valuation primitives (integrator, disability weights, YLL annuity) live in one
# place so 04, 07 and 10 cannot drift apart. See lib_cea.R.
source(file.path(CFG$root, "R", "lib_cea.R"))

sim  <- readRDS(file.path(CFG$dir_derived, "sim_traces.rds"))
prev <- read_tsv(file.path(CFG$dir_derived, "prevalence.tsv"), show_col_types = FALSE)
bg   <- read_tsv(file.path(CFG$dir_derived, "background_mortality.tsv"),
                 show_col_types = FALSE)
mort <- read_tsv(file.path(CFG$dir_derived, "mortality_schedule.tsv"),
                 show_col_types = FALSE)
r <- CFG$discount

# integrate_curve() and dw_vec() are provided by lib_cea.R.

# ---- Per (scenario, iter) DALY + cost ---------------------------------------
params <- sim$params
traces <- sim$traces
tt <- (0:CFG$n_cycles) * CFG$cycle_years      # time points in years (0..10)

# wide occupancy: cycle x state per scenario/iter
calc_one <- function(df, severity, country, setting, gap, it) {
  p  <- params[params$iter == it, ]
  # Shared PSA columns are expanded here to the state-specific names expected
  # by the accounting code.
  pp <- as.list(p)
  dw <- dw_vec(pp, severity)
  cost_yr <- CFG$cost_annual[[country]] * pp$cost_mult

  occ <- df %>% select(cycle, state, prob) %>%
    pivot_wider(names_from = state, values_from = prob) %>% arrange(cycle)
  # Deaths entering each cycle receive the discounted residual-LE annuity.
  ctry <- country
  resid_le <- mort %>% filter(.data$country == .env$ctry) %>% arrange(cycle) %>%
    pull(residual_le)
  annuity  <- discounted_annuity(resid_le, r)
  reward <- value_trace(
    tt, occ, dw, annuity, cost_yr, pp$untreated_cost_frac,
    CFG$gap_base[[setting]], gap, pp$case_find_cost, r)
  as_tibble_row(reward)
}

per_iter <- traces %>%
  group_by(country, setting, severity, gap, iter) %>%
  group_modify(~ calc_one(.x, .y$severity, .y$country, .y$setting, .y$gap, .y$iter)) %>%
  ungroup()

# ---- Summarise per scenario (mean + 95% CI) ----------------------------------
ci <- function(x) c(mean = mean(x),
                    lo = quantile(x, .025, names = FALSE),
                    hi = quantile(x, .975, names = FALSE))

per_person <- per_iter %>%
  group_by(country, setting, severity, gap) %>%
  summarise(
    daly_mean = mean(daly), daly_lo = quantile(daly, .025),
    daly_hi = quantile(daly, .975),
    cost_mean = mean(cost), cost_lo = quantile(cost, .025),
    cost_hi = quantile(cost, .975),
    yld_mean = mean(yld), yll_mean = mean(yll),
    .groups = "drop"
  )
write_tsv(per_person, file.path(CFG$dir_tab, "per_person_dalys.tsv"))

# ---- CEA: incremental vs the status-quo (worst) gap, common random numbers ---
gap_base_of <- function(setting) CFG$gap_base[[setting]]

cea <- per_iter %>%
  group_by(country, setting, severity, iter) %>%
  group_modify(function(d, key) {
    base <- d %>% filter(gap == gap_base_of(key$setting))
    inc <- incremental_cea(base$daly, base$cost, d$daly, d$cost,
                           CFG$cet_oc[[key$country]])
    bind_cols(d, as_tibble(inc) %>% select(daly_averted, incr_cost))
  }) %>%
  ungroup()   # base gap retained: daly_averted=0, icer=NaN -> "n/a (base case)"

# Keep paired draw-level incremental outcomes so combined-cohort uncertainty is
# calculated from combined draws, not by averaging marginal quantiles.
write_tsv(cea %>% select(country, setting, severity, gap, iter,
                         daly_averted, incr_cost),
          file.path(CFG$dir_tab, "cea_draws.tsv"))

cea_summary <- cea %>%
  group_by(country, setting, severity, gap) %>%
  summarise(
    # NB: compute quantiles from the per-iter vectors BEFORE collapsing the
    # columns to their means, otherwise the CI bounds degenerate to the mean.
    da_lo = quantile(daly_averted, .025), da_hi = quantile(daly_averted, .975),
    ic_lo = quantile(incr_cost, .025),    ic_hi = quantile(incr_cost, .975),
    icer  = mean(incr_cost) / mean(daly_averted),
    # V2: incremental net monetary benefit (INB) with CI - the recommended
    # summary when the averted-DALY denominator can cross zero (an ICER point
    # estimate is then uninterpretable). INB = lambda * dDALY - dCost; >0 means
    # cost-effective at willingness-to-pay lambda. Reported at the
    # opportunity-cost threshold (0.5x GDP) and at 1x GDP.
    inb_oc      = mean(CFG$cet_oc[country[1]] * daly_averted - incr_cost),
    inb_oc_lo   = quantile(CFG$cet_oc[country[1]] * daly_averted - incr_cost, .025),
    inb_oc_hi   = quantile(CFG$cet_oc[country[1]] * daly_averted - incr_cost, .975),
    inb_1x      = mean(CFG$gdp_pc[country[1]] * daly_averted - incr_cost),
    inb_1x_lo   = quantile(CFG$gdp_pc[country[1]] * daly_averted - incr_cost, .025),
    inb_1x_hi   = quantile(CFG$gdp_pc[country[1]] * daly_averted - incr_cost, .975),
    daly_averted = mean(daly_averted),
    incr_cost    = mean(incr_cost),
    .groups = "drop"
  ) %>%
  mutate(
    gdp_pc = CFG$gdp_pc[country],
    thr_1x = gdp_pc, thr_3x = 3 * gdp_pc, thr_oc = CFG$cet_oc[country],
    verdict = case_when(
      !is.finite(icer)        ~ "n/a (base case)",
      icer < 0                ~ "dominant (cost-saving)",
      icer < gdp_pc           ~ "highly cost-effective (<1x GDP)",
      icer < 3 * gdp_pc       ~ "cost-effective (<3x GDP)",
      TRUE                    ~ "not cost-effective (>3x GDP)"
    ),
    # V2: stricter opportunity-cost verdict (0.5x GDP proxy).
    verdict_oc = case_when(
      !is.finite(icer)  ~ "n/a (base case)",
      icer < 0          ~ "dominant (cost-saving)",
      icer < thr_oc     ~ "cost-effective (<opportunity-cost threshold)",
      TRUE              ~ "above opportunity-cost threshold"
    )
  )
write_tsv(cea_summary, file.path(CFG$dir_tab, "cea_results.tsv"))

# ---- Per-draw bridged CEA (for CE-plane quadrant table, Singh follow-up) ------
# Persist the per-iteration incremental cost and averted DALYs at the fully
# bridged gap, so the 1000-draw cloud can be classified into cost-effectiveness
# plane quadrants (dominant / cost-saving = incr_cost<0 & daly_averted>0).
psa_bridged <- cea %>%
  filter(gap == map_dbl(setting, ~ CFG$gap_intervention[[.x]])) %>%
  select(country, setting, severity, iter, daly_averted, incr_cost)
write_tsv(psa_bridged, file.path(CFG$dir_tab, "psa_draws_bridged.tsv"))

# ---- Population-level burden averted -----------------------------------------
# V2 FIX: scale each severity cohort by its SEVERITY-SPECIFIC share of PWE, not
# by the total PWE. Previously the severe per-person figure was multiplied by the
# entire PWE population, overstating the severe national burden ~2.8x.
sev_w <- tibble(severity = names(CFG$severity_frac),
                w = as.numeric(CFG$severity_frac))
pop_burden <- cea_summary %>%
  left_join(prev %>% select(country, prevalence_per_1000, population),
            by = "country") %>%
  left_join(sev_w, by = "severity") %>%
  mutate(
    pwe_total = population * prevalence_per_1000 / 1000,  # all PWE
    pwe       = pwe_total * w,                            # severity-specific PWE
    pop_daly_averted = daly_averted * pwe,
    pop_incr_cost    = incr_cost * pwe
  ) %>%
  select(country, setting, severity, gap, pwe_total, pwe,
         pop_daly_averted, pop_incr_cost, icer, verdict, verdict_oc)
write_tsv(pop_burden, file.path(CFG$dir_tab, "population_burden.tsv"))

# ---- Combined national burden (severe + less-severe cohorts summed) ----------
# Each cohort is already scaled to its own slice of PWE, so the combined national
# total is the simple sum across severities (no further weighting).
pop_combined <- pop_burden %>%
  group_by(country, gap) %>%
  summarise(
    pwe = first(pwe_total),
    pop_daly_averted = sum(pop_daly_averted),
    pop_incr_cost    = sum(pop_incr_cost),
    .groups = "drop"
  ) %>%
  mutate(daly_averted_pp = pop_daly_averted / pwe,
         icer = pop_incr_cost / pop_daly_averted)
write_tsv(pop_combined, file.path(CFG$dir_tab, "population_burden_combined.tsv"))

# ---- Probabilistic cost-effectiveness (handles CIs that cross zero) ----------
# A per-person averted-DALY CI can cross zero in tail draws (a high spontaneous-
# remission + modest untreated-SMR draw routes untreated people into Remission,
# SMR 0.95, rather than Treated, SMR 2.54). Rather than lean on the point estimate
# we report, per scenario, the probability that bridging (i) averts any DALYs and
# (ii) is cost-effective at 1x / 3x GDP, via net monetary benefit
# NMB = lambda * daly_averted - incr_cost (cost-effective when NMB > 0).
cea_gdp <- cea %>% mutate(gdp = CFG$gdp_pc[country], oc = CFG$cet_oc[country])

# The programme reaches the whole prevalent population. Combine severity strata
# within each paired draw so its Monte Carlo error is available from the same
# probability table as the subgroup results.
cea_mixed <- cea_gdp %>%
  mutate(weight = CFG$severity_frac[severity]) %>%
  group_by(country, setting, gap, iter) %>%
  summarise(
    daly_averted = sum(weight * daly_averted),
    incr_cost = sum(weight * incr_cost),
    gdp = first(gdp), oc = first(oc),
    .groups = "drop"
  ) %>%
  mutate(severity = "mixed", .before = gap)
cea_gdp_all <- bind_rows(cea_gdp, cea_mixed)

prob_ce <- cea_gdp_all %>%
  group_by(country, setting, severity, gap) %>%
  summarise(
    p_benefit = mean(daly_averted > 0),
    p_ce_oc   = mean(oc * daly_averted        - incr_cost > 0),
    p_ce_1x   = mean(gdp * daly_averted       - incr_cost > 0),
    p_ce_3x   = mean(3 * gdp * daly_averted   - incr_cost > 0),
    .groups = "drop"
  )
write_tsv(prob_ce, file.path(CFG$dir_tab, "prob_cost_effective.tsv"))

# ---- Monte Carlo error on the headline quantities ---------------------------
# The probabilities above are proportions of CFG$n_psa independent draws, so they
# carry simulation error of their own. This matters where a probability sits near
# a decision boundary: a value of 0.49 at n=1000 cannot be distinguished from
# 0.50. A Wilson score interval is used because the normal approximation gives a
# zero-width interval when a probability reaches 0 or 1, which several do.
wilson <- function(p, n, z = 1.96) {
  denom  <- 1 + z^2 / n
  centre <- (p + z^2 / (2 * n)) / denom
  half   <- z / denom * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))
  c(lo = max(0, centre - half), hi = min(1, centre + half))
}
n_draws <- CFG$n_psa
mc_prob <- prob_ce %>%
  pivot_longer(starts_with("p_"), names_to = "quantity", values_to = "estimate") %>%
  rowwise() %>%
  mutate(
    mc_se    = sqrt(estimate * (1 - estimate) / n_draws),
    mc_lo    = wilson(estimate, n_draws)[["lo"]],
    mc_hi    = wilson(estimate, n_draws)[["hi"]]
  ) %>%
  ungroup()

# Simulation error on the mean incremental outcomes is the usual sd/sqrt(n).
mc_mean <- cea_gdp_all %>%
  group_by(country, setting, severity, gap) %>%
  summarise(across(c(daly_averted, incr_cost),
                   list(estimate = mean,
                        mc_se = ~ sd(.x) / sqrt(n_draws)),
                   .names = "{.col}__{.fn}"),
            .groups = "drop") %>%
  pivot_longer(contains("__"), names_to = c("quantity", ".value"), names_sep = "__") %>%
  mutate(mc_lo = estimate - 1.96 * mc_se, mc_hi = estimate + 1.96 * mc_se)

mc_error <- bind_rows(mc_prob, mc_mean) %>%
  mutate(n_draws = n_draws) %>%
  arrange(country, severity, gap, quantity) %>%
  select(country, setting, severity, gap, quantity, n_draws,
         estimate, mc_se, mc_lo, mc_hi)
write_tsv(mc_error, file.path(CFG$dir_tab, "mc_error.tsv"))

message("[04_cea] Monte Carlo error, headline severe-cohort probabilities:")
print(mc_error %>%
  filter(severity == "severe", quantity == "p_ce_oc",
         (setting == "LMIC" & gap == CFG$gap_intervention$LMIC) |
         (setting == "HIC"  & gap == CFG$gap_intervention$HIC)) %>%
  transmute(country, p_ce_oc = round(estimate, 3), mc_se = round(mc_se, 4),
            mc_95 = sprintf("%.3f-%.3f", mc_lo, mc_hi)))

# Cost-effectiveness acceptability curves (CEAC) for the headline intervention
# scenarios (severe cohort, bridged gap), with willingness-to-pay expressed as a
# multiple of each country's own GDP per capita so the three are comparable.
ceac <- cea_gdp %>%
  filter(severity == "severe",
         (setting == "LMIC" & gap == CFG$gap_intervention$LMIC) |
         (setting == "HIC"  & gap == CFG$gap_intervention$HIC)) %>%
  crossing(wtp_ratio = seq(0, 3, by = 0.1)) %>%
  group_by(country, wtp_ratio) %>%
  summarise(p_ce = mean(wtp_ratio * gdp * daly_averted - incr_cost > 0),
            .groups = "drop")
write_tsv(ceac, file.path(CFG$dir_tab, "ceac.tsv"))

# ---- Mortality table: realised death probabilities by model cycle and state --
mortality <- mort %>%
  rowwise() %>%
  mutate(
    Untreated  = rate_to_prob(bg_rate * CFG$smr$Untreated,  CFG$cycle_years),
    Treated    = rate_to_prob(bg_rate * CFG$smr$Treated,    CFG$cycle_years),
    Controlled = rate_to_prob(bg_rate * CFG$smr$Controlled, CFG$cycle_years),
    Remission  = rate_to_prob(bg_rate * CFG$smr$Remission,  CFG$cycle_years)
  ) %>% ungroup() %>%
  pivot_longer(c(Untreated, Treated, Controlled, Remission),
               names_to = "state", values_to = "p_death_cycle") %>%
  mutate(
    p_death_annual = 1 - (1 - p_death_cycle)^(1 / CFG$cycle_years),
    smr = unlist(CFG$smr[state]),
    state = factor(state, levels = CFG$states)
  ) %>%
  arrange(country, cycle, state) %>%
  select(country, cycle, age_mid, state, smr, bg_rate,
         p_death_annual, p_death_cycle, residual_le)
write_tsv(mortality, file.path(CFG$dir_tab, "mortality_table.tsv"))

message("[04_cea] headline (severe cohort, intervention vs status-quo gap):")
print(cea_summary %>%
  filter(severity == "severe",
         (setting == "LMIC" & gap == CFG$gap_intervention$LMIC) |
         (setting == "HIC"  & gap == CFG$gap_intervention$HIC)) %>%
  transmute(country, gap, daly_averted = round(daly_averted, 3),
            incr_cost = round(incr_cost, 1), icer = round(icer, 1), verdict))
message("[04_cea] realized annual death probability per state:")
print(mortality %>% filter(cycle %in% c(1, CFG$n_cycles)) %>%
      transmute(country, age_mid, state, annual = round(p_death_annual, 4)) %>%
      pivot_wider(names_from = state, values_from = annual), n = 10)
message("[04_cea] tables written to output/tables/")
