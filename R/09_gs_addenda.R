# ==============================================================================
# 09_gs_addenda.R
#
# Willem M. (Wim) Otte, w.m.otte@umcutrecht.nl
#
# Addenda requested in Gagandeep Singh's 21-Jun-2026 review round. Reads existing
# model outputs only (no re-simulation) and emits new tables/figures:
#
#   population_daly_reduction.tsv : baseline vs bridged TOTAL epilepsy DALYs and
#                                   % reduction per country (Singh C502, C315).
#                                   Model-internal 10-year-horizon totals; the
#                                   GBD-2021 annual denominator is a separate,
#                                   still-to-be-reconciled framing (see memo).
#   state_occupancy_snapshots.tsv : mean state occupancy per cycle, status-quo vs
#                                   bridged, all countries/cohorts (C513, C596, C901).
#   figS4_combined_dalys.png      : combined (severe+less-severe) DALYs averted
#                                   per person across the gap sweep, by country (C315).
#   figS5_cohort_ecuador.png       : Ecuador state occupancy, both strata.
#   figS6_cohort_uk.png            : UK state occupancy, both strata.
# ==============================================================================

if (!exists("CFG")) {
  cand <- c("R/00_config.R", "00_config.R",
            file.path(Sys.getenv("CAE_ROOT", "."), "R/00_config.R"))
  source(cand[file.exists(cand)][1])
}
if (!exists("theme_cae")) source(file.path(CFG$root, "R", "plot_theme.R"))
suppressPackageStartupMessages({ library(tidyverse); library(scales) })

sim        <- readRDS(file.path(CFG$dir_derived, "sim_traces.rds"))
per_person <- read_tsv(file.path(CFG$dir_tab, "per_person_dalys.tsv"), show_col_types = FALSE)
prev       <- read_tsv(file.path(CFG$dir_derived, "prevalence.tsv"), show_col_types = FALSE)
cea        <- read_tsv(file.path(CFG$dir_tab, "cea_results.tsv"), show_col_types = FALSE)
psa        <- read_tsv(file.path(CFG$dir_tab, "psa_draws_bridged.tsv"), show_col_types = FALSE)
cea_draws  <- read_tsv(file.path(CFG$dir_tab, "cea_draws.tsv"), show_col_types = FALSE)

base_gap_of <- function(setting) CFG$gap_base[[setting]]
itv_gap_of  <- function(setting) CFG$gap_intervention[[setting]]
sev_w <- tibble(severity = names(CFG$severity_frac), w = as.numeric(CFG$severity_frac))

# ---- (1) Population DALY reduction: baseline vs bridged (C502, C315) ----------
# Each cohort's per-person 10y DALYs at the status-quo gap and at the bridged gap,
# scaled by its severity-specific slice of the national PWE population, summed.
pp <- per_person %>%
  left_join(prev %>% select(country, setting, prevalence_per_1000, population),
            by = c("country", "setting")) %>%
  left_join(sev_w, by = "severity") %>%
  mutate(pwe = population * prevalence_per_1000 / 1000 * w)

reduction <- pp %>%
  group_by(country, setting) %>%
  summarise(
    pwe_total   = sum(pwe[gap == itv_gap_of(setting[1])]),  # same total either gap
    daly_base   = sum(daly_mean[gap == base_gap_of(setting[1])] * pwe[gap == base_gap_of(setting[1])]),
    daly_bridged= sum(daly_mean[gap == itv_gap_of(setting[1])]  * pwe[gap == itv_gap_of(setting[1])]),
    .groups = "drop"
  ) %>%
  mutate(daly_reduction = daly_base - daly_bridged,
         pct_reduction  = 100 * daly_reduction / daly_base)
write_tsv(reduction, file.path(CFG$dir_tab, "population_daly_reduction.tsv"))

# A previous version divided discounted 10-year cohort DALYs by ten and compared
# the result with an annual GBD estimate. Those quantities are not commensurate:
# the model has no incident replenishment and discounts lifetime YLL at death.
# The comparison is therefore retired rather than cosmetically relabelled.
legacy_gbd <- file.path(CFG$dir_tab, "gbd_denominator_reconciliation.tsv")
if (file.exists(legacy_gbd)) unlink(legacy_gbd)

# ---- (2) State-occupancy snapshots (C513, C596, C901) ------------------------
snap <- sim$traces %>%
  filter((setting == "LMIC" & gap %in% c(CFG$gap_base$LMIC, CFG$gap_intervention$LMIC)) |
         (setting == "HIC"  & gap %in% c(CFG$gap_base$HIC,  CFG$gap_intervention$HIC))) %>%
  group_by(country, setting, severity, gap, cycle, state) %>%
  summarise(prob = mean(prob), .groups = "drop") %>%
  mutate(year = cycle * CFG$cycle_years,
         scenario = if_else(gap %in% c(CFG$gap_base$LMIC, CFG$gap_base$HIC),
                            "status_quo", "bridged")) %>%
  arrange(country, severity, scenario, cycle, state)
write_tsv(snap, file.path(CFG$dir_tab, "state_occupancy_snapshots.tsv"))

# ---- (3) Fig 7: combined (both cohorts) DALYs averted per person, by country --
# Weighted average of the two cohorts' per-person averted DALYs (severity_frac).
comb_draws <- cea_draws %>%
  left_join(sev_w, by = "severity") %>%
  group_by(country, setting, gap, iter) %>%
  summarise(daly_averted_pp = weighted.mean(daly_averted, w), .groups = "drop")
comb <- comb_draws %>%
  group_by(country, setting, gap) %>%
  summarise(daly_averted_pp = mean(daly_averted_pp),
            da_lo = quantile(daly_averted_pp, 0.025),
            da_hi = quantile(daly_averted_pp, 0.975), .groups = "drop") %>%
  mutate(country = factor(country, levels = c("Nigeria", "Ecuador", "UK")))
figS4 <- ggplot(comb, aes(gap, daly_averted_pp, colour = country)) +
  geom_ribbon(aes(ymin = da_lo, ymax = da_hi, fill = country), alpha = .12,
              colour = NA, show.legend = FALSE) +
  geom_line(linewidth = 1) + geom_point(size = 2.2) +
  scale_colour_manual(values = country_colours) +
  scale_fill_manual(values = country_colours) +
  scale_x_reverse(breaks = seq(0, 1, 0.1), labels = percent_format(1)) +
  scale_y_continuous(breaks = number_ticks(6)) +
  labs(x = "Residual treatment gap (% of people with epilepsy still untreated)",
       y = "DALYs averted per person (all epilepsy, both cohorts)",
       colour = NULL) +
  theme_cae(base_size = 12)
save_fig(figS4, "figS4_combined_dalys", width = 8.5, height = 5.2)

# ---- (4) Cohort evolution for Ecuador and the UK (Singh follow-up) ------------
# State occupancy is severity-independent because severity changes disability
# weights only, not transitions or mortality. Main figure 1 therefore already
# represents both Nigerian strata, and one panel pair per remaining country is
# sufficient.
pal_state <- state_colours
# Status-quo/bridged labels are read from the setting-specific gaps.
cohort_fig <- function(cty) {
  st       <- CFG$countries$setting[CFG$countries$country == cty]
  pc <- function(x) paste0(format(100 * x, drop0trailing = TRUE), "%")
  base_lab <- sprintf("Status quo (%s untreated at entry)", pc(CFG$gap_base[[st]]))
  brid_lab <- sprintf("Reduced gap (%s untreated at entry)", pc(CFG$gap_intervention[[st]]))
  occ <- snap %>%
    filter(country == cty, severity == "less_severe") %>%
    mutate(state = factor(state, levels = CFG$states),
           scenario = factor(if_else(scenario == "status_quo", base_lab, brid_lab),
                             levels = c(base_lab, brid_lab)))
  ggplot(occ, aes(year, prob, fill = state)) +
    geom_area(alpha = .9) +
    facet_wrap(~ scenario) +
    scale_fill_manual(values = pal_state, name = NULL) +
    scale_y_continuous(labels = percent_format(1)) +
    scale_x_continuous(breaks = seq(0, CFG$horizon_y, 2)) +
    labs(x = "Years from model entry", y = "Share of cohort") +
    theme_cae(base_size = 12)
}
save_fig(cohort_fig("Ecuador"), "figS5_cohort_ecuador", width = 9, height = 4.8)
save_fig(cohort_fig("UK"), "figS6_cohort_uk", width = 9, height = 4.8)

# ---- (5) National cohort sizes in absolute numbers (Singh follow-up) ----------
# Report the modelled cohort as counts, not only percentages: national people
# with epilepsy split into the severe (13.7%) and less-severe (86.3%) strata.
cohort_sizes <- prev %>%
  distinct(country, population, prevalence_per_1000) %>%
  transmute(country,
            pwe_total     = population * prevalence_per_1000 / 1000,
            n_severe      = pwe_total * CFG$severity_frac[["severe"]],
            n_less_severe = pwe_total * CFG$severity_frac[["less_severe"]]) %>%
  mutate(across(c(pwe_total, n_severe, n_less_severe), round))
write_tsv(cohort_sizes, file.path(CFG$dir_tab, "cohort_sizes.tsv"))

# ---- (6) Cost-effectiveness-plane quadrant distribution (Singh follow-up) -----
# Classify each of the 1000 PSA draws (severe cohort, fully bridged gap) by the
# sign of incremental cost and averted DALYs. The dominant/cost-saving quadrant
# (lower cost AND more health) answers "is bridging ever cost-saving?".
ce_quadrants <- psa %>%
  filter(severity == "severe") %>%
  mutate(quadrant = case_when(
    incr_cost >  0 & daly_averted >  0 ~ "NE_costlier_moreDALYs",
    incr_cost <= 0 & daly_averted >  0 ~ "SE_dominant_costsaving",
    incr_cost >  0 & daly_averted <= 0 ~ "NW_dominated",
    TRUE                               ~ "SW_costsaving_fewerDALYs")) %>%
  count(country, quadrant) %>%
  group_by(country) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup() %>%
  select(-n) %>%
  pivot_wider(names_from = quadrant, values_from = pct, values_fill = 0)
# Guarantee all four columns exist even if a quadrant is empty in every country.
for (q in c("NE_costlier_moreDALYs", "SE_dominant_costsaving",
            "NW_dominated", "SW_costsaving_fewerDALYs")) {
  if (!q %in% names(ce_quadrants)) ce_quadrants[[q]] <- 0
}
ce_quadrants <- ce_quadrants %>%
  select(country, NE_costlier_moreDALYs, SE_dominant_costsaving,
         NW_dominated, SW_costsaving_fewerDALYs) %>%
  mutate(across(-country, ~ round(.x, 1)))
write_tsv(ce_quadrants, file.path(CFG$dir_tab, "ce_quadrants.tsv"))

# ---- Console readouts to fill the manuscript blanks --------------------------
message("[09_addenda] population DALY reduction (model 10y-horizon totals):")
print(reduction %>% transmute(country,
        daly_base = round(daly_base), daly_bridged = round(daly_bridged),
        daly_reduction = round(daly_reduction), pct = round(pct_reduction, 1)))
message("[09_addenda] spontaneous-remission occupancy at year 10, status quo:")
print(snap %>% filter(scenario == "status_quo", cycle == CFG$n_cycles,
                      state == "Remission") %>%
        transmute(country, severity, remission_pct = round(100 * prob, 1)))
message("[09_addenda] state occupancy at years 4 & 10 (Nigeria severe, both scenarios):")
print(snap %>% filter(country == "Nigeria", severity == "severe",
                      year %in% c(4, 10)) %>%
        transmute(scenario, year, state, pct = round(100 * prob, 1)) %>%
        pivot_wider(names_from = state, values_from = pct))
message("[09_addenda] national cohort sizes (absolute):")
print(cohort_sizes)
message("[09_addenda] CE-plane quadrant distribution (% of 1000 draws, severe, bridged):")
print(ce_quadrants)
message("[09_addenda] wrote population_daly_reduction, state_occupancy_snapshots, cohort_sizes, ce_quadrants, figS4, figS5, figS6")
