# ==============================================================================
# 05_figures_tables.R
#
# Willem M. (Wim) Otte, w.m.otte@umcutrecht.nl
#
# Publication figures (ggplot2) + parameter table.
#
# Outputs (output/figures/, output/diagnostics/, output/tables/):
#   fig1_state_occupancy.png : cohort state-occupancy over time (status-quo vs
#                              bridged gap), Nigeria severe cohort
#   fig2_dalys_averted.png   : per-person DALYs averted vs residual gap, by ctry
#   fig4_ceac.png            : cost-effectiveness acceptability curves
#   output/diagnostics/ce_plane.png : cost-effectiveness plane (PSA cloud) +
#                              thresholds, not a manuscript figure (see below)
#   table_parameters.tsv     : model parameter table with provenance keys
#
# fig* and figS* basenames carry the figure number they hold in the manuscript
# and the appendix respectively. Renumber the figures and rename the files.
# ==============================================================================

if (!exists("CFG")) {
  cand <- c("R/00_config.R", "00_config.R",
            file.path(Sys.getenv("CAE_ROOT", "."), "R/00_config.R"))
  source(cand[file.exists(cand)][1])
}
if (!exists("theme_cae")) source(file.path(CFG$root, "R", "plot_theme.R"))
suppressPackageStartupMessages({ library(tidyverse); library(scales) })

sim <- readRDS(file.path(CFG$dir_derived, "sim_traces.rds"))
per_person <- read_tsv(file.path(CFG$dir_tab, "per_person_dalys.tsv"),
                       show_col_types = FALSE)
cea <- read_tsv(file.path(CFG$dir_tab, "cea_results.tsv"), show_col_types = FALSE)

theme_set(theme_cae(base_size = 12))
pal_state <- state_colours

# ---- Fig 1: state occupancy over time (Nigeria severe; status-quo vs bridged)-
occ <- sim$traces %>%
  filter(country == "Nigeria", severity == "severe",
         gap %in% c(CFG$gap_base$LMIC, CFG$gap_intervention$LMIC)) %>%
  group_by(gap, cycle, state) %>%
  summarise(prob = mean(prob), .groups = "drop") %>%
  mutate(year = cycle * CFG$cycle_years,
         # Status-quo panel LEFT (the baseline the reader starts from), bridged RIGHT.
         scenario = factor(
           if_else(gap == CFG$gap_base$LMIC,
                   sprintf("Status quo (%.0f%% untreated at entry)", 100 * CFG$gap_base$LMIC),
                   sprintf("Bridged (%.0f%% untreated at entry)", 100 * CFG$gap_intervention$LMIC)),
           levels = c(sprintf("Status quo (%.0f%% untreated at entry)", 100 * CFG$gap_base$LMIC),
                      sprintf("Bridged (%.0f%% untreated at entry)", 100 * CFG$gap_intervention$LMIC))))

# Stack states from healthiest (bottom) to sickest (top) is unconventional for an
# epilepsy cohort; we keep clinical reading order (Untreated -> ... -> Death) so the
# growing grey Death band sits where the eye expects accumulating mortality.
occ <- occ %>% mutate(state = factor(state, levels = CFG$states))

fig1 <- ggplot(occ, aes(year, prob, fill = state)) +
  geom_area(alpha = .9) +
  facet_wrap(~ scenario) +
  scale_fill_manual(values = pal_state, name = NULL) +
  scale_y_continuous(labels = percent_format(1)) +
  scale_x_continuous(breaks = seq(0, CFG$horizon_y, 2)) +
  labs(x = "Years from model entry", y = "Share of cohort")
save_fig(fig1, "fig1_state_occupancy", width = 9, height = 4.8)

# ---- Fig 2: DALYs averted vs residual gap ------------------------------------
sev_lab <- c(severe = "Severe", less_severe = "Less severe")
fig2_df <- cea %>%
  filter(daly_averted >= 0 | gap %in% unlist(CFG$gap_intervention)) %>%
  mutate(Severity = factor(sev_lab[severity], levels = c("Severe", "Less severe")))
fig2 <- ggplot(fig2_df, aes(gap, daly_averted, colour = country, linetype = Severity)) +
  geom_ribbon(aes(ymin = da_lo, ymax = da_hi, fill = country),
              alpha = .12, colour = NA, show.legend = FALSE) +
  geom_line(linewidth = .9) + geom_point(size = 2) +
  scale_colour_manual(values = country_colours) +
  scale_fill_manual(values = country_colours) +
  scale_x_reverse(breaks = seq(0, 1, 0.1), labels = percent_format(1)) +
  scale_y_continuous(breaks = number_ticks(6)) +
  labs(x = "Residual treatment gap (% of people with epilepsy still untreated)",
       y = "DALYs averted per person", colour = NULL, fill = NULL, linetype = NULL)
save_fig(fig2, "fig2_dalys_averted", width = 8.5, height = 5.2)

# ---- CE plane: not a manuscript figure ---------------------------------------
# Cut from the manuscript as unintuitive; table 2 carries the same result. Kept
# because the analysis still uses it, and deliberately left outside the numbered
# fig*/figS* sequence so it never collides with a manuscript figure.
params <- sim$params
# rebuild per-iter incremental points for the intervention scenarios
ce_pts <- read_tsv(file.path(CFG$dir_tab, "per_person_dalys.tsv"),
                   show_col_types = FALSE)  # means only; build cloud below
# Use per-iter from a quick recompute saved by 04 if present; else means.
ce <- cea %>%
  filter((setting == "LMIC" & gap == CFG$gap_intervention$LMIC) |
         (setting == "HIC"  & gap == CFG$gap_intervention$HIC)) %>%
  mutate(Severity = factor(sev_lab[severity], levels = c("Severe", "Less severe")),
         country = factor(country, levels = c("Nigeria", "Ecuador", "UK")))

# Per-country threshold lines (own GDP-per-capita) for facetting.
thr_lines <- ce %>% distinct(country, gdp_pc, thr_oc) %>%
  mutate(`Opp. cost (0.5x GDP)` = thr_oc, `1x GDP` = gdp_pc, `3x GDP` = 3 * gdp_pc) %>%
  pivot_longer(c(`Opp. cost (0.5x GDP)`, `1x GDP`, `3x GDP`),
               names_to = "Threshold", values_to = "slope") %>%
  mutate(Threshold = factor(Threshold,
         levels = c("Opp. cost (0.5x GDP)", "1x GDP", "3x GDP")))

ce_plane <- ggplot(ce, aes(daly_averted, incr_cost)) +
  geom_abline(data = thr_lines,
              aes(slope = slope, intercept = 0, linetype = Threshold),
              colour = "gray50") +
  geom_errorbar(aes(ymin = ic_lo, ymax = ic_hi, colour = country), width = 0, alpha = .5,
                show.legend = FALSE) +
  geom_errorbar(aes(xmin = da_lo, xmax = da_hi, colour = country),
                orientation = "y", width = 0, alpha = .5, show.legend = FALSE) +
  geom_point(aes(shape = Severity, colour = country), size = 3) +
  facet_wrap(~ country, scales = "free") +
  scale_colour_manual(values = country_colours) +
  scale_linetype_manual(values = c(`Opp. cost (0.5x GDP)` = "dotted",
                                    `1x GDP` = "dashed", `3x GDP` = "longdash")) +
  scale_x_continuous(breaks = number_ticks(5)) +
  scale_y_continuous(breaks = number_ticks(5)) +
  expand_limits(x = 0, y = 0) +
  guides(colour = "none") +
  labs(x = "DALYs averted per person", y = "Incremental cost per person (USD)",
       shape = NULL, linetype = NULL)
save_fig(ce_plane, "ce_plane", width = 9.5, height = 4.8,
         output_dir = CFG$dir_diag)

# ---- Fig 4: cost-effectiveness acceptability curves (CEAC) -------------------
# Probability bridging is cost-effective as willingness-to-pay rises, expressed
# as a multiple of each country's GDP per capita. Directly answers "the averted
# CI crosses zero" by showing how often bridging actually pays off.
ceac <- read_tsv(file.path(CFG$dir_tab, "ceac.tsv"), show_col_types = FALSE) %>%
  mutate(country = factor(country, levels = c("Nigeria", "Ecuador", "UK")))
# The 0.5x opportunity-cost line is the threshold the paper actually judges
# against, so it has to be on the plot. Linetypes match the CE plane: dotted for
# opportunity cost, dashed for 1x, longdash for 3x.
fig4 <- ggplot(ceac, aes(wtp_ratio, p_ce, colour = country)) +
  geom_vline(xintercept = c(0.5, 1, 3), linetype = c(3, 2, 5), colour = "grey60") +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = country_colours) +
  scale_x_continuous(breaks = c(0, 0.5, 1, 2, 3),
                     labels = c("0", "0.5x", "1x", "2x", "3x")) +
  scale_y_continuous(labels = percent_format(1), limits = c(0, 1)) +
  annotate("text", x = 0.5, y = 0.04, label = "0.5x GDP", size = 3,
           hjust = -0.1, colour = "grey40") +
  labs(x = "Willingness-to-pay per DALY averted (multiple of GDP per capita)",
       y = "Proportion of draws cost-effective", colour = NULL)
save_fig(fig4, "fig4_ceac", width = 8, height = 5)

# ---- Parameter table ---------------------------------------------------------
# Every row carries a Basis, because a triangular limit can be a source
# confidence interval, a prespecified plausible range, a deliberate stress test,
# or an author assumption, and those four have different epistemic weight. The
# earlier table said "values marked as assumptions" and then marked none.
# Abbreviations are spelled out. The appendix table has no room for telegram style.
#
# Formatting notes. make_supplement_tables.py converts a hyphen to an en rule only
# when it sits BETWEEN two digits, so a percentage range ("35%-70%") would keep its
# hyphen and trip the style gate. Percentage ranges therefore use the word "to".
# Currency is pre-formatted with the journal's space thousands separator, because
# the python side does not regroup digits inside free text.
usd <- function(x, digits = 2) {
  sprintf("US$%s", formatC(x, format = "f", digits = digits, big.mark = " "))
}
pc_trim <- function(x) paste0(format(100 * x, drop0trailing = TRUE), "%")
param_tbl <- tribble(
  ~Parameter, ~Value, ~Basis, ~Source,

  "Health states", paste(CFG$states, collapse = ", "),
    "Structural assumption", "Model design",
  "Cycle length", sprintf("%d years", CFG$cycle_years),
    "Structural assumption", "Model design",
  "Cycles and horizon", sprintf("%d cycles, %d years", CFG$n_cycles, CFG$horizon_y),
    "Structural assumption", "Model design",
  "Index age at entry", sprintf("%d years", CFG$index_age),
    "Structural assumption", "Adult index cohort, not an age distribution",
  "Annual discount rate", percent(CFG$discount, .1),
    "Prespecified range 0-6% in one-way analysis",
    "NICE convention, applied to all three settings for comparability",

  "Treatment gap at entry, status quo to scenario",
    sprintf("LMIC %s to %s, HIC %s to %s",
            pc_trim(CFG$gap_base$LMIC), pc_trim(CFG$gap_intervention$LMIC),
            pc_trim(CFG$gap_base$HIC), pc_trim(CFG$gap_intervention$HIC)),
    "Fixed scenario endpoint, not varied",
    "Scenario definition, not measured national coverage",

  "Disability weight, untreated severe epilepsy", sprintf("%.3f", CFG$dw$Untreated_severe),
    sprintf("Source uncertainty interval %.3f-%.3f", dw_bounds("Untreated_severe")[1],
            dw_bounds("Untreated_severe")[2]),
    "Salomon et al., 2015; defined as at least one seizure per month",
  "Disability weight, untreated less-severe epilepsy",
    sprintf("%.3f", CFG$dw$Untreated_less_severe),
    sprintf("Source uncertainty interval %.3f-%.3f", dw_bounds("Untreated_less_severe")[1],
            dw_bounds("Untreated_less_severe")[2]),
    "Salomon et al., 2015; defined as one to eleven seizures per year",
  "Disability weight, treated with continuing seizures", sprintf("%.3f", CFG$dw$Treated),
    "Structural assumption, sharing the less-severe draw",
    "Not a source estimate for treated epilepsy with ongoing seizures",
  "Disability weight, controlled and remission", sprintf("%.3f", CFG$dw$Controlled),
    sprintf("Source uncertainty interval %.3f-%.3f", dw_bounds("Controlled")[1],
            dw_bounds("Controlled")[2]),
    "Salomon et al., 2015; treated, seizure-free",

  "Standardised mortality ratio, untreated", sprintf("%.2f", CFG$smr$Untreated),
    sprintf("Source 95%% confidence interval %.1f-%.1f", CFG$smr_sens$Untreated[1],
            CFG$smr_sens$Untreated[2]),
    "Carpio et al., 2005; newly diagnosed Ecuador cohort, not an untreated subgroup",
  "Standardised mortality ratio, treated with continuing seizures",
    sprintf("%.2f", CFG$smr$Treated),
    sprintf("Source 95%% confidence interval %.2f-%.2f", CFG$smr_sens$Treated[1],
            CFG$smr_sens$Treated[2]),
    "Mohanraj et al., 2006; newly diagnosed cohort, no response to treatment",
  "Standardised mortality ratio, controlled and remission",
    sprintf("%.2f", CFG$smr$Controlled),
    sprintf("Source 95%% confidence interval %.2f-%.2f", CFG$smr_sens$SeizureFree[1],
            CFG$smr_sens$SeizureFree[2]),
    "Mohanraj et al., 2006; newly diagnosed cohort, seizure-free on treatment",

  "Spontaneous remission from untreated",
    sprintf("%s per cycle (about 30%% cumulative over 10 years)",
            percent(CFG$spont_remission_prob, .1)),
    sprintf("Prespecified range %s to %s per cycle",
            percent(CFG$spont_remission_sens[1], .1),
            percent(CFG$spont_remission_sens[2], .1)),
    "Nicoletti et al., 2009, converted from a cumulative outcome",
  "Cumulative seizure control by 10 years, Nigeria and Ecuador",
    percent(CFG$control_cum_10y$LMIC),
    sprintf("Prespecified range %s to %s", percent(CFG$control_cum_sens$LMIC[1]),
            percent(CFG$control_cum_sens$LMIC[2])),
    sprintf("Scenario input, converted to %s per 2-year cycle",
            percent(CFG$control_prob$LMIC, .1)),
  "Cumulative seizure control by 10 years, UK", percent(CFG$control_cum_10y$HIC),
    sprintf("Prespecified range %s to %s", percent(CFG$control_cum_sens$HIC[1]),
            percent(CFG$control_cum_sens$HIC[2])),
    sprintf("Scenario input, converted to %s per 2-year cycle",
            percent(CFG$control_prob$HIC, .1)),
  "Composite loss-to-care and recurrence, scenario only",
    percent(CFG$relapse_scenario), "Structural scenario, base case is zero",
    "Mixes recurrence, discontinuation, and loss to care in one route",

  "Severe share of the epilepsy population", percent(CFG$severity_frac[["severe"]], .1),
    "Fixed; 9.2%, 19.0%, and 36.3% in scenarios",
    "Sultana et al., 2021, population-based prevalence of drug resistance, used as a proxy for higher disability",

  "Annual direct health-care cost, Nigeria", usd(CFG$cost_annual[["Nigeria"]]),
    "Prespecified multiplier 0.75-1.25", "Begley et al., 2022, in 2019 US dollars",
  "Annual direct health-care cost, Ecuador", usd(CFG$cost_annual[["Ecuador"]]),
    "Prespecified multiplier 0.75-1.25",
    sprintf(paste("Begley et al., 2022, lower-middle-income group mean in 2019 US",
                  "dollars. Ecuador is upper-middle income, whose group mean is %s"),
            usd(CFG$cost_annual_ecu_umic)),
  "Annual direct health-care cost, UK", usd(CFG$cost_annual[["UK"]]),
    "Prespecified multiplier 0.75-1.25", "Begley et al., 2022, in 2019 US dollars",
  "Untreated care cost, share of treated care", percent(CFG$cost_untreated_frac),
    sprintf("Sampled %s to %s, stressed to %s in the one-way analysis",
            percent(CFG$cost_untreated_frac_psa[1]),
            percent(CFG$cost_untreated_frac_psa[2]),
            percent(CFG$cost_untreated_frac_sens[2])),
    "Author assumption. The upper bound is a stress test, not a likely value",
  "Case-finding and initiation cost, per additional person treated",
    usd(CFG$cost_case_finding, 0),
    sprintf("Prespecified range %s to %s", usd(CFG$cost_case_finding_sens[1], 0),
            usd(CFG$cost_case_finding_sens[2], 0)),
    "Author assumption, covering identification and initiation only",

  "Prevalence per 1000, Nigeria", "9.8 (8.6-11.1)", "Source range, not propagated",
    "Watila et al., 2021",
  "Prevalence per 1000, Ecuador", "7.5 (6.0-9.0)", "Source range, not propagated",
    "Placencia et al. and the GBD band",
  "Prevalence per 1000, UK", "9.37 (8.85-12.08)", "Source range, not propagated",
    "Wigglesworth et al., 2023. The band spans the four UK nations",
  "National population, Nigeria, Ecuador, UK",
    sprintf("%s / %s / %s",
            formatC(CFG$population[["Nigeria"]], format = "d", big.mark = " "),
            formatC(CFG$population[["Ecuador"]], format = "d", big.mark = " "),
            formatC(CFG$population[["UK"]], format = "d", big.mark = " ")),
    "Fixed, not propagated",
    "Approximate 2024 national totals, used only as a scaling input",

  "GDP per capita, Nigeria, Ecuador, UK",
    sprintf("%s / %s / %s", usd(CFG$gdp_pc[["Nigeria"]], 0),
            usd(CFG$gdp_pc[["Ecuador"]], 0), usd(CFG$gdp_pc[["UK"]], 0)),
    "Fixed, not varied", "World Bank, current US dollars",
  "Threshold proxy per DALY averted, Nigeria, Ecuador, UK",
    sprintf("%s / %s / %s", usd(CFG$cet_oc[["Nigeria"]], 0),
            usd(CFG$cet_oc[["Ecuador"]], 0), usd(CFG$cet_oc[["UK"]], 0)),
    "Fixed, not varied",
    sprintf(paste("%.1f times GDP per capita, a prespecified proxy and not a",
                  "country-specific opportunity cost (Woods 2016, Ochalek 2018)"),
            CFG$cet_oc_frac),

  "Background mortality and residual life expectancy", "Age-varying",
    "Fixed", "UN World Population Prospects 2024, interpolated at cycle midpoints",
  "Probabilistic distributions", "Triangular",
    "See the Basis column for each parameter",
    "Each sampled set is applied to both coverage scenarios",
  "Probabilistic iterations", as.character(CFG$n_psa), "Model design",
    "Monte Carlo error is reported in table S9"
)
write_tsv(param_tbl, file.path(CFG$dir_tab, "table_parameters.tsv"))

message("[05_figures] wrote fig1, fig2, fig4 (CEAC), ce_plane + parameter table")
