# ==============================================================================
# 07_sensitivity.R
#
# Willem M. (Wim) Otte, w.m.otte@umcutrecht.nl
#
# Deterministic sensitivity analyses around the headline ICER (bridging the
# treatment gap, Nigeria severe cohort: status-quo gap -> intervention gap).
#
#   (a) one-way analysis (tornado) over key parameters at their plausible bounds,
#       including the untreated-care cost fraction, the only input whose range
#       could flip the sign of the incremental cost
#   (b) threshold analysis: ICER vs cumulative 10-year seizure control
#   (c) ICER vs discount rate
#
# Outputs (output/tables/, output/figures/):
#   sensitivity_oneway.tsv  : low/high ICER per parameter (sorted by swing)
#   sensitivity_threshold.tsv
#   fig3_tornado.png (manuscript figure 3)
#   figS3_threshold_control.png (appendix figure S3)
# ==============================================================================

if (!exists("CFG")) {
  cand <- c("R/00_config.R", "00_config.R",
            file.path(Sys.getenv("CAE_ROOT", "."), "R/00_config.R"))
  source(cand[file.exists(cand)][1])
}
source(file.path(CFG$root, "R", "02_markov_model.R"))
source(file.path(CFG$root, "R", "lib_cea.R"))   # shared trace-reward functions
if (!exists("theme_cae")) source(file.path(CFG$root, "R", "plot_theme.R"))
suppressPackageStartupMessages({ library(tidyverse); library(scales) })

mort <- read_tsv(file.path(CFG$dir_derived, "mortality_schedule.tsv"),
                 show_col_types = FALSE)

# ---- Base parameter set for the focal scenario (Nigeria severe) --------------
focal_country  <- "Nigeria"
bg_rates_focal <- mort %>% filter(country == focal_country) %>% arrange(cycle) %>% pull(bg_rate)
resid_le_focal <- mort %>% filter(country == focal_country) %>% arrange(cycle) %>% pull(residual_le)
cost_focal     <- CFG$cost_annual[[focal_country]]

base <- list(
  bg_rates     = bg_rates_focal,
  smr_U = CFG$smr$Untreated, smr_T = CFG$smr$Treated,
  smr_C = CFG$smr$Controlled, smr_R = CFG$smr$Remission,
  spont_rem    = CFG$spont_remission_prob,
  control_prob = CFG$control_prob$LMIC,
  # Treated keeps the 0.263 band in both cohorts, as in 04's dw_vec().
  dw_U  = CFG$dw$Untreated_severe, dw_T = CFG$dw$Treated,
  dw_C  = CFG$dw$Controlled, dw_R = CFG$dw$Remission,
  discount = CFG$discount, cost_yr = cost_focal, resid_le = resid_le_focal,
  relapse = CFG$relapse_base, cost_cf = CFG$cost_case_finding,
  untr_frac = CFG$cost_untreated_frac
)

# Deterministic DALY + cost for one gap, given a parameter list `q`.
daly_cost <- function(q, gap) {
  p <- list(bg_rates = q$bg_rates, bg_rate = q$bg_rates[1], initial_gap = gap,
            smr = list(Untreated = q$smr_U, Treated = q$smr_T,
                       Controlled = q$smr_C, Remission = q$smr_R),
            treat_uptake = CFG$background_treatment_uptake, spont_rem = q$spont_rem,
            control_prob = q$control_prob,
            relapse = if (is.null(q$relapse)) 0 else q$relapse)
  full <- run_cohort(p)
  tr <- full %>% select(cycle, state, prob) %>%
    pivot_wider(names_from = state, values_from = prob) %>% arrange(cycle)
  tt <- (0:CFG$n_cycles) * CFG$cycle_years
  dw <- c(Untreated = q$dw_U, Treated = q$dw_T, Controlled = q$dw_C,
          Remission = q$dw_R, Death = 0)
  ann <- discounted_annuity(q$resid_le, q$discount)
  uf <- if (is.null(q$untr_frac)) CFG$cost_untreated_frac else q$untr_frac
  cf <- if (is.null(q$cost_cf)) CFG$cost_case_finding else q$cost_cf
  value_trace(
    tt, tr, dw, ann, q$cost_yr, uf, CFG$gap_base$LMIC, gap, cf, q$discount
  )[c("daly", "cost")]
}

# ICER for bridging (status-quo gap -> intervention gap) under parameter list q.
icer_of <- function(q) {
  sq  <- daly_cost(q, CFG$gap_base$LMIC)
  itv <- daly_cost(q, CFG$gap_intervention$LMIC)
  inc <- incremental_cea(sq[["daly"]], sq[["cost"]],
                         itv[["daly"]], itv[["cost"]],
                         CFG$cet_oc[[focal_country]])
  c(icer = inc$icer, averted = inc$daly_averted,
    incr_cost = inc$incr_cost, inb = inc$inb)
}

base_icer <- icer_of(base)[["icer"]]
base_res <- icer_of(base)
base_inb <- base_res[["inb"]]
message(sprintf("[07_sens] base-case ICER (Nigeria severe) = $%.0f / DALY", base_icer))

# ---- (a) One-way / tornado ---------------------------------------------------
# Each row: parameter, the field(s) it sets, low value, high value, label.
#
# `fields` is a list-column of character vectors, because some rows move more
# than one field. Controlled and Remission are one seizure-free construct and
# share a single PSA draw (`smr_SF`). An earlier version of this table moved one
# half of a coupled SMR pair alone, so the tornado silently analysed a model the
# headline never ran. Any parameter the PSA couples must be coupled here too.
#
# Labels are submission-facing: they name the parameter in full and say what the
# limits are (source interval, prespecified plausible range, or stress test).
vary <- tribble(
  ~param,                                     ~fields,          ~low,   ~high,  ~label,
  "Cumulative 10-year seizure control",       list("control_prob"),
    1-(1-CFG$control_cum_sens$LMIC[1])^(1/CFG$n_cycles),
    1-(1-CFG$control_cum_sens$LMIC[2])^(1/CFG$n_cycles),
    "Cumulative 10-year seizure control 35-70% (prespecified range)",
  "SMR, untreated",                           list("smr_U"),     2.0,   10.0,
    "SMR assigned to untreated 2.0-10.0 (Carpio 95% CI)",
  "SMR, treated with continuing seizures",    list("smr_T"),
    CFG$smr_sens$Treated[1], CFG$smr_sens$Treated[2],
    "SMR treated with continuing seizures 1.84-3.44 (Mohanraj 95% CI)",
  "SMR, seizure-free (controlled and remission)", list(c("smr_C", "smr_R")),
    CFG$smr_sens$SeizureFree[1], CFG$smr_sens$SeizureFree[2],
    "SMR seizure-free, controlled and remission 0.68-1.29 (Mohanraj 95% CI)",
  "Spontaneous remission from untreated",     list("spont_rem"),
    CFG$spont_remission_sens[1], CFG$spont_remission_sens[2],
    sprintf("Spontaneous remission %.1f-%.1f%% per cycle (prespecified range)",
            100*CFG$spont_remission_sens[1], 100*CFG$spont_remission_sens[2]),
  "Annual discount rate",                     list("discount"),  0.0,   0.06,
    "Annual discount rate 0-6% (prespecified range)",
  "Disability weight, severe untreated",      list("dw_U"),      0.375, 0.710,
    "Disability weight, severe untreated 0.375-0.710 (GBD uncertainty interval)",
  "Annual health-care cost",                  list("cost_yr"),
    cost_focal*0.75, cost_focal*1.25,
    "Annual health-care cost 75-125% of base value (prespecified range)",
  "Case-finding and initiation cost",         list("cost_cf"),
    CFG$cost_case_finding_sens[1], CFG$cost_case_finding_sens[2],
    sprintf("Case-finding and initiation cost US$%.0f-%.0f per person (assumption)",
            CFG$cost_case_finding_sens[1], CFG$cost_case_finding_sens[2]),
  "Untreated care cost, share of treated",    list("untr_frac"),
    CFG$cost_untreated_frac_sens[1], CFG$cost_untreated_frac_sens[2],
    sprintf("Untreated care %.0f-%.0f%% of treated care cost (stress test)",
            100*CFG$cost_untreated_frac_sens[1], 100*CFG$cost_untreated_frac_sens[2])
)
# Set every field in `fields` to `val`, so coupled constructs move together.
# `fields` arrives either as a bare string (the threshold and relapse callers) or,
# under rowwise(), as a length-1 list wrapping a character vector. unlist() flattens
# both; without it `q[[c("smr_T","smr_C")]]` would recursively index instead of
# assigning two fields.
set_field <- function(q, fields, val) {
  for (f in unlist(fields)) q[[f]] <- val
  q
}

# Guard the coupling. Controlled and Remission share a single PSA draw
# (`smr_SF`), so the tornado row that moves their SMR must move both fields.
# Moving one alone silently analyses a model the headline never ran.
.smr_sf_fields <- unlist(vary$fields[vary$param == "SMR, seizure-free (controlled and remission)"])
stopifnot(setequal(.smr_sf_fields, c("smr_C", "smr_R")))

oneway <- vary %>%
  rowwise() %>%
  mutate(
    result_low = list(icer_of(set_field(base, fields, low))),
    result_high = list(icer_of(set_field(base, fields, high))),
    icer_low  = result_low[["icer"]],
    icer_high = result_high[["icer"]],
    inb_low = CFG$cet_oc[[focal_country]] * result_low[["averted"]] - result_low[["incr_cost"]],
    inb_high = CFG$cet_oc[[focal_country]] * result_high[["averted"]] - result_high[["incr_cost"]]
  ) %>%
  ungroup() %>%
  mutate(swing = abs(inb_high - inb_low),
         inb_min = pmin(inb_low, inb_high),
         inb_max = pmax(inb_low, inb_high)) %>%
  arrange(desc(swing))
write_tsv(oneway %>% select(param, label, low, high, inb_low, inb_high,
                            icer_low, icer_high, swing),
          file.path(CFG$dir_tab, "sensitivity_oneway.tsv"))

# Tornado plot
tor <- oneway %>% mutate(param = factor(param, levels = rev(param)))
ny <- length(levels(tor$param))

# The "cost-saving (dominant)" note labels the region left of $0. Draw it only
# when a bar actually reaches there. Otherwise it points at an empty region and,
# because it is right-aligned on a negative x, hangs off the left panel edge and
# is clipped. When it is drawn, the panel needs extra room on the left to hold
# it, which the default 5% expansion does not give.
fig3 <- ggplot(tor) +
  geom_rect(aes(xmin = inb_min, xmax = inb_max,
                ymin = as.integer(param) - .4, ymax = as.integer(param) + .4),
            fill = unname(bc_palette["skyblue"]), colour = "grey30") +
  geom_vline(xintercept = 0, linetype = 3, colour = "gray50") +
  geom_vline(xintercept = base_inb, linetype = 2, colour = "gray40") +
  scale_y_continuous(breaks = seq_along(levels(tor$param)),
                     labels = levels(tor$param),
                     expand = expansion(add = c(.6, .6))) +
  scale_x_continuous(labels = dollar_format(), breaks = number_ticks(6),
                     expand = expansion(mult = c(.05, .05))) +
  labs(x = "Incremental net benefit at 0.5x GDP per capita (USD per person)", y = NULL) +
  theme_cae(base_size = 11)
save_fig(fig3, "fig3_tornado", width = 8.5, height = 4.7)

# ---- (b) Threshold: ICER vs cumulative 10-year seizure control ---------------
thr <- tibble(control_cum = seq(0.30, 0.75, by = 0.05)) %>%
  rowwise() %>%
  mutate(control_prob = 1 - (1 - control_cum)^(1 / CFG$n_cycles),
         res = list(icer_of(set_field(base, "control_prob", control_prob)))) %>%
  mutate(icer = res[["icer"]], averted = res[["averted"]]) %>%
  ungroup() %>% select(control_cum, control_prob, icer, averted)
write_tsv(thr, file.path(CFG$dir_tab, "sensitivity_threshold.tsv"))

base_cp <- CFG$control_cum_10y$LMIC
gdp_ng  <- CFG$gdp_pc[["Nigeria"]]
oc_ng   <- CFG$cet_oc[["Nigeria"]]

# Redesigned Figure 5. The old version was two flat threshold lines far above a
# near-flat ICER curve - visually empty. Instead we answer the question directly
# in two stacked panels that share the x-axis (seizure-control effectiveness):
#   (A) the health gain (DALYs averted/person) it buys, and
#   (B) the incremental NET MONETARY BENEFIT at three thresholds, with the
#       cost-effective region (net benefit > 0) shaded. INB = WTP*averted - cost,
#       and cost = icer*averted, so INB = averted*(WTP - icer).
pa <- "A  Health gain"
pb <- "B  Value for money (net monetary benefit)"
thr_lvls <- c("Opportunity cost (0.5x GDP)", "1x GDP", "3x GDP")

panelA <- thr %>% transmute(control_prob = control_cum, panel = pa, series = "DALYs averted / person",
                            value = averted)
panelB <- thr %>%
  transmute(control_prob = control_cum, icer, averted,
            `Opportunity cost (0.5x GDP)` = averted * (oc_ng     - icer),
            `1x GDP`                      = averted * (gdp_ng    - icer),
            `3x GDP`                      = averted * (3 * gdp_ng - icer)) %>%
  pivot_longer(all_of(thr_lvls), names_to = "series", values_to = "value") %>%
  mutate(panel = pb, series = factor(series, levels = thr_lvls)) %>%
  select(control_prob, panel, series, value)
figS3_df <- bind_rows(panelA, mutate(panelB, series = as.character(series)))

ce_zone <- tibble(panel = pb, xmin = -Inf, xmax = Inf, ymin = 0, ymax = Inf)
zero_ln <- tibble(panel = pb, y = 0)
series_cols <- c("DALYs averted / person"       = unname(bc_palette["green"]),
                 "Opportunity cost (0.5x GDP)"  = unname(bc_palette["orange"]),
                 "1x GDP"                        = unname(bc_palette["skyblue"]),
                 "3x GDP"                        = unname(bc_palette["blue"]))

figS3 <- ggplot(figS3_df, aes(control_prob, value, colour = series)) +
  geom_rect(data = ce_zone, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = unname(bc_palette["green"]), alpha = 0.10) +
  geom_hline(data = zero_ln, aes(yintercept = y), linetype = 2, colour = "gray50") +
  geom_vline(xintercept = base_cp, linetype = 2, colour = "gray50") +
  geom_line(linewidth = 1) +
  geom_point(size = 1.8) +
  annotate("text", x = base_cp, y = -Inf, label = "base case (50%)",
           angle = 90, vjust = -0.4, hjust = -0.05, size = 2.9, colour = "grey40") +
  facet_wrap(~panel, ncol = 1, scales = "free_y", strip.position = "top") +
  scale_x_continuous(labels = percent_format(1),
                     expand = expansion(mult = c(0.02, 0.04))) +
  scale_y_continuous(labels = label_number(big.mark = " "), breaks = number_ticks(6)) +
  scale_colour_manual(values = series_cols, name = NULL) +
  labs(x = "Cumulative seizure control by 10 years in LMICs",
       y = NULL) +
  theme_cae(base_size = 11) +
  theme(strip.text = element_text(face = "bold", hjust = 0))
save_fig(figS3, "figS3_threshold_control", width = 8.0, height = 6.4)

# ---- (c) Relapse / secondary-gap scenario (probabilistic) --------------------
# Bound the optimism of the no-relapse structure: re-compute the headline ICER
# with a per-cycle relapse to Untreated from Treated/Controlled/Remission.
# Each scenario is run deterministically (point estimate) AND over the same
# 1,000 PSA draws used in the main analysis, so the table carries 95% CIs.
psa_params <- readRDS(file.path(CFG$dir_derived, "sim_traces.rds"))$params

# Map one row of sampled PSA parameters onto the focal (Nigeria severe) param
# list, for a given per-cycle relapse probability.
focal_iter <- function(row, relapse) {
  modifyList(base, list(
    smr_U = row$smr_U, smr_T = row$smr_T, smr_C = row$smr_SF, smr_R = row$smr_SF,
    spont_rem = row$spont_rem,
    control_prob = 1 - (1 - row$control_cum_LMIC)^(1 / CFG$n_cycles),
    dw_U = row$dw_sev_UT, dw_T = row$dw_less_UT,
    dw_C = row$dw_control, dw_R = row$dw_control,
    cost_yr = cost_focal * row$cost_mult,
    cost_cf = row$case_find_cost,
    untr_frac = row$untreated_cost_frac,
    relapse = relapse))
}
q025 <- function(x) unname(quantile(x, 0.025, na.rm = TRUE))
q975 <- function(x) unname(quantile(x, 0.975, na.rm = TRUE))

# Singh C898: the relapse ICER barely moves because numerator and denominator
# shrink together; the more informative question is what SHARE of PSA draws stay
# below each threshold. We add the probability cost-effective at the
# opportunity-cost, 1x and 3x GDP thresholds (Nigeria) for each relapse level.
gdp_ng_r <- CFG$gdp_pc[["Nigeria"]]
oc_ng_r  <- CFG$cet_oc[["Nigeria"]]
relapse_tbl <- map_dfr(c(CFG$relapse_base, CFG$relapse_scenario), function(rl) {
  det <- icer_of(set_field(base, "relapse", rl))                 # point estimate
  psa <- pmap_dfr(psa_params, function(...) {                    # PSA distribution
    r <- list(...); v <- icer_of(focal_iter(r, rl))
    tibble(averted = v[["averted"]], incr_cost = v[["incr_cost"]], icer = v[["icer"]])
  })
  tibble(
    relapse      = rl,
    scenario     = if_else(rl == 0, "No relapse (base case)",
                           sprintf("Relapse %.0f%%/cycle", 100 * rl)),
    averted      = det[["averted"]],   averted_lo = q025(psa$averted),   averted_hi = q975(psa$averted),
    incr_cost    = det[["incr_cost"]], ic_lo      = q025(psa$incr_cost), ic_hi      = q975(psa$incr_cost),
    icer         = det[["icer"]],      icer_lo    = q025(psa$icer),      icer_hi    = q975(psa$icer),
    p_ce_oc      = mean(oc_ng_r     * psa$averted - psa$incr_cost > 0),
    p_ce_1x      = mean(gdp_ng_r    * psa$averted - psa$incr_cost > 0),
    p_ce_3x      = mean(3 * gdp_ng_r * psa$averted - psa$incr_cost > 0))
})
write_tsv(relapse_tbl, file.path(CFG$dir_tab, "sensitivity_relapse.tsv"))
message("[07_sens] relapse scenario (point + 95% PSA CI):")
print(relapse_tbl %>% mutate(across(where(is.numeric), ~round(.x, 2))))

message("[07_sens] one-way (sorted by swing):")
print(oneway %>% transmute(param, icer_low = round(icer_low),
                           icer_high = round(icer_high),
                           inb_low = round(inb_low), inb_high = round(inb_high)))
message("[07_sens] wrote sensitivity tables + fig3_tornado + figS3_threshold_control")
