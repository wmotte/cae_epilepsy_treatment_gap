# ==============================================================================
# 11_additional_sensitivity.R
#
# Willem M. (Wim) Otte, w.m.otte@umcutrecht.nl
#
# Additional sensitivity analyses for inputs whose evidence is weak or contested.
# Results appear in appendix tables S12-S16.
#
# Every scenario re-uses the same 1000 paired draws as the base case, so each
# row is directly comparable to the headline. A scenario that changes one
# distribution re-maps the existing draw to the new distribution at the same
# quantile (remap_tri), keeping its pairing with every other parameter.
#
#   Mortality in untreated epilepsy
#     SMR-U-NGUGI   Ngugi 2014, Kilifi, SMR 6.5 (95% CI 5.0-8.3). Whole active
#                   convulsive epilepsy cohort, about half non-adherent.
#     SMR-U-LEVIRA  Levira 2017, weighted median SMR 2.6 of the higher-quality
#                   LMIC population-based studies. Fixed value.
#     SMR-U-EQ-T    Untreated equal to treated with continuing seizures (2.54,
#                   shared draw). The untreated SMR is probably not far from
#                   that of uncontrolled seizures. Mortality then falls
#                   only when seizures stop.
#     SMR-BEST      Source 95% CI bounds most favourable to treatment:
#                   untreated 10.0, treated 1.84, seizure-free 0.68.
#     SMR-WORST     Least favourable, with treated capped at untreated so that
#                   medication never raises mortality: untreated 2.0, treated
#                   2.0, seizure-free 1.29.
#     SMR-U-DECLINE SMRs fall with follow-up, so a 3-year SMR overstates
#                   excess mortality over 10 years. The untreated SMR falls
#                   linearly from its sampled value at entry to 2.6 (Levira
#                   median) at year 10, evaluated at each cycle midpoint.
#   Untreated SMR grid: the untreated SMR fixed at values from 2.0 to
#     10.0, all other inputs at their paired draws, and the break-even SMR at
#     which the mean ICER equals 0.5 x GDP and at which half of the draws are
#     cost-effective. Written to sensitivity_smr_grid.tsv / smr_breakeven.tsv.
#   UK treatment gap
#     UK-GAP-30     Wider gap, 30% -> 10%, read as a broad management gap.
#     UK-GAP-5      5% -> 1.25%, mid-range of European estimates.
#   Drug resistance (Sultana 2021)
#     DRE-9, DRE-19 Severe (drug-resistant) share at the bounds of the
#                   population-based 95% CI (9.2%, 19.0%) instead of 13.7%.
#     DRE-36        Clinic-based prevalence 36.3%.
#                   All three change the mixed population only.
#     UK-CTRL-85    UK cumulative seizure control 85% (80-90%) instead of 70%.
#   Costs
#     UK-INIT-500   UK case-finding and initiation US$500 (250-750) per person.
#
# A 10 000-draw re-run of the headline checks that 1000 draws are enough.
#
# THIS SCRIPT CHANGES NOTHING IN THE BASE CASE.
#
# Outputs (output/tables/):
#   sensitivity_additional.tsv scenario x country x cohort, with ICER 95% UIs
#   sensitivity_smr_grid.tsv   untreated SMR grid x country x cohort
#   smr_breakeven.tsv          break-even untreated SMR per country x cohort
#   mc_10000.tsv               headline probabilities at 1000 and 10 000 draws
# ==============================================================================

if (!exists("CFG")) {
  cand <- c("R/00_config.R", "00_config.R",
            file.path(Sys.getenv("CAE_ROOT", "."), "R/00_config.R"))
  source(cand[file.exists(cand)][1])
}
source(file.path(CFG$root, "R", "02_markov_model.R"))
source(file.path(CFG$root, "R", "lib_cea.R"))
suppressPackageStartupMessages(library(tidyverse))

mort   <- read_tsv(file.path(CFG$dir_derived, "mortality_schedule.tsv"),
                   show_col_types = FALSE)
params <- readRDS(file.path(CFG$dir_derived, "sim_traces.rds"))$params

tt <- (0:CFG$n_cycles) * CFG$cycle_years
r  <- CFG$discount

CTX <- set_names(CFG$countries$country) %>%
  map(function(ctry) {
    m <- mort[mort$country == ctry, ]
    m <- m[order(m$cycle), ]
    list(bg_rates = m$bg_rate, annuity_le = discounted_annuity(m$residual_le, r))
  })

trace_matrix <- function(p) {
  s  <- CFG$states
  x  <- setNames(c(p$initial_gap, 1 - p$initial_gap, 0, 0, 0), s)
  tr <- matrix(0, CFG$n_cycles + 1L, length(s), dimnames = list(NULL, s))
  tr[1, ] <- x
  for (k in seq_len(CFG$n_cycles)) {
    pk <- p; pk$bg_rate <- p$bg_rates[k]
    if (!is.null(p$smr_U_path)) pk$smr$Untreated <- p$smr_U_path[k]
    x <- as.numeric(x %*% build_transition_matrix(pk)); names(x) <- s
    tr[k + 1L, ] <- x
  }
  tr
}

# ---- Scenario definitions ----------------------------------------------------
# Each scenario is a function of (row, country) returning the resolved inputs
# it overrides. Anything it does not return keeps the base-case draw.
tri_U  <- c(CFG$smr_sens$Untreated[1], CFG$smr$Untreated, CFG$smr_sens$Untreated[2])
tri_HIC <- c(CFG$control_cum_sens$HIC[1], CFG$control_cum_10y$HIC,
             CFG$control_cum_sens$HIC[2])
tri_CF <- c(CFG$cost_case_finding_sens[1], CFG$cost_case_finding,
            CFG$cost_case_finding_sens[2])

SCENARIOS <- list(
  list(id = "BASE", group = "Base case",
       label = "Base case", f = function(row, ctry) list()),

  list(id = "SMR-U-NGUGI", group = "Untreated mortality",
       label = "Untreated SMR 6.5 (5.0-8.3), Ngugi 2014",
       f = function(row, ctry) list(smr_U = remap_tri(row$smr_U, tri_U, c(5.0, 6.5, 8.3)))),
  list(id = "SMR-U-LEVIRA", group = "Untreated mortality",
       label = "Untreated SMR 2.6, Levira 2017 weighted median",
       f = function(row, ctry) list(smr_U = 2.6)),
  list(id = "SMR-U-EQ-T", group = "Untreated mortality",
       label = "Untreated SMR equal to treated with continuing seizures",
       f = function(row, ctry) list(smr_U = row$smr_T)),
  list(id = "SMR-BEST", group = "Untreated mortality",
       label = "Best case: untreated 10.0, treated 1.84, seizure-free 0.68",
       f = function(row, ctry) list(smr_U = 10.0, smr_T = 1.84, smr_SF = 0.68)),
  list(id = "SMR-WORST", group = "Untreated mortality",
       label = "Worst case: untreated 2.0, treated 2.0, seizure-free 1.29",
       f = function(row, ctry) list(smr_U = 2.0, smr_T = 2.0, smr_SF = 1.29)),
  list(id = "SMR-U-DECLINE", group = "Untreated mortality",
       label = "Untreated SMR declining linearly to 2.6 by year 10",
       f = function(row, ctry) {
         mid <- (seq_len(CFG$n_cycles) - 0.5) * CFG$cycle_years
         list(smr_U_path = row$smr_U + (2.6 - row$smr_U) * mid / (CFG$n_cycles * CFG$cycle_years))
       }),

  list(id = "UK-GAP-30", group = "UK treatment gap", countries = "UK",
       label = "UK gap 30% to 10% (wider gap)",
       f = function(row, ctry) list(g0 = 0.30, g1 = 0.10)),
  list(id = "UK-GAP-5", group = "UK treatment gap", countries = "UK",
       label = "UK gap 5% to 1.25%",
       f = function(row, ctry) list(g0 = 0.05, g1 = 0.0125)),

  list(id = "DRE-9", group = "Drug resistance",
       label = "Severe (drug-resistant) share 9.2%, lower population bound",
       mixed_only = TRUE, sev_w = 0.092, f = function(row, ctry) list()),
  list(id = "DRE-19", group = "Drug resistance",
       label = "Severe (drug-resistant) share 19.0%, upper population bound",
       mixed_only = TRUE, sev_w = 0.190, f = function(row, ctry) list()),
  list(id = "DRE-36", group = "Drug resistance",
       label = "Severe (drug-resistant) share 36.3%, clinic-based",
       mixed_only = TRUE, sev_w = 0.363, f = function(row, ctry) list()),
  list(id = "UK-CTRL-85", group = "Drug resistance", countries = "UK",
       label = "UK cumulative seizure control 85% (80-90%)",
       f = function(row, ctry) list(control_cum = remap_tri(row$control_cum_HIC,
                                                            tri_HIC, c(0.80, 0.85, 0.90)))),

  list(id = "UK-INIT-500", group = "Costs", countries = "UK",
       label = "UK case-finding and initiation US$500 (250-750)",
       f = function(row, ctry) list(case_find = remap_tri(row$case_find_cost,
                                                          tri_CF, c(250, 500, 750))))
)

# ---- Value one paired draw ---------------------------------------------------
value_pair <- function(row, ctry, severity, ov) {
  setting <- CFG$countries$setting[CFG$countries$country == ctry]
  ctx <- CTX[[ctry]]
  g0 <- ov$g0 %||% CFG$gap_base[[setting]]
  g1 <- ov$g1 %||% CFG$gap_intervention[[setting]]
  smr_U  <- ov$smr_U  %||% row$smr_U
  smr_T  <- ov$smr_T  %||% row$smr_T
  smr_SF <- ov$smr_SF %||% row$smr_SF
  control_cum <- ov$control_cum %||%
    (if (setting == "LMIC") row$control_cum_LMIC else row$control_cum_HIC)
  case_find <- ov$case_find %||% row$case_find_cost
  p <- list(bg_rates = ctx$bg_rates,
            smr = list(Untreated = smr_U, Treated = smr_T,
                       Controlled = smr_SF, Remission = smr_SF),
            treat_uptake = CFG$background_treatment_uptake,
            spont_rem = row$spont_rem,
            control_prob = 1 - (1 - control_cum)^(1 / CFG$n_cycles),
            relapse = CFG$relapse_base,
            smr_U_path = ov$smr_U_path)
  dw <- dw_vec(row, severity)
  cost_yr <- CFG$cost_annual[[ctry]] * row$cost_mult
  arm <- function(gap) {
    p$initial_gap <- gap
    value_trace(tt, trace_matrix(p), dw, ctx$annuity_le, cost_yr,
                row$untreated_cost_frac, g0, gap, case_find, r)
  }
  sq <- arm(g0); itv <- arm(g1)
  c(da = sq[["daly"]] - itv[["daly"]], ic = itv[["cost"]] - sq[["cost"]])
}

summarise_draws <- function(da, ic, ctry) {
  oc  <- CFG$cet_oc[[ctry]]; gdp <- CFG$gdp_pc[[ctry]]
  inb <- oc * da - ic
  tibble(daly_averted = mean(da),
         da_lo = unname(quantile(da, .025)), da_hi = unname(quantile(da, .975)),
         incr_cost = mean(ic),
         ic_lo = unname(quantile(ic, .025)), ic_hi = unname(quantile(ic, .975)),
         icer = mean(ic) / mean(da),
         icer_lo = unname(quantile(ic / da, .025)),
         icer_hi = unname(quantile(ic / da, .975)),
         p_da_pos = mean(da > 0),
         inb_oc = mean(inb),
         inb_lo = unname(quantile(inb, .025)), inb_hi = unname(quantile(inb, .975)),
         p_ce_oc = mean(inb > 0),
         p_ce_1x = mean(gdp * da - ic > 0),
         p_ce_3x = mean(3 * gdp * da - ic > 0))
}

draws <- transpose(as.list(params))

# Cache per-draw severe and less-severe results by the overrides they depend
# on, so DRE-15/DRE-20 re-weight the base-case draws without re-running them.
run_country <- function(sc, ctry, draws) {
  n <- length(draws)
  res <- map(c(severe = "severe", less_severe = "less_severe"), function(sv) {
    m <- vapply(draws, function(row) value_pair(row, ctry, sv, sc$f(row, ctry)),
                c(da = 0, ic = 0))
    list(da = m["da", ], ic = m["ic", ])
  })
  ws <- sc$sev_w %||% CFG$severity_frac[["severe"]]
  mixed <- list(da = ws * res$severe$da + (1 - ws) * res$less_severe$da,
                ic = ws * res$severe$ic + (1 - ws) * res$less_severe$ic)
  out <- bind_rows(
    summarise_draws(res$severe$da, res$severe$ic, ctry) %>% mutate(severity = "severe"),
    summarise_draws(res$less_severe$da, res$less_severe$ic, ctry) %>%
      mutate(severity = "less_severe"),
    summarise_draws(mixed$da, mixed$ic, ctry) %>% mutate(severity = "mixed"))
  if (isTRUE(sc$mixed_only)) out <- filter(out, severity == "mixed")
  out %>% mutate(country = ctry, .before = 1)
}

results <- map_dfr(SCENARIOS, function(sc) {
  message("[11_additional] ", sc$id)
  ctries <- sc$countries %||% CFG$countries$country
  map_dfr(ctries, ~ run_country(sc, .x, draws)) %>%
    mutate(scenario_id = sc$id, group = sc$group, label = sc$label, .before = 1)
})
write_tsv(results, file.path(CFG$dir_tab, "sensitivity_additional.tsv"))

# ---- Untreated SMR grid and break-even value ---------------------------------
# The untreated SMR is fixed at each grid value in every draw, so the grid shows
# how the verdict moves with this one input while every other input keeps its
# paired draw. The break-even SMR is found by root-finding on the same draws.
SMR_GRID <- c(2.0, 2.54, 3.0, 3.5, 4.0, 4.5, 5.0, 6.3, 8.0, 10.0)
grid_sc <- function(v) list(id = "SMR-GRID", f = function(row, ctry) list(smr_U = v))
smr_grid <- map_dfr(SMR_GRID, function(v) {
  message("[11_additional] SMR grid ", v)
  map_dfr(CFG$countries$country, ~ run_country(grid_sc(v), .x, draws)) %>%
    mutate(smr_U = v, .before = 1)
})
write_tsv(smr_grid, file.path(CFG$dir_tab, "sensitivity_smr_grid.tsv"))

# Per-draw values at a fixed untreated SMR, returned for both cohorts and the
# mixed population, so the root-finder re-uses one evaluation per SMR value.
eval_smr <- function(v, ctry) {
  res <- map(c(severe = "severe", less_severe = "less_severe"), function(sv) {
    m <- vapply(draws, function(row) value_pair(row, ctry, sv, list(smr_U = v)),
                c(da = 0, ic = 0))
    list(da = m["da", ], ic = m["ic", ])
  })
  ws <- CFG$severity_frac[["severe"]]
  res$mixed <- list(da = ws * res$severe$da + (1 - ws) * res$less_severe$da,
                    ic = ws * res$severe$ic + (1 - ws) * res$less_severe$ic)
  res
}
breakeven <- map_dfr(CFG$countries$country, function(ctry) {
  oc <- CFG$cet_oc[[ctry]]
  cache <- new.env()
  get <- function(v) {
    key <- sprintf("%.6f", v)
    if (is.null(cache[[key]])) cache[[key]] <- eval_smr(v, ctry)
    cache[[key]]
  }
  map_dfr(c("severe", "less_severe", "mixed"), function(sv) {
    f_icer <- function(v) { x <- get(v)[[sv]]; oc * mean(x$da) - mean(x$ic) }
    f_pce  <- function(v) { x <- get(v)[[sv]]; mean(oc * x$da - x$ic > 0) - 0.5 }
    root <- function(f) {
      lo <- f(1.0); hi <- f(12)
      if (lo >= 0) return(-Inf)            # cost-effective even at SMR 1.0
      if (hi < 0) return(Inf)              # not cost-effective even at SMR 12
      uniroot(f, c(1.0, 12), tol = 1e-3)$root
    }
    tibble(country = ctry, severity = sv,
           smr_icer_eq_threshold = root(f_icer),
           smr_pce_50 = root(f_pce))
  })
})
write_tsv(breakeven, file.path(CFG$dir_tab, "smr_breakeven.tsv"))

# ---- Monte Carlo stability: 10 000 draws for the headline --------------------
# Fresh draws from the same sampler with a distinct seed, so the 1000-draw
# headline is not a subset of this run.
set.seed(CFG$seed + 1L)
big <- transpose(as.list(sample_psa_params(10000L)))
mc <- map_dfr(CFG$countries$country, function(ctry) {
  run_country(SCENARIOS[[1]], ctry, big) %>% mutate(n_draws = 10000L)
}) %>%
  bind_rows(results %>% filter(scenario_id == "BASE") %>%
              select(-scenario_id, -group, -label) %>% mutate(n_draws = CFG$n_psa)) %>%
  mutate(mc_se = sqrt(p_ce_oc * (1 - p_ce_oc) / n_draws),
         mc_lo = pmax(0, p_ce_oc - 1.96 * mc_se),
         mc_hi = pmin(1, p_ce_oc + 1.96 * mc_se)) %>%
  select(country, severity, n_draws, daly_averted, incr_cost, icer, p_ce_oc,
         mc_se, mc_lo, mc_hi) %>%
  arrange(country, severity, n_draws)
write_tsv(mc, file.path(CFG$dir_tab, "mc_10000.tsv"))

# ---- Console summary ---------------------------------------------------------
message("\n[11_additional] Probability cost-effective at 0.5x GDP (severe | mixed)")
print(results %>%
  filter(severity %in% c("severe", "mixed")) %>%
  transmute(scenario_id, country, severity, da = round(daly_averted, 2),
            icer = round(icer), p = sprintf("%.1f%%", 100 * p_ce_oc)) %>%
  pivot_wider(names_from = severity, values_from = c(da, icer, p)), n = 60)
message("\n[11_additional] Monte Carlo stability")
print(mc %>% mutate(across(c(p_ce_oc, mc_lo, mc_hi), ~ sprintf("%.1f%%", 100 * .x))))
message("\n[11_additional] Break-even untreated SMR")
print(breakeven)
message("[11_additional] wrote sensitivity_additional.tsv, sensitivity_smr_grid.tsv, ",
        "smr_breakeven.tsv, mc_10000.tsv")
