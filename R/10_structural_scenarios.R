# ==============================================================================
# 10_structural_scenarios.R
#
# Willem M. (Wim) Otte, w.m.otte@umcutrecht.nl
#
# Structural (not parameter) uncertainty: what happens to the headline when the
# model's contested source-to-state mappings are replaced by defensible
# alternatives?
#
# The probabilistic analysis in 03/04 samples *around* the base-case mappings.
# It cannot tell a reader what the answer would be if a mapping is wrong. Four
# independent reviews converged on the same five mappings, and each was checked
# against the primary source:
#
#   DW-RETAIN  Salomon 2015 defines 0.552 and 0.263 by SEIZURE FREQUENCY
#              (>=1/month vs 1-11/year), not by treatment status. The base case
#              drops a severe person from 0.552 to 0.263 the moment medication
#              starts, before any modelled clinical response. This scenario keeps
#              the entry weight until an explicit transition to Controlled.
#
#   SMR-SF1    Mohanraj 2006 (PMID 16713919) reports "no increase in risk
#              observed in patients who were seizure free". Since the 2026-09
#              revision the base case already uses its seizure-free SMR, 0.95
#              (0.68-1.29), for Controlled and Remission. This scenario sets
#              both to exactly 1.0, removing any benefit below population rates.
#
#   SMR-NOMORT Carpio's 6.3 is from a newly diagnosed cohort with no untreated
#              subgroup. Its contrast with Mohanraj's 2.54 is therefore not an
#              identified causal effect of starting ASM. This scenario removes
#              any direct mortality benefit of medication: Treated keeps the
#              ratio assigned to Untreated, and mortality improves only on
#              reaching seizure freedom.
#
#   ECU-UMIC   Begley 2022 (PMID 35195894) reports 2019 USD. US$354.20 is its
#              LOWER-middle-income group mean; Ecuador is upper-middle income,
#              whose group mean in the same table is US$2048.20.
#
#   YLL-LT     The base case charges each death the discounted annuity of the
#              life-table expectation e_x. The expected discounted life-years is
#              the integral of exp(-ru)S(u). The annuity is concave, so the
#              convention is an upper bound. This scenario integrates the
#              survival curve instead.
#
# CONSERVATIVE stacks DW-RETAIN + SMR-SF1 + SMR-NOMORT.
#
# The script also reports the MIXED economic population (13.7% severe / 86.3% less
# severe, weighted WITHIN each paired draw). The programme treats a whole
# prevalent cohort, so the mixed cohort -- not the severe subgroup -- is the
# population a health ministry would decide over.
#
# THIS SCRIPT CHANGES NOTHING IN THE BASE CASE. It re-uses the same 1000 paired
# draws, so every scenario is directly comparable to the headline.
#
# Outputs (output/tables/):
#   structural_scenarios.tsv        : scenario x country x cohort decision table
#   structural_scenarios_yll.tsv    : annuity vs life-table YLL comparison
#   daly_decomposition.tsv          : DALYs averted split into YLD and YLL
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

tt      <- (0:CFG$n_cycles) * CFG$cycle_years
r       <- CFG$discount
t_mid   <- tt[-1] - CFG$cycle_years / 2

# ---- Life tables, for the survival-integral YLL --------------------------
life_tables <- set_names(CFG$countries$country) %>%
  map(~ read_tsv(file.path(CFG$dir_derived, paste0("life_table__", .x, ".tsv")),
                 show_col_types = FALSE))

# Discounted life-years remaining at each cycle midpoint age, per country. The
# base case uses discounted_annuity(e_x); YLL-LT uses this instead.
age_mids <- CFG$index_age + t_mid

# ---- Per-country context, resolved once --------------------------------------
# The inner loop runs ~84 000 times across scenarios; a dplyr::filter in there
# dominates the runtime. Everything that depends only on the country is hoisted.
CTX <- set_names(CFG$countries$country) %>%
  imap(function(ctry, .) {
    m <- mort[mort$country == ctry, ]
    m <- m[order(m$cycle), ]
    list(
      bg_rates    = m$bg_rate,
      residual_le = m$residual_le,
      annuity_le  = discounted_annuity(m$residual_le, r),
      lt_le       = map_dbl(age_mids,
                            ~ discounted_le_lifetable(life_tables[[ctry]], .x, r)))
  })

# ---- Lean trace: numeric matrix, no tibble round-trip ------------------------
trace_matrix <- function(p) {
  s  <- CFG$states
  x  <- setNames(c(p$initial_gap, 1 - p$initial_gap, 0, 0, 0), s)
  tr <- matrix(0, CFG$n_cycles + 1L, length(s), dimnames = list(NULL, s))
  tr[1, ] <- x
  for (k in seq_len(CFG$n_cycles)) {
    pk <- p; pk$bg_rate <- p$bg_rates[k]
    x <- as.numeric(x %*% build_transition_matrix(pk)); names(x) <- s
    tr[k + 1L, ] <- x
  }
  tr
}

# ---- Value one arm ----------------------------------------------------------
# `sc` is the scenario definition; `row` one PSA draw (a plain list).
value_arm <- function(row, ctry, setting, severity, gap, sc) {
  ctx <- CTX[[ctry]]

  smr_U <- row$smr_U
  smr_T <- if (sc$no_direct_mortality_benefit) row$smr_U else row$smr_T
  smr_C <- if (sc$seizure_free_smr_one) 1.0 else row$smr_SF
  smr_R <- if (sc$seizure_free_smr_one) 1.0 else row$smr_SF

  control_cum <- if (setting == "LMIC") row$control_cum_LMIC else row$control_cum_HIC

  tr <- trace_matrix(list(
    bg_rates     = ctx$bg_rates,
    initial_gap  = gap,
    smr          = list(Untreated = smr_U, Treated = smr_T,
                        Controlled = smr_C, Remission = smr_R),
    treat_uptake = CFG$background_treatment_uptake,
    spont_rem    = row$spont_rem,
    control_prob = 1 - (1 - control_cum)^(1 / CFG$n_cycles),
    relapse      = CFG$relapse_base
  ))

  dw <- dw_vec(row, severity, treated_keeps_entry_dw = sc$treated_keeps_entry_dw)
  le_disc <- if (sc$yll_life_table) ctx$lt_le else ctx$annuity_le
  cost_yr <- (if (sc$ecuador_umic_cost && ctry == "Ecuador")
                CFG$cost_annual_ecu_umic else CFG$cost_annual[[ctry]]) * row$cost_mult
  value_trace(
    tt, tr, dw, le_disc, cost_yr, row$untreated_cost_frac,
    CFG$gap_base[[setting]], gap, row$case_find_cost, r
  )[c("yld", "yll", "daly", "cost")]
}

# ---- Scenario definitions ---------------------------------------------------
sc_default <- list(treated_keeps_entry_dw = FALSE, seizure_free_smr_one = FALSE,
                   no_direct_mortality_benefit = FALSE, ecuador_umic_cost = FALSE,
                   yll_life_table = FALSE)
mk <- function(...) modifyList(sc_default, list(...))

SCENARIOS <- list(
  `Base case`                              = mk(),
  `DW-RETAIN: treated keeps entry weight`  = mk(treated_keeps_entry_dw = TRUE),
  `SMR-SF1: seizure-free SMR = 1`          = mk(seizure_free_smr_one = TRUE),
  `SMR-NOMORT: no direct mortality benefit`= mk(no_direct_mortality_benefit = TRUE),
  `ECU-UMIC: Ecuador upper-middle cost`    = mk(ecuador_umic_cost = TRUE),
  `YLL-LT: life-table survival YLL`        = mk(yll_life_table = TRUE),
  `CONSERVATIVE: DW-RETAIN + SMR-SF1 + SMR-NOMORT` =
      mk(treated_keeps_entry_dw = TRUE, seizure_free_smr_one = TRUE,
         no_direct_mortality_benefit = TRUE)
)

# ---- Run every scenario over the shared draws --------------------------------
sev_w <- CFG$severity_frac
# Plain list-of-lists: the inner loop touches these ~84 000 times, and tibble
# column extraction there is the single hottest cost in the script.
draws <- transpose(as.list(params))
n_it  <- length(draws)

summarise_draws <- function(da, ic, yld, yll, oc) {
  inb <- oc * da - ic
  tibble(daly_averted = mean(da), yld_averted = mean(yld), yll_averted = mean(yll),
         incr_cost = mean(ic), icer = mean(ic) / mean(da),
         inb_oc = mean(inb),
         inb_lo = unname(quantile(inb, .025)), inb_hi = unname(quantile(inb, .975)),
         p_ce_oc = mean(inb > 0), p_benefit = mean(da > 0))
}

run_scenario <- function(sc_name, sc) {
  message("[10_structural] ", sc_name)
  out <- map_dfr(seq_len(nrow(CFG$countries)), function(i) {
    ctry    <- CFG$countries$country[i]
    setting <- CFG$countries$setting[i]
    g0 <- CFG$gap_base[[setting]]; g1 <- CFG$gap_intervention[[setting]]
    oc <- CFG$cet_oc[[ctry]]

    sevs <- c("severe", "less_severe")
    acc <- list()
    for (sv in sevs)
      acc[[sv]] <- list(da = numeric(n_it), ic = numeric(n_it),
                        yld = numeric(n_it), yll = numeric(n_it))

    for (j in seq_len(n_it)) {
      row <- draws[[j]]
      for (sv in sevs) {
        sq  <- value_arm(row, ctry, setting, sv, g0, sc)
        itv <- value_arm(row, ctry, setting, sv, g1, sc)
        acc[[sv]]$da[j]  <- sq[["daly"]] - itv[["daly"]]
        acc[[sv]]$yld[j] <- sq[["yld"]]  - itv[["yld"]]
        acc[[sv]]$yll[j] <- sq[["yll"]]  - itv[["yll"]]
        acc[[sv]]$ic[j]  <- itv[["cost"]] - sq[["cost"]]
      }
    }

    # Mixed economic population: weighted WITHIN each paired draw, so combined
    # uncertainty comes from combined draws rather than averaged quantiles.
    ws <- sev_w[["severe"]]; wl <- sev_w[["less_severe"]]
    mixed <- list(
      da  = ws * acc$severe$da  + wl * acc$less_severe$da,
      ic  = ws * acc$severe$ic  + wl * acc$less_severe$ic,
      yld = ws * acc$severe$yld + wl * acc$less_severe$yld,
      yll = ws * acc$severe$yll + wl * acc$less_severe$yll)

    bind_rows(
      summarise_draws(acc$severe$da, acc$severe$ic,
                      acc$severe$yld, acc$severe$yll, oc) %>%
        mutate(severity = "severe", .before = 1),
      summarise_draws(acc$less_severe$da, acc$less_severe$ic,
                      acc$less_severe$yld, acc$less_severe$yll, oc) %>%
        mutate(severity = "less_severe", .before = 1),
      summarise_draws(mixed$da, mixed$ic, mixed$yld, mixed$yll, oc) %>%
        mutate(severity = "mixed", .before = 1)
    ) %>% mutate(country = ctry, .before = 1)
  })
  out %>% mutate(scenario = sc_name, .before = 1)
}

results <- imap_dfr(SCENARIOS, ~ run_scenario(.y, .x))
write_tsv(results, file.path(CFG$dir_tab, "structural_scenarios.tsv"))

# ---- DALY decomposition (base case): how much of the gain is mortality? ------
decomp <- results %>%
  filter(scenario == "Base case") %>%
  transmute(country, severity, daly_averted, yld_averted, yll_averted,
            yll_share = yll_averted / daly_averted)
write_tsv(decomp, file.path(CFG$dir_tab, "daly_decomposition.tsv"))

# ---- YLL valuation comparison ------------------------------------------------
yll_cmp <- tibble(
  country = rep(CFG$countries$country, each = length(age_mids)),
  age     = rep(age_mids, times = nrow(CFG$countries))) %>%
  rowwise() %>%
  mutate(
    residual_le = mort$residual_le[mort$country == country &
                                   mort$age_mid == age][1],
    annuity_of_ex     = discounted_annuity(residual_le, r),
    lifetable_integral = discounted_le_lifetable(life_tables[[country]], age, r)) %>%
  ungroup() %>%
  mutate(ratio = annuity_of_ex / lifetable_integral)
write_tsv(yll_cmp, file.path(CFG$dir_tab, "structural_scenarios_yll.tsv"))

# ---- Console summary ---------------------------------------------------------
message("\n[10_structural] Decision table: MIXED prevalent cohort (13.7% severe / 86.3% less severe)")
print(results %>%
  filter(severity == "mixed") %>%
  transmute(scenario, country,
            daly_averted = round(daly_averted, 2),
            incr_cost = round(incr_cost),
            icer = round(icer),
            inb_oc = round(inb_oc),
            p_ce_oc = sprintf("%.1f%%", 100 * p_ce_oc)) %>%
  arrange(country, scenario), n = 40)

message("\n[10_structural] Severe cohort, probability of cost-effectiveness at the 0.5x GDP proxy")
print(results %>%
  filter(severity == "severe") %>%
  select(scenario, country, p_ce_oc) %>%
  mutate(p_ce_oc = sprintf("%.1f%%", 100 * p_ce_oc)) %>%
  pivot_wider(names_from = country, values_from = p_ce_oc), n = 20)

message("\n[10_structural] Base-case DALY decomposition (share of averted DALYs from mortality)")
print(decomp %>% mutate(yll_share = sprintf("%.1f%%", 100 * yll_share)), n = 20)

message("\n[10_structural] YLL valuation: annuity of e_x vs discounted life-table survival")
print(yll_cmp %>%
  group_by(country) %>%
  summarise(mean_annuity = round(mean(annuity_of_ex), 2),
            mean_lifetable = round(mean(lifetable_integral), 2),
            mean_ratio = round(mean(ratio), 3)))

message("[10_structural] wrote structural_scenarios.tsv, daly_decomposition.tsv, structural_scenarios_yll.tsv")
