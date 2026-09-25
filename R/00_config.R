# ==============================================================================
# 00_config.R
# Central configuration for the epilepsy treatment-gap cost-effectiveness model.
#
# Willem M. (Wim) Otte, w.m.otte@umcutrecht.nl
#
# Every parameter, source path, and assumption used by the model lives here so
# that the rest of the pipeline (01-09) contains no hard-coded magic numbers.
# The generated parameter table with sources is output/tables/table_parameters.tsv.
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

# ---- Paths -------------------------------------------------------------------
# Resolve project root by searching upward for a marker dir, so the pipeline is
# standalone and path-independent (works from any wd, source() or Rscript).
CFG <- list()
.find_root <- function() {
  # 1. honour explicit override
  env <- Sys.getenv("CAE_ROOT", "")
  if (nzchar(env) && dir.exists(file.path(env, "R"))) return(normalizePath(env))
  # 2. script path from Rscript --file=
  args <- commandArgs(FALSE)
  fa <- grep("^--file=", args, value = TRUE)
  start <- if (length(fa)) dirname(sub("^--file=", "", fa[1])) else getwd()
  d <- normalizePath(start, mustWork = FALSE)
  for (i in 1:6) {
    if (dir.exists(file.path(d, "R")) &&
        file.exists(file.path(d, "R", "00_config.R"))) return(d)
    parent <- dirname(d)
    if (parent == d) break
    d <- parent
  }
  normalizePath(getwd(), mustWork = FALSE)
}
CFG$root <- .find_root()
CFG$dir_raw      <- file.path(CFG$root, "data", "raw")
CFG$dir_derived  <- file.path(CFG$root, "data", "derived")
CFG$dir_fig      <- file.path(CFG$root, "output", "figures")
CFG$dir_diag     <- file.path(CFG$root, "output", "diagnostics")
CFG$dir_tab      <- file.path(CFG$root, "output", "tables")
CFG$wpp_xlsx     <- file.path(CFG$dir_raw,
  "WPP2024_MORT_F07_1_ABRIDGED_LIFE_TABLE_BOTH_SEXES.xlsx")

for (d in c(CFG$dir_derived, CFG$dir_fig, CFG$dir_diag, CFG$dir_tab)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ---- Reproducibility ---------------------------------------------------------
CFG$seed   <- 20260520L
CFG$n_psa  <- 1000L         # Monte Carlo PSA iterations. Tail-probability claims
                            # need a stable MC estimate (MC SE on a 0.9
                            # proportion is ~3% at n=100, ~1% at n=1000)

# ---- Time horizon ------------------------------------------------------------
CFG$cycle_years <- 2L       # cycle length
CFG$n_cycles    <- 5L       # 5 cycles
CFG$horizon_y   <- CFG$cycle_years * CFG$n_cycles   # 10-year horizon
CFG$discount    <- 0.035    # annual discount rate (3.5%)
CFG$index_age   <- 30       # explicit starting age of the representative cohort

# ---- Health states (fixed order) ---------------------------------------------
# 1 Untreated  : symptomatic, no/insufficient ASM (the treatment gap)
# 2 Treated    : on ASM, still seizing
# 3 Controlled : on ASM, seizure-free
# 4 Remission  : spontaneous, off ASM, seizure-free
# 5 Death      : absorbing
CFG$states <- c("Untreated", "Treated", "Controlled", "Remission", "Death")

# ---- Disability weights (GBD 2013 Salomon / GBD mapping) -------
# Severity is a cohort attribute: the "severe" cohort enters Untreated with the
# severe-epilepsy DW; the "less-severe" cohort with the moderate DW.
CFG$dw <- list(
  Untreated_severe      = 0.552,  # GBD2013 severe epilepsy (UI 0.375-0.710)
  Untreated_less_severe = 0.263,  # GBD less-severe epilepsy
  Treated               = 0.263,  # on ASM, still seizing (shares the 0.263 band)
  Controlled            = 0.049,  # on ASM, seizure-free
  Remission             = 0.049,  # off ASM, seizure-free (shares Controlled draw)
  Death                 = 0.0
)
# Severity split of the epilepsy population. Drug resistance is used as a PROXY
# for the higher-seizure-frequency / higher-disability cohort; the two constructs
# are related but not identical (see manuscript Limitations). The base case is the population/community-based
# DRE prevalence of Sultana et al. 2021, 13.7% (95% CI 9.2-19.0). The clinic-based
# value of 36.3% (30.4-42.4) and the population CI bounds are scenarios in
# 11_additional_sensitivity.R. Used to
# weight severe and less-severe cohorts into a combined population-level burden.
CFG$severity_frac <- c(severe = 0.137, less_severe = 0.863)

# Plausible limits used by the triangular PSA distributions. Read by the sampler
# (03) and by the boundary checks (06), so a bound cannot drift between the two.
# The 0.263 band is 0.173-0.367. The 2026-05-20 standalone rebuild introduced
# 0.179-0.359 as a fresh literal with no source note, and the later manuscript
# variants then inherited it from the code. The author's accepted manuscript and
# the variant-B builders both carry 0.173-0.367, sourced to Salomon 2015, and the
# two flanking bands are identical in every document.
CFG$dw_ci <- tribble(
  ~param,                  ~lower, ~upper,
  "Untreated_severe",       0.375,  0.710,
  "Untreated_less_severe",  0.173,  0.367,
  "Treated",                0.173,  0.367,
  "Controlled",             0.031,  0.072,
  "Remission",              0.031,  0.072
)
# Bounds lookup: dw_bounds("Controlled") -> c(lower, upper).
dw_bounds <- function(param) {
  row <- CFG$dw_ci[CFG$dw_ci$param == param, ]
  stopifnot(nrow(row) == 1L)
  c(row$lower, row$upper)
}

# ---- Standardized mortality ratios ----------------
# The treated and seizure-free SMRs come from ONE source,
# Mohanraj et al. 2006 (PMID 16713919), stratified by response in its newly
# diagnosed cohort:
#   Treated    2.54  (95% CI 1.84-3.44)  on ASM, did not respond (42 deaths,
#                    16.5 expected, n=318).
#   Controlled 0.95  (95% CI 0.68-1.29)  entered remission on ASM (41 deaths,
#                    43.1 expected, n=462).
#   Remission  0.95  Seizure-free off medication. No separate estimate exists,
#                    so it shares the seizure-free value and its PSA draw.
# Transporting these from a newly diagnosed Scottish cohort to a prevalent
# cohort in three settings is an assumption. The paper's chronic cohort (SMR
# 2.04) was not stratified by seizure control.
#
# Untreated  6.3   Carpio et al. 2005. The Ecuador component followed a newly
#                  diagnosed cohort and reported an overall SMR of 6.3 (95% CI
#                  2.0-10.0). It did not define an untreated subgroup. Assignment
#                  of that estimate to Untreated is a structural model mapping.
#                  Alternatives (Ngugi 2014 6.5, Levira 2017 median 2.6, equal
#                  to Treated, best and worst case) are run in
#                  11_additional_sensitivity.R.
CFG$smr <- list(
  Untreated  = 6.3,
  Treated    = 2.54,
  Controlled = 0.95,
  Remission  = 0.95,
  Death      = NA_real_
)
# SMR limits for the triangular PSA, each the source's published 95% CI.
# Controlled and Remission are one construct (seizure-free) and share one draw.
CFG$smr_sens <- list(
  Untreated   = c(2.0,  10.0),
  Treated     = c(1.84,  3.44),
  SeizureFree = c(0.68,  1.29)
)

# ---- Transition parameters ---------------------------------------------------
# Spontaneous remission from Untreated.
# Nicoletti 2009 reports ~30-44% *cumulative* 5-year remission over a ~10-year
# follow-up in untreated rural Bolivia. Applying 0.30 per 2-year cycle would
# compound to ~83% spontaneous remission over the horizon, which is implausible.
# We therefore treat ~30% as a cumulative 10-year figure and convert to a per-cycle probability:
#   annual rate r = -log(1 - 0.30)/10 = 0.0357; p(2y) = 1 - exp(-2r) = 0.069.
CFG$spont_remission_prob <- 0.069  # Untreated -> Remission share (per 2y cycle)
# Sensitivity range from the 20%-44% cumulative-10y span -> per-cycle 0.044-0.110.
CFG$spont_remission_sens <- c(0.044, 0.110)

# Published seizure-control proportions are cumulative outcomes, not probabilities
# that can be reapplied unchanged every 2 years. We therefore specify 10-year
# cumulative control and convert it to a conditional per-cycle probability.
CFG$control_cum_10y <- list(LMIC = 0.50, HIC = 0.70)
CFG$control_cum_sens <- list(LMIC = c(0.35, 0.70), HIC = c(0.55, 0.80))
CFG$control_prob <- map(CFG$control_cum_10y,
  ~ 1 - (1 - .x)^(1 / CFG$n_cycles))

# Optional relapse (secondary treatment gap). Base case = 0 (no relapse). A
# scenario analysis (07) sets a per-cycle probability of moving
# Treated/Controlled/Remission back to Untreated, to bound the optimism of the
# no-relapse structure.
CFG$relapse_base     <- 0.00
CFG$relapse_scenario <- 0.15   # per-2y-cycle relapse to Untreated (scenario)

# ---- Countries ---------------------------------------------------------------
# setting drives control_prob; iso3 maps to WPP life table.
CFG$countries <- tribble(
  ~country,  ~iso3, ~setting, ~prevalence_per_1000, ~prev_lo, ~prev_hi,
  "Nigeria", "NGA", "LMIC",   9.8,                  8.6,      11.1,   # Watila 2021
  "Ecuador", "ECU", "LMIC",   7.5,                  6.0,      9.0,    # Placencia/GBD band
  "UK",      "GBR", "HIC",    9.37,                 8.85,     12.08   # Wigglesworth 2023 (CPRD); pt 9.37 (95%CI 9.34-9.40), band = nation range England 8.85 - NI 12.08
)

# Total national populations (for population-level burden scaling).
CFG$population <- c(
  Nigeria = 240e6,   # ~2024
  Ecuador = 18e6,    # ~2024
  UK      = 69e6     # ~2024
)
CFG$population_ssa <- 1.273663761e9   # Sub-Saharan Africa (burden framing)

# WPP reference year for background mortality.
CFG$wpp_year <- 2023L

# ---- Costs (direct annual cost per person on ASM, 2019 USD) ---
# Applied to Treated + Controlled (ongoing ASM + care). Untreated incurs a small
# fraction (informal/crisis care); Remission and Death = 0.
#
# PRICE YEAR. Begley et al. 2022 (PMID 35195894) state that country costs were
# "extracted and adjusted to generate an average cost per person in 2019 US
# dollars", using the country GDP deflator and 2019 purchasing-power-parity
# exchange rates. Earlier text in this repository asserted that no common price
# year could be reconstructed from the source. That was wrong. These are 2019 USD.
#
# ECUADOR. US$354.20 is Begley's *lower-middle-income* (World Bank category 2)
# group mean, i.e. the average of India and Nigeria. Ecuador is upper-middle
# income, and the same source table reports a category-3 mean of US$2048.20.
# The base case therefore understates Ecuadorian cost by roughly 5.8x, which
# biases the Ecuadorian ICER favourably; the +/-25% PSA multiplier cannot span
# that gap. Structural scenario ECU-UMIC in 10_structural_scenarios.R re-runs
# Ecuador at US$2048.20.
CFG$cost_annual <- c(
  Nigeria = 654.59,   # Begley 2022, country estimate (2019 USD)
  Ecuador = 354.20,   # Begley 2022, lower-middle-income group mean (2019 USD)
  UK      = 2320.76   # Begley 2022, country estimate (2019 USD)
)
# Upper-middle-income group mean from the same source table, used by the Ecuador
# cost structural scenario.
CFG$cost_annual_ecu_umic <- 2048.20
CFG$cost_untreated_frac <- 0.10   # documented assumption: untreated formal-care cost
# This fraction is the only input whose plausible range can change the SIGN of
# the incremental cost, so it belongs in the one-way analysis. The bounds are
# deliberately generous: 0 is "no formal care reaches the untreated at all", and
# 1 is the extreme where crisis, injury, and emergency care for untreated
# epilepsy costs as much as full ASM care. Bridging stays cost-incurring across
# the whole range; incremental cost only reaches zero near 1.43, which no cost
# source supports.
CFG$cost_untreated_frac_sens <- c(0, 1.00)
# The probabilistic analysis samples a narrower, plausible range. The one-way
# bound above is a deliberate stress test up to parity with full ASM care, which
# no cost source supports as a *likely* value, so carrying it into the PSA would
# let an unsupported extreme dominate the cost distribution.
CFG$cost_untreated_frac_psa <- c(0, 0.30)
# Annual care costs are varied by a shared multiplier in the PSA.
CFG$cost_mult_sens <- c(0.75, 1.25)

# One-off case-finding / treatment-initiation cost for each additional person
# covered at model entry.
CFG$cost_case_finding      <- 50      # USD per newly treated person (assumption)
CFG$cost_case_finding_sens <- c(0, 150)

# ---- Cost-effectiveness thresholds (GDP per capita, USD) ---------------------
# 1x and 3x GDP per capita (WHO-CHOICE convention; see table_parameters.tsv).
# NB: the 1x/3x-GDP thresholds are widely criticised as too permissive
# (Bertram et al., 2016) and have been effectively retired by WHO. We retain
# them for comparability but ALSO report an opportunity-cost-based threshold.
CFG$gdp_pc <- c(
  Nigeria = 2000,
  Ecuador = 6500,
  UK      = 49000
)
# Opportunity-cost (health-system) threshold, as an explicit proxy of 0.5x GDP
# per capita. Empirical work (Woods et al., 2016; Ochalek et al., 2018) places
# supply-side thresholds for low-income settings well below 1x GDP, commonly
# ~0.5x or lower. Used as a stricter, more defensible benchmark alongside WHO.
CFG$cet_oc_frac <- 0.5
CFG$cet_oc <- CFG$cet_oc_frac * CFG$gdp_pc

# ---- Treatment-gap scenario sweep -------------------------------------------
# Treatment gap is a stock: the fraction of the prevalent cohort untreated at
# model entry. Bridging reallocates the additional covered share from Untreated
# to Treated at time zero. It is not reapplied as a transition every cycle.
CFG$gap_sweep <- list(
  LMIC = seq(0.20, 0.80, by = 0.10),    # 20-80% residual gap (7 levels)
  HIC  = seq(0.025, 0.10, by = 0.025)   # 2.5-10% residual gap (4 levels)
)
# Base-case representative gaps for headline numbers.
# The UK status quo is 10%, the upper bound for high-income countries in Meyer 2010, reduced by
# the same 75% relative reduction as the LMIC scenario (80% -> 20%). A wider
# 30% -> 10% and a 5% -> 1.25% scenario are run in 11_additional_sensitivity.R.
CFG$gap_base <- list(LMIC = 0.80, HIC = 0.10)   # status-quo gap
CFG$gap_intervention <- list(LMIC = 0.20, HIC = 0.025)  # reduced gap

# No background treatment initiation is imposed after model entry because a
# comparable country-specific rate was not available. This makes the estimand an
# immediate, sustained coverage change and is tested as a structural limitation.
CFG$background_treatment_uptake <- 0

# ---- Helpers -----------------------------------------------------------------
# Convert an annual rate to a probability over t years.
rate_to_prob <- function(rate, t = CFG$cycle_years) 1 - exp(-rate * t)
# Convert a per-cycle probability to an annual rate.
prob_to_rate <- function(prob, t = CFG$cycle_years) -log(1 - prob) / t

set.seed(CFG$seed)
message(sprintf("[00_config] root=%s | %d states | %d-y horizon | %d PSA",
                CFG$root, length(CFG$states), CFG$horizon_y, CFG$n_psa))
