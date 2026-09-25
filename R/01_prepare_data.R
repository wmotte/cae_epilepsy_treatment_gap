# ==============================================================================
# 01_prepare_data.R
#
# Willem M. (Wim) Otte, w.m.otte@umcutrecht.nl
#
# Derive country-specific background mortality from the UN WPP2024 abridged
# life table (both sexes) for Nigeria (NGA), Ecuador (ECU) and the UK (GBR).
#
# Outputs (data/derived/):
#   life_table__<country>.tsv  : abridged life table (age, n, m(x,n), q(x,n))
#   mortality_schedule.tsv     : age-varying background death rate and residual
#                                life expectancy for each model cycle
#   prevalence.tsv             : prevalence + population per country (from config)
#
# The model follows an explicit representative cohort starting at age 30. Rates
# and residual life expectancy are interpolated from the abridged life table at
# each cycle midpoint. This avoids applying a crude all-age rate to young adults.
# ==============================================================================

if (!exists("CFG")) {
  cand <- c("R/00_config.R", "00_config.R",
            file.path(Sys.getenv("CAE_ROOT", "."), "R/00_config.R"))
  source(cand[file.exists(cand)][1])
}

suppressPackageStartupMessages({ library(tidyverse); library(readxl) })

# The full WPP2024 xlsx is 131 MB (too large for GitHub). We therefore commit a
# compact extract (3 countries, target year) and read that on a clean checkout;
# if the full xlsx is present locally we (re)build the extract from it.
extract_tsv <- file.path(CFG$dir_raw, "wpp2024_lifetable_extract.tsv")

if (file.exists(CFG$wpp_xlsx)) {
  message("[01_prepare_data] reading full WPP2024 xlsx (text mode) ...")
  raw <- read_excel(
    CFG$wpp_xlsx, sheet = "Estimates", skip = 16,
    col_types = "text", .name_repair = "minimal", guess_max = 5000
  )
  # Column positions are stable in the WPP layout (see 00 inspection).
  wpp <- tibble(
    iso3  = raw[[6]],
    year  = suppressWarnings(as.integer(raw[[11]])),
    age   = suppressWarnings(as.numeric(raw[[12]])),
    n     = suppressWarnings(as.numeric(raw[[13]])),
    mx    = suppressWarnings(as.numeric(raw[[14]])),  # central death rate m(x,n)
    qx    = suppressWarnings(as.numeric(raw[[15]])),  # probability of dying
    Lx    = suppressWarnings(as.numeric(raw[[19]])),  # person-years lived L(x,n)
    ex    = suppressWarnings(as.numeric(raw[[22]]))   # expectation of life e(x)
  ) %>%
    filter(iso3 %in% CFG$countries$iso3, year == CFG$wpp_year,
           !is.na(age), !is.na(mx))
  write_tsv(wpp, extract_tsv)        # committable, ~60 rows
} else {
  message("[01_prepare_data] full xlsx absent; reading committed extract ...")
  stopifnot(file.exists(extract_tsv))
  wpp <- read_tsv(extract_tsv, show_col_types = FALSE)
}

iso2country <- setNames(CFG$countries$country, CFG$countries$iso3)
wpp <- wpp %>% mutate(country = iso2country[iso3])

stopifnot(all(CFG$countries$iso3 %in% wpp$iso3))

# ---- Per-country abridged life-table TSV ------------------------------------
walk(CFG$countries$iso3, function(code) {
  ctry <- iso2country[[code]]
  lt <- wpp %>% filter(iso3 == code) %>%
    transmute(Age = age, Interval = n, Death_Rate = mx, q = qx,
              Person_Years = Lx, Life_Exp = ex, Sex = "Both") %>%
    arrange(Age)
  write_tsv(lt, file.path(CFG$dir_derived, sprintf("life_table__%s.tsv", ctry)))
})

# ---- Age-varying mortality schedule for the representative cohort ------------
cycle_mid_age <- CFG$index_age + (seq_len(CFG$n_cycles) - 0.5) * CFG$cycle_years
mortality_schedule <- wpp %>%
  group_by(country, iso3) %>%
  group_modify(~ tibble(
    cycle = seq_len(CFG$n_cycles),
    age_mid = cycle_mid_age,
    bg_rate = approx(.x$age, .x$mx, xout = cycle_mid_age, rule = 2)$y,
    residual_le = approx(.x$age, .x$ex, xout = cycle_mid_age, rule = 2)$y
  )) %>%
  ungroup() %>%
  arrange(match(country, CFG$countries$country), cycle)
write_tsv(mortality_schedule,
          file.path(CFG$dir_derived, "mortality_schedule.tsv"))

# Retain a compact age-30 table for human inspection and backwards-compatible
# downstream diagnostics. It is no longer the source of transition mortality.
bg <- mortality_schedule %>%
  filter(cycle == 1) %>%
  select(country, iso3, bg_rate, residual_le)
write_tsv(bg, file.path(CFG$dir_derived, "background_mortality.tsv"))

# ---- Prevalence + population table (from config) -----------------------------
prev <- CFG$countries %>%
  mutate(population = CFG$population[country],
         gdp_pc     = CFG$gdp_pc[country],
         cost_annual = CFG$cost_annual[country]) %>%
  left_join(bg %>% select(country, bg_rate), by = "country")

write_tsv(prev, file.path(CFG$dir_derived, "prevalence.tsv"))

message("[01_prepare_data] background annual death rate per country:")
print(mortality_schedule %>%
        select(country, cycle, age_mid, bg_rate) %>%
        mutate(bg_rate = round(bg_rate, 5)), n = 30)
message(sprintf("[01_prepare_data] wrote %d life tables + mortality schedule + prevalence",
                nrow(CFG$countries)))
