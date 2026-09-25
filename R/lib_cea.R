# ==============================================================================
# lib_cea.R
#
# Willem M. (Wim) Otte, w.m.otte@umcutrecht.nl
#
# Shared valuation primitives for the cost-effectiveness accounting.
#
# Willem M. (Wim) Otte, w.m.otte@umcutrecht.nl
#
# 04_dalys_costs_cea.R (probabilistic accounting), 07_sensitivity.R
# (deterministic one-way / relapse analyses) and 10_structural_scenarios.R
# (structural scenario matrix) all value the same cohort trace. Before this file
# existed each script carried its own copy of the integrator, the disability-
# weight vector and the YLL annuity. The copies drifted at least twice: once the
# disability weights were reverted in 04 only, leaving the tornado running a
# different model than the headline while every gate stayed green.
#
# Anything that turns a cohort trace into DALYs or costs belongs here, and
# nowhere else.
# ==============================================================================

if (!exists("CFG")) {
  cand <- c("R/00_config.R", "00_config.R",
            file.path(Sys.getenv("CAE_ROOT", "."), "R/00_config.R"))
  source(cand[file.exists(cand)][1])
}

# ---- Within-cycle integration -----------------------------------------------
# Composite Simpson 1/3 needs an even number of intervals. The 10-year horizon
# has five 2-year cycles, so Simpson is applied to the leading even block
# (years 0-8) and the trailing interval (years 8-10) is closed with the
# trapezoid rule. The manuscript and appendix state this hybrid explicitly
# rather than calling the whole thing Simpson's rule.
integrate_curve <- function(t, y) {
  stopifnot(length(t) == length(y), length(t) >= 2L)
  n <- length(t) - 1L                     # number of intervals
  steps <- diff(t)
  h <- steps[1]
  stopifnot(h > 0, max(abs(steps - h)) < 1e-12)
  if (n == 1L) return((h / 2) * (y[1] + y[2]))
  m <- if (n %% 2L == 0L) n else n - 1L   # leading even block
  s <- y[1] + y[m + 1L] + 4 * sum(y[seq(2, m, by = 2)])
  if (m >= 4L) s <- s + 2 * sum(y[seq(3, m - 1L, by = 2)])
  area <- (h / 3) * s
  if (m < n) area <- area + (h / 2) * (y[n] + y[n + 1])   # trapezoid tail
  area
}

# ---- Trace rewards -----------------------------------------------------------
# These functions are deliberately small. They are used by every analysis that
# turns state occupancy into health or cost outcomes, and each has a closed-form
# executable check in 06_model_checks.R.
state_column <- function(trace, state) {
  if (is.data.frame(trace)) as.numeric(trace[[state]]) else as.numeric(trace[, state])
}

discounted_yld <- function(t, trace, dw, discount) {
  occ <- as.matrix(trace[, names(dw), drop = FALSE])
  discount_at_t <- exp(-discount * t)
  integrate_curve(t, as.numeric(occ %*% dw) * discount_at_t)
}

discounted_yll <- function(t, trace, life_years_remaining, discount) {
  stopifnot(length(life_years_remaining) == length(t) - 1L)
  new_deaths <- diff(state_column(trace, "Death"))
  death_midpoints <- t[-1] - diff(t) / 2
  sum(new_deaths * life_years_remaining * exp(-discount * death_midpoints))
}

discounted_care_cost <- function(t, trace, annual_cost,
                                 untreated_cost_fraction, discount) {
  on_cost <- state_column(trace, "Treated") + state_column(trace, "Controlled") +
             state_column(trace, "Untreated") * untreated_cost_fraction
  integrate_curve(t, on_cost * annual_cost * exp(-discount * t))
}

initiation_cost <- function(comparator_gap, intervention_gap, unit_cost) {
  max(comparator_gap - intervention_gap, 0) * unit_cost
}

value_trace <- function(t, trace, dw, life_years_remaining, annual_cost,
                        untreated_cost_fraction, comparator_gap,
                        intervention_gap, initiation_unit_cost, discount) {
  yld <- discounted_yld(t, trace, dw, discount)
  yll <- discounted_yll(t, trace, life_years_remaining, discount)
  care_cost <- discounted_care_cost(t, trace, annual_cost,
                                    untreated_cost_fraction, discount)
  detect_cost <- initiation_cost(comparator_gap, intervention_gap,
                                 initiation_unit_cost)
  c(yld = yld, yll = yll, daly = yld + yll,
    cost = care_cost + detect_cost,
    care_cost = care_cost, detect_cost = detect_cost)
}

incremental_cea <- function(comparator_daly, comparator_cost,
                            intervention_daly, intervention_cost, threshold) {
  daly_averted <- comparator_daly - intervention_daly
  incremental_cost <- intervention_cost - comparator_cost
  data.frame(
    daly_averted = daly_averted,
    incr_cost = incremental_cost,
    icer = incremental_cost / daly_averted,
    inb = threshold * daly_averted - incremental_cost
  )
}

# ---- Disability weights ------------------------------------------------------
# Severity is a cohort attribute that acts on the Untreated weight only. Treated
# means "on ASM, still seizing"; in the base case both severity cohorts carry the
# 0.263 band in that state. `treated_keeps_entry_dw = TRUE` is the structural
# scenario in which a severe person who starts medication retains the severe
# weight until an explicit transition to Controlled occurs.
dw_vec <- function(p, severity, treated_keeps_entry_dw = FALSE) {
  entry <- if (severity == "severe") p$dw_sev_UT else p$dw_less_UT
  treated <- if (treated_keeps_entry_dw) entry else p$dw_less_UT
  remission <- if (is.null(p$dw_remission)) p$dw_control else p$dw_remission
  c(Untreated  = entry,
    Treated    = treated,
    Controlled = p$dw_control,
    Remission  = remission,
    Death      = 0)
}

# ---- Years of life lost ------------------------------------------------------
# Conventional DALY valuation: a death at age x is charged the discounted
# annuity of the life-table *expectation* e_x,
#
#     a(e_x) = (1 - exp(-r e_x)) / r.
#
# This is the GBD-style convention. It is NOT the expected discounted survival
# E[ integral_0^T exp(-ru) du ] = integral_0^inf exp(-ru) S(u) du. Because the
# annuity is concave in its argument, Jensen's inequality gives a(E[T]) >= E[a(T)],
# so the convention is an upper bound on expected discounted life-years.
# discounted_le_lifetable() below computes the survival-curve quantity so the two
# can be compared numerically (structural scenario YLL-LT).
discounted_annuity <- function(L, r) {
  if (r == 0) return(L)
  (1 - exp(-r * L)) / r
}

# Expected discounted life-years remaining at exact age `age`, integrating the
# abridged life table's piecewise-constant hazard:
#
#   integral_0^inf exp(-r u) S(u) du
#
# Within an interval of width w starting at elapsed time t0 with survival S0 and
# constant hazard m, the contribution is closed-form:
#
#   S0 exp(-r t0) (1 - exp(-(r+m) w)) / (r + m).
#
# `lt` is a life table with columns Age, Interval, Death_Rate (central rate m_x).
discounted_le_lifetable <- function(lt, age, r) {
  lt <- lt[order(lt$Age), ]
  # Expand the abridged table to per-interval (start, end, hazard), giving the
  # open final interval a width implied by its own hazard (1/m years).
  starts <- lt$Age
  widths <- lt$Interval
  m      <- lt$Death_Rate
  last   <- length(starts)
  if (!is.finite(widths[last]) || is.na(widths[last]) || widths[last] <= 0)
    widths[last] <- if (m[last] > 0) 1 / m[last] else 1
  ends <- starts + widths

  total <- 0
  surv  <- 1
  elapsed <- 0
  for (i in seq_along(starts)) {
    if (ends[i] <= age) next                      # entirely before entry age
    lo <- max(starts[i], age)
    w  <- ends[i] - lo
    if (w <= 0) next
    mi <- m[i]
    denom <- r + mi
    contrib <- if (denom > 0) {
      surv * exp(-r * elapsed) * (1 - exp(-denom * w)) / denom
    } else {
      surv * exp(-r * elapsed) * w
    }
    total <- total + contrib
    surv  <- surv * exp(-mi * w)
    elapsed <- elapsed + w
    if (surv < 1e-10) break
  }
  total
}


# ---- PSA sampling ------------------------------------------------------------
# Triangular distribution with lower limit, most likely value, and upper limit.
# The reported bounds are treated as plausible limits, not posterior intervals.
rtri <- function(n, lo, mode, hi) {
  stopifnot(lo <= mode, mode <= hi, lo < hi)
  u <- runif(n)
  qtri(u, lo, mode, hi)
}
qtri <- function(u, lo, mode, hi) {
  cut <- (mode - lo) / (hi - lo)
  ifelse(u < cut,
         lo + sqrt(u * (hi - lo) * (mode - lo)),
         hi - sqrt((1 - u) * (hi - lo) * (hi - mode)))
}
ptri <- function(x, lo, mode, hi) {
  ifelse(x <= mode,
         (x - lo)^2 / ((hi - lo) * (mode - lo)),
         1 - (hi - x)^2 / ((hi - lo) * (hi - mode)))
}
# Re-map a triangular draw onto a different triangular distribution at the same
# quantile. Sensitivity scenarios that change one distribution then keep their
# pairing with every other sampled parameter.
remap_tri <- function(x, from, to) {
  qtri(ptri(x, from[1], from[2], from[3]), to[1], to[2], to[3])
}

# One row per Monte Carlo draw. Shared clinical constructs use shared draws.
# Bounds come from CFG so the sampler and the verification checks cannot disagree.
sample_psa_params <- function(n) {
  tibble(
    iter        = seq_len(n),
    smr_U       = rtri(n, CFG$smr_sens$Untreated[1],
                       CFG$smr$Untreated, CFG$smr_sens$Untreated[2]),
    smr_T       = rtri(n, CFG$smr_sens$Treated[1],
                       CFG$smr$Treated,   CFG$smr_sens$Treated[2]),
    # Controlled and Remission are one seizure-free construct with one draw.
    smr_SF      = rtri(n, CFG$smr_sens$SeizureFree[1],
                       CFG$smr$Controlled, CFG$smr_sens$SeizureFree[2]),
    spont_rem   = rtri(n, CFG$spont_remission_sens[1],
                       CFG$spont_remission_prob, CFG$spont_remission_sens[2]),
    control_cum_LMIC = rtri(n, CFG$control_cum_sens$LMIC[1],
                            CFG$control_cum_10y$LMIC, CFG$control_cum_sens$LMIC[2]),
    control_cum_HIC  = rtri(n, CFG$control_cum_sens$HIC[1],
                            CFG$control_cum_10y$HIC, CFG$control_cum_sens$HIC[2]),
    dw_sev_UT   = rtri(n, dw_bounds("Untreated_severe")[1],
                       CFG$dw$Untreated_severe, dw_bounds("Untreated_severe")[2]),
    dw_less_UT  = rtri(n, dw_bounds("Untreated_less_severe")[1],
                       CFG$dw$Untreated_less_severe,
                       dw_bounds("Untreated_less_severe")[2]),
    dw_control  = rtri(n, dw_bounds("Controlled")[1],
                       CFG$dw$Controlled, dw_bounds("Controlled")[2]),
    cost_mult   = rtri(n, CFG$cost_mult_sens[1], 1.00, CFG$cost_mult_sens[2]),
    case_find_cost = rtri(n, CFG$cost_case_finding_sens[1],
                          CFG$cost_case_finding, CFG$cost_case_finding_sens[2]),
    untreated_cost_frac = rtri(n, CFG$cost_untreated_frac_psa[1],
                               CFG$cost_untreated_frac,
                               CFG$cost_untreated_frac_psa[2])
  )
}

message(paste("[lib_cea] loaded integration, trace-reward, disability-weight,",
              "YLL, incremental CEA, and PSA sampling functions"))
