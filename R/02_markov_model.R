# ==============================================================================
# 02_markov_model.R
#
# Willem M. (Wim) Otte, w.m.otte@umcutrecht.nl
#
# 5-state Markov cohort state-transition model (tidyverse / base-matrix core).
#
# States (fixed order, see 00_config.R):
#   Untreated -> Treated -> Controlled -> Remission -> Death(absorbing)
#
# Exposes:
#   build_transition_matrix(p)  -> 5x5 matrix, rows sum to 1 (asserted)
#   run_cohort(p)               -> tidy tibble: cycle, state, probability
#
# `p` is a named list of resolved per-cycle parameters produced by the caller
# (03_run_simulations.R), so this file is pure and side-effect free.
# ==============================================================================

if (!exists("CFG")) {
  cand <- c("R/00_config.R", "00_config.R",
            file.path(Sys.getenv("CAE_ROOT", "."), "R/00_config.R"))
  source(cand[file.exists(cand)][1])
}
suppressPackageStartupMessages(library(tidyverse))

# Build the per-cycle transition matrix from resolved parameters.
#   p$bg_rate        annual background death rate
#   p$smr            named vector of SMRs (Untreated/Treated/Controlled/Remission)
#   p$treat_uptake   background per-cycle treatment initiation (base case 0)
#   p$spont_rem      per-cycle spontaneous remission share (Untreated survivors)
#   p$control_prob   per-cycle prob a Treated survivor becomes Controlled
build_transition_matrix <- function(p) {
  s  <- CFG$states
  M  <- matrix(0, length(s), length(s), dimnames = list(s, s))
  pd <- function(state) rate_to_prob(p$bg_rate * p$smr[[state]], CFG$cycle_years)
  # V2: optional per-cycle relapse to Untreated (secondary treatment gap).
  # Defaults to 0 so the base case is unchanged from the no-relapse structure.
  relapse <- if (is.null(p$relapse)) 0 else p$relapse

  # --- Untreated ---
  pdU <- pd("Untreated")
  ptU <- (1 - pdU) * p$treat_uptake                  # -> Treated
  prU <- (1 - pdU - ptU) * p$spont_rem               # -> Remission (spontaneous)
  M["Untreated", "Death"]     <- pdU
  M["Untreated", "Treated"]   <- ptU
  M["Untreated", "Remission"] <- prU
  M["Untreated", "Untreated"] <- 1 - pdU - ptU - prU

  # --- Treated (on ASM, seizing) ---
  pdT <- pd("Treated")
  pcT <- (1 - pdT) * p$control_prob                  # -> Controlled
  puT <- (1 - pdT - pcT) * relapse                   # -> Untreated (relapse)
  M["Treated", "Death"]      <- pdT
  M["Treated", "Controlled"] <- pcT
  M["Treated", "Untreated"]  <- puT
  M["Treated", "Treated"]    <- 1 - pdT - pcT - puT

  # --- Controlled (on ASM, seizure-free) ---
  pdC <- pd("Controlled")
  puC <- (1 - pdC) * relapse                         # -> Untreated (relapse)
  M["Controlled", "Death"]      <- pdC
  M["Controlled", "Untreated"]  <- puC
  M["Controlled", "Controlled"] <- 1 - pdC - puC

  # --- Remission (off ASM, seizure-free) ---
  pdR <- pd("Remission")
  puR <- (1 - pdR) * relapse                         # -> Untreated (relapse)
  M["Remission", "Death"]     <- pdR
  M["Remission", "Untreated"] <- puR
  M["Remission", "Remission"] <- 1 - pdR - puR

  # --- Death (absorbing) ---
  M["Death", "Death"] <- 1

  # invariants
  stopifnot(all(M >= -1e-9), all(abs(rowSums(M) - 1) < 1e-9))
  M
}

# Run the cohort forward from an initial prevalence treatment gap. Scenario
# coverage is applied once, before cycle 1: gap are Untreated and 1-gap Treated.
# `p$bg_rates` can vary by cycle as the representative cohort ages.
run_cohort <- function(p) {
  s <- CFG$states
  initial_gap <- if (is.null(p$initial_gap)) 1 else p$initial_gap
  stopifnot(initial_gap >= 0, initial_gap <= 1)
  x0 <- setNames(c(initial_gap, 1 - initial_gap, 0, 0, 0), s)
  bg_rates <- if (!is.null(p$bg_rates)) p$bg_rates else rep(p$bg_rate, CFG$n_cycles)
  stopifnot(length(bg_rates) == CFG$n_cycles)

  trace <- matrix(0, CFG$n_cycles + 1L, length(s), dimnames = list(NULL, s))
  trace[1, ] <- x0
  x <- x0
  for (k in seq_len(CFG$n_cycles)) {
    pk <- p
    pk$bg_rate <- bg_rates[k]
    M <- build_transition_matrix(pk)
    x <- as.numeric(x %*% M); names(x) <- s
    trace[k + 1L, ] <- x
  }

  as_tibble(trace) %>%
    mutate(cycle = row_number() - 1L) %>%
    pivot_longer(all_of(s), names_to = "state", values_to = "prob") %>%
    mutate(state = factor(state, levels = s))
}

message("[02_markov_model] loaded build_transition_matrix() + run_cohort()")
