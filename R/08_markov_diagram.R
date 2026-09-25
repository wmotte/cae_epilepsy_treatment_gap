# ==============================================================================
# 08_markov_diagram.R
#
# Willem M. (Wim) Otte, w.m.otte@umcutrecht.nl
#
# Supplementary figures describing the five-state Markov model:
#   figS1_markov_diagram.png    state-transition diagram (publication layout)
#   figS2_transition_matrix.png the corresponding per-cycle transition matrix,
#                               numeric, status-quo vs bridged (Nigeria)
# Pure ggplot2 (no extra dependencies). Arrow/parameter labels and the matrix
# numbers are produced from CFG + the actual model code, so the figures cannot
# drift out of sync with the analysis.
# ==============================================================================

if (!exists("CFG")) {
  cand <- c("R/00_config.R", "00_config.R",
            file.path(Sys.getenv("CAE_ROOT", "."), "R/00_config.R"))
  source(cand[file.exists(cand)][1])
}
source(file.path(CFG$root, "R", "02_markov_model.R"))   # build_transition_matrix()
if (!exists("theme_cae")) source(file.path(CFG$root, "R", "plot_theme.R"))
suppressPackageStartupMessages({ library(tidyverse) })

# ==============================================================================
# FIGURE S1 — state-transition diagram
# ==============================================================================
# Layout: care pathway runs left->right along the top (Untreated -> Treated ->
# Controlled); spontaneous Remission sits below Untreated; Death is the
# absorbing sink at bottom centre.
# Disability-weight captions are generated from CFG$dw so the figure cannot
# disagree with the model. A hand-typed earlier version showed "DW 0.263-0.552"
# for Treated (the code uses a single 0.263) and "DW 0" for Remission (the code
# uses 0.049), so figure S1 depicted a model the results never came from.
# Untreated carries two categorical weights, one per severity cohort, so they are
# listed separately instead of as a continuous range.
.dw_lab <- function(x) sub("^0", "0", sprintf("%.3f", x))
dw_untreated  <- sprintf("DW %s severe; %s less severe",
                         .dw_lab(CFG$dw$Untreated_severe),
                         .dw_lab(CFG$dw$Untreated_less_severe))
dw_treated    <- sprintf("DW %s (both cohorts)", .dw_lab(CFG$dw$Treated))
dw_controlled <- sprintf("DW %s", .dw_lab(CFG$dw$Controlled))
dw_remission  <- sprintf("DW %s", .dw_lab(CFG$dw$Remission))

states <- tribble(
  ~id,          ~x,    ~y,    ~name,         ~desc,                                 ~dw,
  "Untreated",  0.0,   3.0,   "UNTREATED",   "active epilepsy, no or inadequate ASM", dw_untreated,
  "Treated",    4.2,   3.0,   "TREATED",     "on ASM, still seizing",                 dw_treated,
  "Controlled", 8.4,   3.0,   "CONTROLLED",  "on ASM, seizure-free",                  dw_controlled,
  "Remission",  0.0,  -0.2,   "REMISSION",   "off ASM, seizure-free",                 dw_remission,
  "Death",      4.2,  -2.6,   "DEATH",       "absorbing state",                       "DW -"
) %>% mutate(fill = unname(state_colours[id]),
             label = sprintf("%s\n(%s)\n%s", name, desc, dw))

# box half-extent used for trimming arrows to the box edge
BW <- 1.55   # half-width
BH <- 0.62   # half-height

xy <- function(id) states %>% filter(.data$id == !!id) %>% select(x, y) %>% as.list()

# Segment from box edge of `from` to box edge of `to` (rectangular boxes, so
# trim along the line by the box half-extent projected on the direction).
seg <- function(from, to) {
  a <- xy(from); b <- xy(to)
  ang <- atan2(b$y - a$y, b$x - a$x)
  trim <- function(hw, hh) {
    # distance from centre to rectangle edge along angle `ang`
    cx <- abs(cos(ang)); cy <- abs(sin(ang))
    min(if (cx > 1e-6) hw / cx else Inf, if (cy > 1e-6) hh / cy else Inf)
  }
  t <- trim(BW, BH)
  tibble(x = a$x + t * cos(ang), y = a$y + t * sin(ang),
         xend = b$x - t * cos(ang), yend = b$y - t * sin(ang), ang = ang)
}

# label placed at the segment midpoint, offset PERPENDICULAR to the arrow so it
# never sits on top of the line (the masking bug in the first version)
mid <- function(d, off = 0) {
  d %>% mutate(nx = -sin(ang), ny = cos(ang),
               mx = (x + xend) / 2 + off * nx,
               my = (y + yend) / 2 + off * ny)
}

# Parameter labels (in sync with the model via CFG).
lab_uptake  <- "coverage intervention at model entry\n(1 - treatment gap)"
lab_control <- sprintf("achieve control\n%.1f%% (LMIC) / %.1f%% (HIC) per cycle",
                       100 * CFG$control_prob$LMIC, 100 * CFG$control_prob$HIC)
lab_spont   <- sprintf("spontaneous\nremission\n~%.1f%%/cycle",
                       100 * CFG$spont_remission_prob)

flows <- bind_rows(
  mid(seg("Untreated",  "Treated"),    off =  0.95) %>% mutate(lab = lab_uptake, kind = "entry"),
  mid(seg("Treated",    "Controlled"), off =  0.95) %>% mutate(lab = lab_control, kind = "transition"),
  mid(seg("Untreated",  "Remission"),  off = -1.05) %>% mutate(lab = lab_spont, kind = "transition")
)
deaths <- bind_rows(
  mid(seg("Untreated",  "Death"), off = -0.55) %>% mutate(lab = sprintf("SMR %.1f",  CFG$smr$Untreated)),
  mid(seg("Treated",    "Death"), off =  0.45) %>% mutate(lab = sprintf("SMR %.1f",  CFG$smr$Treated)),
  mid(seg("Controlled", "Death"), off =  0.55) %>% mutate(lab = sprintf("SMR %.1f",  CFG$smr$Controlled)),
  mid(seg("Remission",  "Death"), off = -0.50) %>% mutate(lab = sprintf("SMR %.2f", CFG$smr$Remission))
)

# self-loops (stay in same state next cycle) — small arc above each transient box
selfloop <- function(id) {
  a <- xy(id)
  tibble(x = a$x - 0.55, y = a$y + BH + 0.02,
         xend = a$x + 0.55, yend = a$y + BH + 0.02)
}
loops <- bind_rows(lapply(c("Untreated", "Treated", "Controlled", "Remission"), selfloop))

arrow_std   <- arrow(length = unit(0.13, "inches"), type = "closed")
arrow_thin  <- arrow(length = unit(0.09, "inches"), type = "closed")
arrow_open  <- arrow(length = unit(0.12, "inches"), type = "open", angle = 22)

p1 <- ggplot() +
  annotate("text", x = 4.2, y = xy("Untreated")$y + BH + 1.75, size = 2.8,
           colour = "grey40", fontface = "italic",
           label = "composite loss-to-care and recurrence path (scenario only): treated, controlled, or remission to untreated") +
  # self-loops
  geom_curve(data = loops, aes(x, y, xend = xend, yend = yend),
             curvature = -1.4, ncp = 12, linewidth = 0.5, colour = "grey60",
             arrow = arrow_thin) +
  # mortality arrows (thin grey) -> Death
  geom_segment(data = deaths, aes(x, y, xend = xend, yend = yend),
               linewidth = 0.6, colour = "grey50", arrow = arrow_thin) +
  geom_text(data = deaths, aes(mx, my, label = lab), size = 2.7, colour = "grey35") +
  # base transition arrows (thick blue)
  geom_segment(data = filter(flows, kind == "transition"), aes(x, y, xend = xend, yend = yend),
               linewidth = 1.15, colour = unname(bc_palette["blue"]), arrow = arrow_std) +
  geom_segment(data = filter(flows, kind == "entry"), aes(x, y, xend = xend, yend = yend),
               linewidth = 1.0, linetype = "dashed",
               colour = unname(bc_palette["blue"]), arrow = arrow_open) +
  geom_label(data = flows, aes(mx, my, label = lab), size = 2.85, colour = unname(bc_palette["blue"]),
             linewidth = 0, fill = "white", fontface = "bold",
             label.padding = unit(0.12, "lines"), lineheight = 0.92) +
  # state boxes drawn last
  geom_label(data = states, aes(x, y, label = label, fill = fill),
             size = 3.0, label.r = unit(0.30, "lines"), linewidth = 0.5,
             label.padding = unit(0.5, "lines"), lineheight = 0.95, colour = "grey15") +
  # relapse arc drawn LAST so its arrowhead sits above the self-loops
  geom_curve(aes(x = xy("Controlled")$x - 0.3, y = xy("Controlled")$y + BH + 0.95,
                 xend = xy("Untreated")$x + 0.3, yend = xy("Untreated")$y + BH + 0.95),
             curvature = 0.08, ncp = 20, linetype = "dashed", linewidth = 0.7,
             colour = "grey45", arrow = arrow_open) +
  scale_fill_identity() +
  coord_equal(clip = "off") +
  scale_x_continuous(expand = expansion(mult = 0.10)) +
  scale_y_continuous(expand = expansion(mult = 0.10)) +
  theme_void(base_size = 11) +
  theme(plot.margin = margin(12, 18, 12, 18))

save_fig(p1, "figS1_markov_diagram", width = 10.0, height = 6.6)
message("[08_markov_diagram] wrote figS1_markov_diagram.png + .pdf")

# ==============================================================================
# FIGURE S2 — the corresponding per-cycle transition matrix (numeric)
# ==============================================================================
# Built from the SAME build_transition_matrix() used in the analysis, for
# Nigeria at the first cycle midpoint. Status quo and intervention use the same
# transition matrix: their difference is the initial Untreated/Treated vector.
bg_ng <- read_tsv(file.path(CFG$dir_derived, "mortality_schedule.tsv"),
                  show_col_types = FALSE) %>%
  filter(country == "Nigeria", cycle == 1) %>% pull(bg_rate)

p_base <- list(
  bg_rate      = bg_ng,
  smr          = CFG$smr,
  spont_rem    = CFG$spont_remission_prob,
  control_prob = CFG$control_prob$LMIC,
  relapse      = 0,
  treat_uptake = CFG$background_treatment_uptake
)

matrix_to_tiles <- function(M, scenario) {
  as.data.frame(M) %>%
    rownames_to_column("from") %>%
    pivot_longer(-from, names_to = "to", values_to = "p") %>%
    mutate(scenario = scenario,
           from = factor(from, levels = CFG$states),
           to   = factor(to,   levels = CFG$states))
}

M <- build_transition_matrix(p_base)
tiles <- matrix_to_tiles(M, "Per-cycle transition matrix") %>%
  mutate(txt = ifelse(p < 1e-6, "0", formatC(p, format = "f", digits = 3)),
         is_zero = p < 1e-6,
         diag = as.character(from) == as.character(to))

p2 <- ggplot(tiles, aes(to, fct_rev(from))) +
  geom_tile(aes(fill = p), colour = "grey80", linewidth = 0.4) +
  geom_tile(data = filter(tiles, diag), fill = NA, colour = "grey25", linewidth = 0.9) +
  geom_text(aes(label = txt, colour = p > 0.5), size = 3.0, show.legend = FALSE) +
  scale_fill_gradient(low = "white", high = unname(bc_palette["blue"]), limits = c(0, 1),
                      name = "Per-cycle probability",
                      guide = guide_colourbar(title.position = "top",
                                              barwidth = unit(3.2, "cm"))) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey20")) +
  scale_x_discrete(position = "top") +
  coord_equal() +
  labs(x = NULL, y = NULL) +
  theme_cae(base_size = 11) +
  theme(strip.text = element_text(face = "bold", size = 10, hjust = 0),
        axis.text.x = element_text(angle = 30, hjust = 0, face = "bold"),
        axis.text.y = element_text(face = "bold"),
        axis.line = element_blank(), axis.ticks = element_blank(),
        legend.position = "top",
        plot.margin = margin(12, 14, 12, 14))

save_fig(p2, "figS2_transition_matrix", width = 8.5, height = 6.2)
message("[08_markov_diagram] wrote figS2_transition_matrix.png + .pdf")
