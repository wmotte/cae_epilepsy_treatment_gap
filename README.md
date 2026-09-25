# Reducing the epilepsy treatment gap: health benefits and costs

This repository holds the code and data behind our paper:

> Otte WM, Keezer MR, Sander JW, Singh G. Health benefits and costs of reducing
> the epilepsy treatment gap in Nigeria, Ecuador, and the UK: a Markov modelling
> study. *Submitted to The Lancet Neurology*, 2026.

Most people with epilepsy in low-income and middle-income countries never
receive antiseizure medication, even though cheap drugs stop seizures in many of
them. We wanted to know what closing part of that gap would buy in health, what
it would cost, and how sure we can be. We built a Markov cohort model and ran it
for three very different settings: Nigeria, Ecuador and the UK.

Everything in the paper's tables and figures comes out of the scripts here. One
command rebuilds the lot from the raw UN life tables in about four minutes.

## Quick start

You need R (we used 4.6.1). From the repository root:

```bash
Rscript -e 'install.packages("renv"); renv::restore()'   # first time only
Rscript R/run_all.R
```

`renv::restore()` installs the exact package versions recorded in `renv.lock`.
If you would rather use your own library, installing `tidyverse`, `readxl` and
`scales` is enough, although numbers may then differ in the last decimal.

`run_all.R` works from any directory. If you source scripts one at a time from
somewhere else, set `CAE_ROOT` to the repository path first.

Ghostscript is optional. When `cairo_pdf` is not available, the figure code uses
it to embed fonts in the PDF figures. Without either you still get every figure,
with a warning that the fonts were not embedded.

## What the model does

The model follows a cohort of people with active epilepsy, entering at age 30,
over 10 years in five 2-year cycles. Each person is in one of five states:

| State | On medication | Having seizures |
|---|---|---|
| Untreated | no | yes |
| Treated | yes | yes |
| Controlled | yes | no |
| Remission | no | no |
| Death | | |

Untreated people can go into spontaneous remission. Treated people can become
seizure-free. Everyone faces their country's background mortality, taken from
the UN World Population Prospects 2024 and updated for age at every cycle, and
multiplied by a state-specific standardised mortality ratio (SMR). There is no
relapse in the base case. A relapse scenario is reported separately.

The treatment gap is the untreated share of the cohort when the model starts.
The intervention moves part of that share into Treated once, at entry. In
Nigeria and Ecuador the gap falls from 80% to 20%. In the UK it falls from 10% to
2.5%, the same relative reduction.

We count disability-adjusted life-years (DALYs), made of years lived with
disability plus years of life lost, and direct health-care costs. Both are
discounted at 3.5% a year. We then compare cost per DALY averted with a threshold
of half of GDP per capita, which approximates the health a health system gives
up elsewhere when it spends the money. We also report 1 and 3 times GDP per
capita, the older WHO yardsticks, so the results can be compared with earlier
studies.

Uncertainty comes from 1000 draws of every uncertain input, each from a
triangular distribution. The status-quo and intervention arms share the same
draws, so the difference between them is not blurred by sampling noise. The seed
is fixed at 20260520.

We run two cohorts in each country. The severe cohort has at least one seizure a
month when untreated, and the less-severe cohort has between one and eleven a
year. We use drug resistance as a stand-in for severity and weight the two
cohorts 13.7% to 86.3%, following the population-based estimate of Sultana and
colleagues (2021).

## Main inputs

| Input | Base case | Source |
|---|---|---|
| Disability weights | 0.552 severe untreated, 0.263 less-severe untreated and treated with seizures, 0.049 seizure-free | Salomon et al., 2015 (GBD 2013) |
| SMR, untreated | 6.3 (95% CI 2.0–10.0) | Carpio et al., 2005 |
| SMR, treated with seizures | 2.54 (1.84–3.44) | Mohanraj et al., 2006 |
| SMR, seizure-free | 0.95 (0.68–1.29) | Mohanraj et al., 2006 |
| Seizure control by 10 years | 50% in Nigeria and Ecuador, 70% in the UK | Scenario input |
| Spontaneous remission by 10 years | 30% | Nicoletti et al., 2009 |
| Annual cost of treated epilepsy | US$654.59 Nigeria, US$354.20 Ecuador, US$2320.76 UK (2019 US$) | Begley et al., 2022 |
| Case finding and starting treatment | US$50 per extra person treated | Our assumption |
| Background mortality | Age-specific, 2023 | UN World Population Prospects 2024 |

The full table, with ranges and notes on each source, is written by the pipeline
to `output/tables/table_parameters.tsv`. All values live in `R/00_config.R`,
which is also where to change them. The comments there explain each choice,
including the weak points. The untreated SMR is the most important of those. No
study has measured it in an untreated group, and the results in Nigeria depend
on it.

## Main results

For people with severe epilepsy, narrowing the gap averted 2.85 DALYs per person
in Nigeria, 1.85 in Ecuador and 0.22 in the UK. At half of GDP per capita it was
good value in Ecuador and the UK. In Nigeria the severe cohort came out at
US$977 per DALY averted against a threshold of US$1000, and the
probability of being cost-effective was 54%. Nigeria's result turned on the
untreated SMR, which had to exceed 5.8 for the severe cohort to be
cost-effective. `output/tables/smr_breakeven.tsv` gives the break-even values
for every cohort.

## Repository layout

```text
R/                 analysis scripts, run in order by run_all.R
data/raw/          UN WPP 2024 life-table extract for the three countries
data/derived/      life tables, mortality schedule and prevalence built from it
output/tables/     every results table, as tab-separated text
output/figures/    the paper's figures, each as a 600 dpi PNG and a vector PDF
output/diagnostics/  a cost-effectiveness plane that is not in the paper
renv.lock          exact R package versions
```

## The scripts

| Script | What it does |
|---|---|
| `00_config.R` | Every parameter, path and seed. Change inputs here. |
| `01_prepare_data.R` | Builds age-specific background mortality and residual life expectancy from the WPP extract |
| `02_markov_model.R` | Builds the transition matrix for each cycle and runs the cohort. Checks that every row sums to one. |
| `03_run_simulations.R` | Runs the 1000 probabilistic draws for every country, severity and gap |
| `04_dalys_costs_cea.R` | Turns cohort traces into DALYs and costs, then ICERs, net benefit and probabilities |
| `05_figures_tables.R` | Figures 1, 2 and 4, the parameter table and the diagnostic plane |
| `06_model_checks.R` | Sixteen automatic checks on the model's structure, sampling and traces |
| `07_sensitivity.R` | One-way (tornado) analysis, seizure-control threshold, relapse scenario and acceptability curves |
| `08_markov_diagram.R` | State diagram and transition matrix figures |
| `09_descriptive_addenda.R` | State occupancy, national cohort sizes, the combined population and the Ecuador and UK cohort figures |
| `10_structural_scenarios.R` | Alternative ways of mapping the evidence to model states, and a split of DALYs into disability and early death |
| `11_additional_sensitivity.R` | Alternative mortality evidence, a declining untreated SMR, other UK gaps, other severity mixes, the SMR grid and break-even, and a 10 000-draw check |
| `lib_cea.R` | Shared code for discounting, disability weights and the integration within cycles |
| `plot_theme.R` | Figure style (Okabe–Ito colours) and the PNG and PDF export |
| `run_all.R` | Runs everything in order and records the software versions |

Integration within cycles uses Simpson's one-third rule over the first eight
years and the trapezoid rule for the last two, because Simpson's rule needs an
even number of intervals.

## Where each table and figure comes from

| In the paper | File |
|---|---|
| Figure 1 | `output/figures/fig1_state_occupancy` |
| Figure 2 | `output/figures/fig2_dalys_averted` |
| Figure 3 | `output/figures/fig3_tornado` |
| Figure 4 | `output/figures/fig4_ceac` |
| Table 2 | `cea_results.tsv` and `prob_cost_effective.tsv` |
| Figures S1 and S2 | `figS1_markov_diagram`, `figS2_transition_matrix` |
| Figure S3 | `figS3_threshold_control` |
| Figure S4 | `figS4_combined_dalys` |
| Figures S5 and S6 | `figS5_cohort_ecuador`, `figS6_cohort_uk` |
| Table S1 | `table_parameters.tsv` |
| Table S2 | `mortality_table.tsv` |
| Table S3 | `cea_results.tsv` and `prob_cost_effective.tsv` |
| Table S4 | `sensitivity_oneway.tsv` |
| Table S5 | `sensitivity_relapse.tsv` |
| Table S6 | `state_occupancy_snapshots.tsv` |
| Table S7 | `cohort_sizes.tsv` |
| Table S8 | `model_checks.tsv` |
| Table S9 | `mc_error.tsv` |
| Table S10 | `structural_scenarios.tsv` |
| Table S11 | `daly_decomposition.tsv` |
| Tables S12 and S14 | `sensitivity_additional.tsv` |
| Table S13 | `mc_10000.tsv` |
| Table S15 | `sensitivity_smr_grid.tsv` |
| Table S16 | `smr_breakeven.tsv` |

Tables are in `output/tables/`. Figures are in `output/figures/`, and each one has
a `.png` and a `.pdf`. The draw-level results are in `cea_draws.tsv` and
`psa_draws_bridged.tsv` if you want to run your own analyses on them.
`session_provenance.txt` lists the R version, platform and package versions of
the run that produced these files.

## Data

`data/raw/wpp2024_lifetable_extract.tsv` is a small extract of the UN World
Population Prospects 2024 abridged life table (both sexes, 2023) for Nigeria,
Ecuador and the UK. The full workbook is 131 MB, too large for GitHub. If you
download it from <https://population.un.org/wpp/> and put it in `data/raw/`,
`01_prepare_data.R` rebuilds the extract from it. Otherwise the committed extract
is used. `data/raw/README.md` has the exact file name.

The simulation traces (`data/derived/sim_traces.rds`, about 7 MB) are not
committed. `run_all.R` writes them on the first run.

Everything else, including costs, disability weights, SMRs and prevalence, comes
from published studies. The values are in `R/00_config.R`, and the sources are
listed in `output/tables/table_parameters.tsv` and in the paper's appendix.

## Reporting

The analysis is reported according to CHEERS 2022. Because it uses Global Burden
of Disease disability weights, it also follows GATHER. Both checklists are in the
paper's appendix.

## Citation

If you use this code or its results, please cite the paper. Until it is
published, cite this repository. `CITATION.cff` has the details, and GitHub shows
them under "Cite this repository".

## Licence

The code is released under the MIT licence (`LICENSE`). The derived data, results
tables and figures are released under CC BY 4.0 (`LICENSE-data`). The UN World
Population Prospects data are © United Nations and are used under CC BY 3.0 IGO.

## Contact

Wim Otte, University Medical Center Utrecht, w.m.otte@umcutrecht.nl.
Questions and bug reports are welcome as GitHub issues.
