# Replication package for "Encouraging Renewable Investment: Risk Sharing Using Auctions" by Konan Hara

## Overview

- The code under `code/` replicates the estimates, counterfactuals, tables, figures, and reported numerical values in the paper.
  - `data_processing/`: processes the source datasets.
  - `estimation/`: estimates the structural model and runs bootstrap and sensitivity calculations.
  - `counterfactual/`: calibrates and solves the counterfactual auction model.
  - `result/`: constructs paper-facing results.
  - `batch/`: contains the five stage-level execution files described below.
- The root `batch.sh` dispatches one replication stage at a time. The full pipeline cannot be run as one unattended command because data processing, counterfactuals, and results run locally, whereas the two estimation stages submit asynchronous Slurm jobs on an HPCC.
- Paper-ready outputs are written to `output/table/`, `output/figure/`, and `output/snt/`. Compact estimation and counterfactual results are written to `output/result/`.
- Large intermediate results and execution logs are stored outside the repository in the data directory specified by `data.txt`.
- The replication code uses R, Julia, Python, and Unix shell scripts.


## Data availability

The analysis uses data from the following publicly available sources:

- Chamber of Electric Energy Commercialization (Câmara de Comercialização de Energia Elétrica, CCEE) consolidated auction results and spot market prices;
- Energy Research Company (Empresa de Pesquisa Energética, EPE) auction participant reports, ten-year expansion plans, and power generation cost report;
- National Electric Energy Agency (Agência Nacional de Energia Elétrica, ANEEL) generation, generation-agent, and corporate-structure data;
- Central Bank of Brazil target interest rates;
- U.S. Energy Information Administration wholesale electricity data; and
- U.S. Bureau of Labor Statistics Consumer Price Index.

All underlying source data are publicly available. The assembled replication data directory referenced by `data.txt` will not be made publicly available until publication of the paper. Researchers interested in access before publication may contact the author.

### Summary of availability

- All underlying source data are publicly available.
- The assembled replication data directory will be released upon publication of the paper.
- The replication data directory contains source files under `original/` and three manually curated firm-linkage inputs under `processed/`.
- `data.txt` is user-specific and is therefore not tracked. Its first line must contain the absolute path to the replication data directory.

### Details on each data source

| Data name | Data files | Location within data directory | Provided | Citation |
| -- | -- | -- | -- | -- |
| Consolidated Auction Results Report | `generation_auction_results_202210.xlsx` | `original/Brazil/` | Yes | CCEE (2022) |
| Auction Participant Data | `generation_auction_participation.xlsx`; source reports in `EPE_auction_report_*` | `original/Brazil/` | Yes | EPE (2022a) |
| Spot Market Electricity Price Data | `monthly_average_price_200101_202212.xls` | `original/Brazil/` | Yes | CCEE (2023) |
| Short-Term Purchase Agreement Auction Descriptions | `CCEE_auction_description_LEE/*.txt` | `original/Brazil/` | Yes | CCEE (2024) |
| Power Generation Costs Report | `epe_cost/*.csv`, image files, and WebPlotDigitizer files | `original/Brazil/` | Yes | EPE (2022b) |
| Ten-Year Energy Expansion Plans | plan PDFs and `marginal_cost_operation_forecast.xlsx` | `original/Brazil/` | Yes | EPE (2026) |
| ANEEL Generation Information System | `siga-comercial-generation.csv` and data dictionary | `original/Brazil/ANEEL/` | Yes | ANEEL (2026a) |
| ANEEL Electric Power Generation Agents | `electricity-generation-agents.csv` and data dictionary | `original/Brazil/ANEEL/` | Yes | ANEEL (2025a) |
| ANEEL Corporate Structure | `corporate-structure.csv` and data dictionary | `original/Brazil/ANEEL/` | Yes | ANEEL (2025b) |
| Central Bank Target Interest Rate | `base_interest_rates_history.xlsx` | `original/Brazil/central_bank/` | Yes | Central Bank of Brazil (2023) |
| U.S. Wholesale Electricity Data | `ICE/*` | `original/US/` | Yes | EIA (2024) |
| U.S. Consumer Price Index | `CPI_2001_2022.xlsx` | `original/US/` | Yes | BLS (2023) |
| Curated firm and parent crosswalks | `generation_auction_firm.csv`, `generation_auction_parent.csv`, and `generation_auction_represent_other.csv` | `processed/` | Yes | Author construction |

The source-file preparation for each data source follows.

#### Consolidated Auction Results Report

- Download “Relatório de Resultados Consolidados dos Leilões” from CCEE (2022).
- Download date: November 2, 2022.
- The original downloaded filename was `Resultado_Consolidado_Outubro_2022.xlsx`.

#### Auction Participant Data

- Download registered- and qualified-project reports for new-energy, reserve-energy, and alternative-source auctions from EPE (2022a).
- Download date for the original compilation: November 2, 2022.
- Extract the reported project counts and capacities and compile them in `generation_auction_participation.xlsx`.

#### Spot Market Electricity Price Data

- Select “Preço Médio Mensal” and download the historical series from CCEE (2023).
- Download date: May 9, 2023.
- The original downloaded filename was `Historico_do_Preco_Medio_Mensal_-_janeiro_de_2001_a_dezembro_de_2022.xls`.

#### Short-Term Purchase Agreement Auction Descriptions

- Copy the auction descriptions and products for existing-energy auctions from CCEE (2024) and translate them into English.
- Download date: November 23, 2024.

#### Power Generation Costs Report

- Extract the median and the left half of the violin plots in Figures 3 and 8 of EPE (2022b) using [WebPlotDigitizer](https://automeris.io/WebPlotDigitizer/).
- The image and WebPlotDigitizer project files are included with the four extracted CSV files under `original/Brazil/epe_cost/`.
- Processing date: May 18, 2023.

#### Ten-Year Energy Expansion Plans

- The EPE plan PDFs are under `original/Brazil/EPE_Ten-Year Energy Expansion Plan/`.
- The marginal-cost forecasts used in the analysis are compiled in `original/Brazil/marginal_cost_operation_forecast.xlsx`.

#### ANEEL Open Data and Firm Linkages

- The generation-system snapshot was downloaded July 28, 2026; the generation-agent and corporate-structure files were downloaded February 17, 2025.
- `code/data_processing/generation_auction_firm.R` validates the supplied manually curated firm and parent crosswalks and constructs the representative parent associated with each bidder. These curated crosswalks reflect ANEEL records and the sources documented in that script.

#### Central Bank Target Interest Rate

- Save the historical interest-rate webpage from the Central Bank of Brazil (2023), import its table into a spreadsheet program, and save it as `base_interest_rates_history.xlsx`.
- Download date: May 25, 2023.

#### U.S. Wholesale Electricity Data

- The supplied files were downloaded from EIA (2024) on November 18, 2024.
- `code/data_processing/spot_us_download.sh` documents the original download commands but is not run in the replication pipeline because the remote source may change.

#### U.S. Consumer Price Index

- Download annual CPI observations from 2001 through 2022, including annual averages.
- Download date: May 25, 2023.


## Computational requirements

### Software requirements

- R 4.5.3
  - `bit64`
  - `data.table`
  - `dplyr`
  - `fixest`
  - `ggplot2`
  - `here`
  - `kableExtra`
  - `lubridate`
  - `msm`
  - `readxl`
  - `zoo`
- Julia 1.12 for local calculations and Julia 1.11.3 on the HPCC
  - `CSV`
  - `DataFrames`
  - `Distributions`
  - `FastGaussQuadrature`
  - `FiniteDifferences`
  - `Interpolations`
  - `Ipopt`
  - `JuMP`
  - `NLopt`
  - `Roots`
  - `Tables`
  - Julia standard libraries used by the code
- Python 3.12.3
  - the bundled `code/tools/xlsx2csv.py` script is used to convert U.S. wholesale-price workbooks
- Bash 3.2 or later
- Gnumeric command-line utility `ssconvert` (on macOS: `brew install gnumeric`)
- Slurm workload manager for the two estimation stages

The supplied `.sbatch` files load `Julia/1.11.3-linux-x86_64` and request AMD24 nodes. Replicators using a different cluster should adjust the Julia module name and the `--constraint` directives.

### Controlled randomness

- `code/estimation/prep.jl` sets the bootstrap, winning-probability, and simulated-distribution seed sequences used by the estimation code.
- `code/estimation/step2.jl`, `step2_dist.jl`, and `step2_sa.jl` reset the relevant seeds before their simulation and bootstrap calculations.
- `code/counterfactual/counterfactual_prep.jl` and `counterfactual_model.jl` use fixed `MersenneTwister` seeds for the primitive, objective-capacity, winning-probability, and allocation simulations.
- Conditional on the documented software and numerical solvers, repeated runs therefore use the same simulation draws.

### Runtime and hardware requirements

The current local production environment is:

- MacBook Pro with Apple M4 Max, 16 CPU cores, and 128 GB memory;
- macOS 26.6;
- four local counterfactual workers by default.

The estimation stages were run on the HPCC provided by Michigan State University's [Institute for Cyber-Enabled Research (ICER)](https://icer.msu.edu/) using Slurm jobs constrained to `amd24`. [ICER's cluster-resource documentation](https://docs.icer.msu.edu/Cluster_Resources/) describes the `amd24` class as 192-core nodes with AMD EPYC 9654 or 9684X processors and 768 GB to 2.3 TB of memory; the class also includes GPU-equipped nodes. This replication does not use GPUs: each Slurm task requests one CPU and 4 GB of memory. Step 1 and Step 2 each include 200 bootstrap array tasks and 18 sensitivity array tasks.

#### Approximate runtime

The following runtimes were measured during a clean replication run. Local stage runtimes are written automatically to `<data directory>/log/replication/runtime.tsv`; HPCC runtimes use elapsed execution time from Slurm accounting rather than queue time.

| Execution file | Environment | Approximate runtime |
| -- | -- | -- |
| `code/batch/data_processing.sh` | Local | About 15 seconds |
| `code/batch/estimation_step1.sh` | HPCC | About 30 minutes |
| `code/batch/estimation_step2.sh` | HPCC | About 10 minutes |
| `code/batch/counterfactual.sh` | Local, four workers | About 9 minutes |
| `code/batch/result.sh` | Local | About 20 seconds |

#### Storage space

Approximate storage use in the current replication environment is:

- original datasets: 898 MB;
- processed datasets: 20 MB;
- estimation and counterfactual intermediate files: 99 MB;
- repository `output/` directory: 3.3 MB; and
- paper source and compiled documents: 4.7 MB.


## Description of code

The shell entry points use Unix paths and have not been tested on Windows.

### Stage-level execution files (`code/batch/`)

- `data_processing.sh`: runs the complete local data-processing sequence.
- `estimation_step1.sh`: validates the processed estimation dataset and submits the Step 1 Slurm dependency chain.
- `estimation_step2.sh`: validates Step 1 outputs and submits the Step 2 Slurm dependency chain.
- `counterfactual.sh`: runs a fresh local counterfactual by default and writes complete results externally and compact CSVs to `output/result/`.
- `result.sh`: runs the equilibrium diagnostics and creates all paper-facing tables, figures, and sentence files.

### Data processing (`code/data_processing/`)

- `generation_auction.R`: cleans the consolidated auction results.
- `generation_auction_participation.R`: cleans registered- and qualified-project counts and capacities.
- `corporate_structure.R`: cleans ANEEL corporate-structure records.
- `generation_data.R`: cleans ANEEL project and generation-agent records.
- `generation_auction_firm.R`: validates curated firm and parent links and constructs the representative bidder parent.
- `generation_auction_lee_technology.R`: assigns technology eligibility to short-term purchase agreement auctions.
- `spot_average.R`: constructs annual and monthly Brazilian wholesale-price series.
- `spot_us_step1.sh`, `spot_us_step2.R`, and `spot_us_average.R`: convert, clean, and summarize U.S. wholesale-price and CPI data.
- `epe_cost_step2.R`: constructs the engineering-cost series from digitized EPE figures.
- `base_interest_rate.R`: cleans the Central Bank of Brazil target-rate history.
- `estimation_prep.R`: creates the final auction-level and bidder-level estimation dataset.

### Estimation Step 1 (`code/estimation/`)

- `step1.sh`: submits the Step 1 dependency chain.
- `step1.sbatch`: runs `step1_spot.jl`, `step1_prespecified.jl`, `step1_cdfdc.jl`, and the original-sample `step1.jl` estimation.
- `step1_bs.sbatch`: runs 200 firm-level bootstrap replications of `step1.jl`.
- `step1_sa.sbatch`: runs 18 Step 1 sensitivity specifications.
- `step1_result.sbatch` and `step1_result.jl`: validate and aggregate the Step 1 bootstrap results.
- `prep.jl`: defines shared data, parameters, simulation controls, and helper functions.

### Estimation Step 2 (`code/estimation/`)

- `step2.sh`: submits the Step 2 dependency chain.
- `step2_dist.sbatch`, `step2_dist.jl`, and `step2_dist_func.jl`: estimate the bid-price, capacity, and prior-started-capacity distributions and their bootstrap sampling distributions.
- `step2.sbatch`, `step2.jl`, and `step2_func.jl`: recover the cost distribution in the original sample and construct the participant primitives used by the counterfactual.
- `step2_bs.sbatch`: runs 200 Step 2 bootstrap replications.
- `step2_sa.sbatch` and `step2_sa.jl`: recover costs under the 18 sensitivity specifications.
- `step2_result.sbatch` and `step2_result.jl`: validate and aggregate the Step 2 results.

### Counterfactual (`code/counterfactual/`)

- `counterfactual.sh`: coordinates the production baseline, eight primitive scenarios, and aggregation using four workers by default.
- `counterfactual_baseline.jl` and `counterfactual_baseline.sh`: calibrate the primary representative auction and registered-project distribution.
- `counterfactual.jl` and `counterfactual_task.sh`: solve one primitive scenario over the cost-reduction and fixed-share grids.
- `counterfactual_model.jl`, `counterfactual_functions.jl`, and `counterfactual_prep.jl`: define the two-risk-type common-upper-support model and numerical routines.
- `counterfactual_result.jl` and `counterfactual_result.sh`: aggregate the complete scenario outputs and call `code/result/counterfactual.jl` to create compact repository CSVs.

### Results (`code/result/`)

- `bne_check.jl`: constructs numerical diagnostics for the sufficient monotone-equilibrium conditions.
- `descriptive.R`: creates descriptive tables, figures, and sentence files.
- `estimation.R`: creates structural-estimation tables, figures, and sentence files.
- `counterfactual.jl`: converts complete external counterfactual results into compact paper-facing CSVs.
- `counterfactual.R`: creates counterfactual tables, figures, and sentence files.


## Instructions to replicators

### Directory setup

1. Place the repository and replication data directory on the local machine. They may be in different locations.
2. Create `data.txt` in the repository root. Its first line must be the absolute path to the data directory, for example:

   ```text
   /path/to/renewable-auction-data
   ```

3. Confirm that the data directory contains `original/` and the three supplied curated crosswalks in `processed/`.
4. Install the software and packages listed above.

### 1. Data processing on the local machine

From the repository root, run:

```bash
./batch.sh data_processing
```

The stage writes processed data under `<data directory>/processed/` and records its log and runtime under `<data directory>/log/replication/`.

### 2. Estimation Step 1 on the HPCC

Synchronize the repository and the data directory's `processed/` contents to the HPCC. Create an HPCC-specific `data.txt` in the synchronized repository root. From that root, run:

```bash
./batch.sh estimation_step1
```

The command prints the IDs for the original estimation, bootstrap array, sensitivity array, and final aggregation job. Wait until the printed Step 1 result job finishes successfully before proceeding.

### 3. Estimation Step 2 on the HPCC

After Step 1 completes, run:

```bash
./batch.sh estimation_step2
```

Wait until the printed Step 2 result job finishes successfully. The final job validates all required bootstrap outputs before aggregating them.

### 4. Transfer estimation outputs to the local machine

Synchronize the following from the HPCC back to the corresponding local locations:

- the repository's `output/result/` directory; and
- the data directory's `intermediate/` directory, including `intermediate/counterfactual_task/participant_primitives.csv`.

The `output/result/` and `intermediate/` directories should come from the same completed estimation run.

### 5. Counterfactual on the local machine

Run:

```bash
./batch.sh counterfactual
```

The stage performs a fresh run by default, uses four parallel workers, stores complete task and aggregate outputs under `<data directory>/intermediate/counterfactual_task/` and `counterfactual_output/`, and writes compact CSVs to `output/result/`.

To resume or reuse an existing baseline and completed scenario outputs, use:

```bash
CF_FRESH=0 ./batch.sh counterfactual
```

The number of local workers can be changed, for example:

```bash
CF_WORKERS=6 ./batch.sh counterfactual
```

### 6. Paper outputs on the local machine

Run:

```bash
./batch.sh result
```

This stage generates numerical equilibrium diagnostics and refreshes `output/table/`, `output/figure/`, and `output/snt/`. Compiling the paper PDFs is separate from the computational replication.

### Runtime recording

- Local elapsed times are appended to `<data directory>/log/replication/runtime.tsv`.
- The HPCC submission logs list all job IDs under `<data directory>/log/replication/estimation_step1/` and `estimation_step2/`.
- For future runtime updates, use Slurm accounting to record elapsed execution time for the submitted jobs, excluding queue time. For example, `sacct -j JOB_ID --format=JobID,State,Elapsed,Start,End` reports the timing of a submitted job or array.


## List of tables and figures

Paper-ready files are in `output/table/` and `output/figure/`. Numerical text inserted with `\input{snt/...}` is in `output/snt/`.

| Figure/Table | Generating code | Output file |
| -- | -- | -- |
| Table 1: Auction summary statistics | `descriptive.R` | `summary_auction.tex` |
| Table 2: Winning-bid summary statistics | `descriptive.R` | `summary_bidder.tex` |
| Table 3: Bid-share regressions | `descriptive.R` | `preliminary_regression.tex` |
| Table 4: First-stage structural estimates | `estimation.R` | `first_stage_estimates.tex` |
| Table 5: Joint risk-aversion and cost distribution | `estimation.R` | `gamma_cost_dist.tex` |
| Table 6: Matched fixed-share comparison | `counterfactual.R` | `counterfactual_policy.tex` |
| Figure 2: Bid shares by price decile | `descriptive.R` | `bidshare_price_bin.png` |
| Figure 4: Participant risk-aversion composition | `counterfactual.R` | `counterfactual_participant_fixed_primary.png`; `counterfactual_participant_choice_primary.png` |
| Figure 5: Counterfactual auction outcomes | `counterfactual.R` | `counterfactual_choice_q_primary.png`; `counterfactual_outside_capacity_mechanisms.png`; `counterfactual_winner_high_gamma_mechanisms.png`; `counterfactual_markup_mechanisms.png` |
| Figure 6: Policymaker certainty equivalent | `counterfactual.R` | `counterfactual_pm_ce_rn_mechanisms.png`; `counterfactual_pm_ce_ra_mechanisms.png` |
| Figure C1: Brazil and U.S. spot prices | `descriptive.R` | `spotprice_annual.png`; `spotprice_monthly.png` |
| Table D1: Wholesale-price process | `estimation.R` | `spotestimates.tex` |
| Table D2: Prespecified auction distributions | `estimation.R` | `prespecdist.tex` |
| Figure D1: Bid-share goodness of fit | `estimation.R` | `bidshare_cdf_gof.png`; `bidshare_price_gof.png`; `bidshare_prior_capacity_gof.png`; `bidshare_mc_gof.png` |
| Table D3: Equilibrium bid, capacity, and prior-capacity distributions | `estimation.R` | `winprobparam.tex` |
| Figure D2: Equilibrium winning probabilities | `estimation.R` | `winprob.png` |
| Figure D3: Linear projection of costs | `estimation.R` | `c_dist_lp_t_l0.png`; `c_dist_lp_l_by_t.png` |
| Table D4: Linear projection of costs | `estimation.R` | `cost_dist_lp.tex` |
| Table D5: Structural sensitivity analysis | `estimation.R` | `sa_rp_cost_markup.tex` |

Figures 1 and 3 in the main paper are conceptual figures maintained in the paper source rather than generated from data. `counterfactual_matched_markup.png` is also generated by `counterfactual.R` for presentation slides.


## References

- ANEEL. 2025a. “Electric Power Generation Agents,” National Electric Energy Agency Open Data.
- ANEEL. 2025b. “Corporate Structure,” National Electric Energy Agency Open Data.
- ANEEL. 2026a. “SIGA—Generation Information System,” National Electric Energy Agency Open Data.
- BLS. 2023. “Consumer Price Index for All Urban Consumers,” data retrieved May 25, 2023, https://data.bls.gov/timeseries/CUUR0000SA0.
- CCEE. 2022. “Consolidated Auction Results Report” (in Portuguese), data retrieved November 2, 2022, https://www.ccee.org.br/acervo-ccee.
- CCEE. 2023. “Price Panel” (in Portuguese), data retrieved May 9, 2023, https://www.ccee.org.br/precos/painel-precos.
- CCEE. 2024. “Auction Schedule and Results” (in Portuguese), data retrieved November 23, 2024, https://www.ccee.org.br/web/guest/mercado/leilao-mercado.
- Central Bank of Brazil. 2023. “Basic Interest Rates—History” (in Portuguese), data retrieved May 25, 2023, https://www.bcb.gov.br/controleinflacao/historicotaxasjuros.
- EIA. 2024. “Wholesale Electricity and Natural Gas Market Data,” data retrieved November 18, 2024, https://www.eia.gov/electricity/wholesale/.
- EPE. 2022a. “Energy Auctions” (in Portuguese), data retrieved November 2, 2022, https://www.epe.gov.br/pt/leiloes-de-energia/leiloes.
- EPE. 2022b. “Power Generation Costs Report 2021,” EPE-DEE-RE-089/2021-r1, Energy Research Company, Brazil.
- EPE. 2026. “Ten-Year Energy Expansion Plans,” Energy Research Company, Brazil, https://www.epe.gov.br/pt/publicacoes-dados-abertos/publicacoes/plano-decenal-de-expansao-de-energia-pde.
