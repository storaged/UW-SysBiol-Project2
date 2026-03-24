# Spatiotemporal Signal Propagation — Course Materials

**Systems Biology · University of Warsaw**

This repository contains all materials for the *Spatiotemporal Signal Propagation* module —
a three-tutorial practical sequence on quantifying how biochemical signals (ERK, FoxO3A)
spread across living cell monolayers in space and time.

> **Data.** The single-cell tracking dataset is distributed exclusively through the
> course Moodle page. Download the data files and place them in the project root before
> running any scripts or notebooks.

---

## Repository structure

```
scripts/                          # Core analysis scripts
├── spatiotemporal_signal_propagation.py   # Single-block graph construction & propagation analysis
└── compare_spatiotemporal_behavior.py     # Multi-block cohort comparison across conditions

notes/                            # Lecture notes and worksheets (PDF + LaTeX source)
├── spatiotemporal_propagation_handout.pdf # Mathematical formalization handout
└── classroom_worksheet.pdf                # 3-hour guided worksheet (Parts I & II)

notebooks-part1/                  # Tutorial 1 — Exploratory data analysis
notebooks-part2/                  # Tutorial 2 — Guided spatiotemporal pipeline (Part I of worksheet)
notebooks-part3/                  # Tutorial 3 — Critical analysis & algorithmic improvements (Part II)
```

---

## Tutorials at a glance

### Tutorial 1 · Exploratory data analysis (`notebooks-part1/`)

An introduction to the dataset and the biological question.
Students load single-cell tracks, visualise ERK and FoxO3A dynamics,
characterise temporal signal patterns, and examine spatial organisation
of cells within imaging sites.

| Notebook | Topic |
|---|---|
| `00_get_data.ipynb` | Data access and preprocessing |
| `01_explore_erk_akt_signaling.ipynb` | ERK/AKT signal distributions |
| `02_temporal_dynamics.ipynb` | Temporal signal dynamics per track |
| `03_spatial_structure.ipynb` | Spatial layout and cell neighbourhoods |
| `04_collective_behavior.ipynb` | Collective activity patterns |

### Tutorial 2 · Guided spatiotemporal pipeline (`notebooks-part2/`)

Students run the full spatiotemporal graph analysis on individual
experiment–site blocks. Each notebook corresponds to one block
of the classroom worksheet (§Part I).

| Notebook | Tasks | Focus |
|---|---|---|
| `Part1_Block1_FirstRun.ipynb` | 1–2 | Run the script; inspect node and edge tables |
| `Part1_Block2_ParameterSensitivity.ipynb` | 3–5 | Effect of spatial radius, future window, jump quantile |
| `Part1_Block3_Comparison.ipynb` | 6–7 | Pairwise condition comparison; cohort-level analysis |

### Tutorial 3 · Critical analysis and improvements (`notebooks-part3/`)

Students identify three mathematical limitations of the baseline approach,
derive improved estimators on paper, and implement them in Python.
Corresponds to §Part II of the classroom worksheet.

| Notebook | Limitation | Improvement |
|---|---|---|
| `Part2_L1_LaggedExposure.ipynb` | L1 — simultaneity | Lagged exposure $E_k^{(\tau)}(t)$; lag curve $\mathrm{RR}^{(\tau)}$ |
| `Part2_L2_GlobalThreshold.ipynb` | L2 — non-comparable thresholds | Global cross-block threshold $\theta^*$ |
| `Part2_L3_DoseResponse.ipynb` | L3 — binary exposure | Fractional exposure $\tilde{E}_k(t)$; dose–response analysis |

---

## Notes and handouts

| File | Description |
|---|---|
| `notes/spatiotemporal_propagation_handout.pdf` | Self-contained mathematical formalization of the pipeline: graph construction, jump detection, exposure and propagation estimators, risk measures |
| `notes/classroom_worksheet.pdf` | 3-hour worksheet with guided tasks (Part I) and open-ended improvement exercises (Part II) |

The LaTeX source (`.tex`) is included alongside each PDF for instructors who wish to adapt the materials.

---

## Getting started

```bash
# 1. Clone the repository
git clone git@github.com:storaged/UW-SysBiol-Project2.git
cd UW-SysBiol-Project2

# 2. Create and activate a virtual environment
bash setup_env.sh
source .venv/bin/activate        # macOS / Linux
# .venv\Scripts\activate         # Windows

# 3. Download the data from Moodle and place it in the project root

# 4. Run a quick single-block analysis to verify the setup
python scripts/spatiotemporal_signal_propagation.py \
    --exp-id 1 --site-id 1 \
    --signal-col ERKKTR_ratio \
    --output-dir analysis_outputs
```

---

## Dependencies

See `requirements.txt`. Core dependencies: `numpy`, `pandas`, `scipy`, `matplotlib`.
All notebooks require `jupyter` (or JupyterLab).

---

## Signal glossary

| Symbol | Meaning |
|---|---|
| $S_k(t)$ | Signal value (ERK-KTR or FoxO3A ratio) for cell $k$ at frame $t$ |
| $\Delta S_k(t)$ | Frame-to-frame signal increment |
| $\mathcal{J}_k(t)$ | Jump indicator: large upward step exceeding threshold $\theta$ |
| $\mathcal{N}_k(t)$ | Spatial neighbourhood of cell $k$ at frame $t$ |
| $E_k(t)$ | Exposure indicator: at least one jumping neighbour at frame $t$ |
| $\mathcal{F}_k(t)$ | Future self-jump: cell $k$ jumps within the next $W$ frames |
| $\mathrm{RR}$ | Relative risk: future jump rate (exposed) / future jump rate (unexposed) |
| $\mathrm{RD}$ | Risk difference: absolute gap between the two rates |
