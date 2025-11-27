# Evolutioanry Trade-offs Between Intergenerational and Transgenerational Fitness Effects

This repository contains the full dataset and analysis scripts associated with the study "Evolutioanry Trade-offs Between Intergenerational and Transgenerational Fitness Effects".  

All analyses were carried out in **R** and the repository provides a complete, reproducible workflow from raw data loading to figure and model output generation.

---

## 📂 Repository structure

```
Multigenerational-Trade-offs/
│
├── /data/                     # Raw CSV files used for analyses
├── /scripts/                  # Analytical scripts for P0, F1 & F3 generations
├── /plots/                    # Exported figures (final publication quality)
├── /model_outputs/            # Model tables exported as .docx files
├── run_all.R                  # Main pipeline to reproduce all analyses
└── README.md                  # Project documentation
```

---

## 🔧 Requirements

Software:
- **R ≥ 4.2**

Recommended IDE:
- **RStudio**

R packages required (automatically loaded in scripts):
```
here, dplyr, tidyr, ggplot2, glmmTMB, DHARMa, emmeans, car,
patchwork, gt, survival, survminer, coxme
```

---

## ▶️ Running the full pipeline

To reproduce the full analysis (stats + plots + exported tables):

1. Clone or download the repository
2. Open the project in **RStudio**
3. Run:

```r
source("run_all.R")
```

When working correctly, you will see:
- Model outputs exported to `/model_outputs/`
- Publication-ready figures exported to `/plots/`
- Console messages indicating pipeline progress

---

## ✨ What the analysis produces

The workflow generates:
| Output type | Location | Details |
|------------|----------|---------|
| Figures | `/plots/` | Age-specific reproduction curves, LRS, Fitness, Survival |
| Model tables | `/model_outputs/` | Mixed models + Cox survival outputs exported as `.docx` |
| data objects | In memory | Used for downstream models & figures |

---

## 👩‍🔬 Authors

**Isaac Harris¹, Elizabeth M. L. Duxbury¹, Tracey Chapman¹, Simone Immler¹, Alexei A. Maklakov¹**

¹ School of Biological Sciences, University of East Anglia,  
Norwich Research Park, Norwich NR4 7TJ, UK

---

## ⚖️ License

This repository is shared for academic transparency and reproducibility.  

---

### 📬 Contact

For questions, issues, or collaboration, please email **Isaac.Harris@uea.ac.uk**.
