# Novel-NPtest

## A New Composite Mann–Whitney Test for Two-Sample Survival Comparisons with Right-Censored Data

This repository contains the R code accompanying the paper

**Hussain, A. and Ahmad, T. (2026).**
*A New Composite Mann–Whitney Test for Two-Sample Survival Comparisons with Right-Censored Data.*
Journal of Statistical Computation and Simulation.

## Abstract

Comparing two survival distributions in the presence of right-censored observations remains a fundamental problem in survival analysis. Classical procedures such as the Log-rank and Wilcoxon tests often exhibit substantial power loss when the underlying alternative departs from the scenario for which they are optimal.

This repository implements a novel distribution-free composite Mann–Whitney test that addresses this limitation. The proposed procedure decomposes the observed data into uncensored and censored components and constructs a composite statistic by combining two independent Mann–Whitney statistics. The resulting test automatically adapts to a wide range of alternatives, including proportional hazards, early differences, late differences, and crossing survival curves, without requiring pre-specified weighting schemes or tuning parameters.

Extensive simulation studies presented in the paper demonstrate that the proposed test maintains the nominal Type I error rate while providing robust and competitive power across diverse censoring and hazard configurations.

---

## Paper

Published article:

https://doi.org/10.1080/00949655.2026.2660176

Preprint version:

https://arxiv.org/abs/2510.05353

---

## Repository Contents

```text
Novel-NPtest/
│
├── R/                # Core functions
├── Simulations/      # Monte Carlo studies
├── Examples/         # Illustrative examples
├── Data/             # Data used in examples
├── Figures/          # Reproducible figures
└── README.md
```

---

## Main Features

* Composite Mann–Whitney test for right-censored survival data
* Fully nonparametric and distribution-free framework
* No tuning parameters or weight selection required
* Applicable under proportional and non-proportional hazards
* Robust performance under crossing survival curves
* Monte Carlo simulation framework
* Reproducible examples from the paper

---

## Installation

Clone the repository:

```bash
git clone https://github.com/touqeerahmadunipd/Novel-NPtest.git
```

Install required R packages:

```r
install.packages(c(
  "survival",
  "MASS",
  "mvtnorm",
  "ggplot2"
))
```

---

## Usage

Load the main functions:

```r
source("R/CompositeMWTest.R")
```

Example:

```r
result <- CompositeMWTest(
  time   = time,
  status = status,
  group  = group
)

print(result)
```

---

## Reproducing the Results

All simulation studies reported in the paper can be reproduced using the scripts contained in the `Simulations` directory.

Examples include:

* Type I error evaluation
* Power under proportional hazards
* Power under early hazard differences
* Power under late hazard differences
* Power under crossing hazards
* Sensitivity to censoring proportions

---

## Citation

If you use this software in your research, please cite:

```bibtex
@article{HussainAhmad2026,
  author  = {Abid Hussain and Touqeer Ahmad},
  title   = {A New Composite Mann--Whitney Test for Two-Sample Survival Comparisons with Right-Censored Data},
  journal = {Journal of Statistical Computation and Simulation},
  year    = {2026},
  doi     = {10.1080/00949655.2026.2660176}
}
```

---

## Author

**Touqeer Ahmad**

DSTrain-MSCA Research Fellow

Department of Mathematics

University of Oslo, Norway

GitHub: https://github.com/touqeerahmadunipd

---

## License

MIT License
