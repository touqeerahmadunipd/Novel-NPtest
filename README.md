# Novel-NPtest

## A New Composite Mann–Whitney Test for Two-Sample Survival Comparisons with Right-Censored Data

This repository contains the R code accompanying the paper

**Hussain, A. and Ahmad, T. (2026).**
*A New Composite Mann–Whitney Test for Two-Sample Survival Comparisons with Right-Censored Data.*
Journal of Statistical Computation and Simulation.

## Abstract

A fundamental challenge in comparing two survival distributions with right-censored data is the selection of an appropriate nonparametric test, as the power of standard tests like the Log-rank and Wilcoxon is highly dependent on the often-unknown nature of the alternative hypothesis. This paper introduces a new, distribution-free two-sample test designed to overcome this limitation. The proposed method is based on a strategic decomposition of the data into uncensored and censored subsets, from which a composite test statistic is constructed as the sum of two independent Wilcoxon statistics. This design allows the test to automatically and inherently adapt to various patterns of difference – including early, late, and crossing hazards – without requiring pre-specified parameters, pre-testing, or complex weighting schemes. An extensive Monte Carlo simulation study demonstrates that the proposed test robustly maintains the nominal Type I error rate. Crucially, its power is highly competitive with the optimal traditional tests in standard scenarios and superior in complex settings with crossing survival curves, while also exhibiting remarkable robustness to high levels of censoring. The test’s power effectively approximates the maximum power achievable by either the Log-rank or Wilcoxon tests across a wide range of alternatives, offering a powerful, versatile, and computationally simple tool for survival analysis.

---

## Paper

Published article:

https://doi.org/10.1080/00949655.2026.2660176

Preprint version:

https://arxiv.org/abs/2510.05353

---

---

## Reproducing the Results

All simulation studies reported in the paper can be reproduced.

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
@article{hussain2026new,
  title={A new composite Mann-Whitney test for two-sample survival comparisons with right-censored data},
  author={Hussain, Abid and Ahmad, Touqeer},
  journal={Journal of Statistical Computation and Simulation},
  pages={1--21},
  year={2026},
  publisher={Taylor \& Francis},
  doi     = {10.1080/00949655.2026.2660176}
}
```

## License

MIT License
