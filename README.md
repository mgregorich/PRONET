
<h1 align="center"> Individual-specific networks for prediction modelling </h1>

## Background

The simulation study examines the use of graph-theoretical features
extracted from individual-specific networks for prognostic prediction
modelling. In particular, the study focuses on individual-specific
network sparsification, which is mostly performed in applied data
analysis by selecting any fixed cut-off value. We propose a flexible
parameterization of the graph-theoretical features generated by a
sequential sparsification and investigate whether it provides an
improvement in terms of the prediction performance of the model for a linear outcome.


<p align="center">
    <img src="./figures/ISN.png" style="width:60%" />
</p>


## Contents

This repository contains the main code, R functions and auxiliary
investigative files to conduct the simulation study and reproduce its
results.

**Folders**

- `main`: contains the code and r functions for the simulation study
- `auxiliary`: contains additional R Markdown reports for examining individual aspects of the simulation study
- `figures`: contains layout graphics

## Usage

The main folder `main` contains all the necessary codes and R functions that are required to carry out the simulation study.

- **x_setup.R** - parameter setup e.g. define number of iterations, sample size, effect sizes, beta-Bernoulli design of networks
- **x_functions.R** - all R functions used in code

- **00_main.R** - main file to carry out the simulation for parameters specified in **setup.R**
- **01_data_generation.R** - data generation of the individual specific networks and the linear outcome
- **02_data_analysis.R** - applies the most common methods of network inference for predictive modelling and our proposed approach to a generated dataset
- **03_summarize_results.R** - summarizes results across iterations
- **04_report_results.Rmd** - R Markdown report containing simulation results and further explanations regarding statistical methods

## Installation

You can install the development version from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("mgregorich/PRONET")
```

## Prerequisites

The code uses the statistical software `R` (>= 4.0)

## Disclaimer

This repository is still in the development phase.
