# Evaluating the Effectiveness of International Travel Controls to Identify Monkeypox Virus Infected Travelers

## Contents of this file

 - Introduction
 - Requirements
 - Usage
 - Data
 - Maintainer

## Introduction

This repository contains the code and data used for the following manuscript: 
Ejima, K., Wang, Y., Endo, A. et al. Evaluating the effectiveness of international travel controls to identify MPXV-infected travelers: a simulation study. BMC Med 23, 473 (2025). https://doi.org/10.1186/s12916-025-04286-6


## Requirements

**Note**: R version 4.3.1.

## Usage

 - The function for simulating within-host viral dynamics, calculating false-negative rates and generating the probability density function of post-entry illness onset is available in mpox_function.R.
 - Plotting script can be found in mpox_figure.R.
 - Result from the Monolix model fitting are located in the Monolix folder.
 - Viral load data are stored in the data folder. 

## Data

The viral load data is sourced from the following publication:
Yang, Y., Niu, S., Shen, C., Yang, L., Song, S., Peng, Y., Xu, Y., Guo, L., Shen, L., Liao, Z., Liu, J., Zhang, S., Cui, Y., Chen, J., Chen, S., Huang, T., Wang, F., Lu, H., & Liu, Y. (2024). Longitudinal viral shedding and antibody response characteristics of men with acute infection of monkeypox virus: a prospective cohort study. Nature communications, 15(1), 4488. https://doi.org/10.1038/s41467-024-48754-8

Brosius, I., Van Dijck, C., Coppens, J., Vandenhove, L., Bangwen, E., Vanroye, F., Verschueren, J., ITM MPOX Study Group, Zange, S., Bugert, J., Michiels, J., Bottieau, E., Soentjens, P., van Griensven, J., Kenyon, C., Ariën, K. K., Van Esbroeck, M., Vercauteren, K., & Liesenborghs, L. (2023). Presymptomatic viral shedding in high-risk mpox contacts: A prospective cohort study. Journal of medical virology, 95(5), e28769. https://doi.org/10.1002/jmv.28769



## Maintainer

Current maintainer:
- Yuqian Wang (yuqian.wang[at]uzh.ch)
