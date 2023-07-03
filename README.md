# Test-retest reliability of functional connectivity in depressed adolescents

[Access the preprint here](https://www.medrxiv.org/content/10.1101/2022.10.11.22280962v1.article-info)

## Authors

[Chris C. Camp](https://orcid.org/0000-0003-2960-8872), [Stephanie Noble](https://orcid.org/0000-0002-4804-5553), [Dustin Scheinost](https://orcid.org/0000-0002-6301-1167), [Argyris Stringaris](https://orcid.org/0000-0002-6264-8377), [Dylan M. Nielson](https://orcid.org/0000-0003-4613-6643)

## Abstract

This repository contains the code and data used in [Test-retest reliability of functional connectivity in depressed adolescents](https://www.medrxiv.org/content/10.1101/2022.10.11.22280962v1.article-info). This work investigates the reliability of fMRI functional connectivity in a longitudinal dataset of adolescents with and without major depressive disorder. We leverage univariate (intraclass correlation coefficient) and multivariate (fingerprinting and discriminability) reliability metrics to conduct a comprehensive analysis of reliability over a one year period. We did not find strong evidence for an association between depression and reliability; both groups had low univariate reliability and high multivariate reliability.  
Preprocessing was conducted in Python, and most analyses were conducted in R. Please reach out to the corresponding author ([Chris C. Camp](chris.camp@yale.edu)) with any questions. 

## Where/how to get the data

*TODO: get data on openneuro and get link to data

## How to install dependencies

R Code: Install [renv.lock](../master/renv.lock) using [Renv](https://rstudio.github.io/renv/articles/renv.html)

# List of steps to reproduce results

### Preprocessing

|Filename|Description|Input|Output|
|---------|------------------------|----------------|----------------|
|[`notebooks/run_fmriprepv21.0.0.ipynb`](../master/notebooks/run_fmriprepv21.0.0.ipynb)|Run fmriprep|[bids data](../master/data/CATD)|fmripreped derivatives|
|[fmriprep group report](https://github.com/transatlantic-comppsych/fmriprep-group-report)|Run fmriprep group report command line tool|fmripreped derivatives|
|[`notebooks/add_additional_confounds.ipynb`](../master/notebooks/add_additional_confounds.ipynb)|Add additional confounds to fmriprep files|fmripreped derivatives|fmripreped derivatives|
|[`notebooks/prep_tables.ipynb`](../master/notebooks/prep_tables.ipynb)|Create tables for subsequent processing|fmripreped derivatives|[summary tables](../master/data/CATD/derivatives/summary_tables)|
|[`notebooks/extract_all_timeseries.ipynb`](../master/notebooks/extract_all_timeseries.ipynb)|Extract all the timeseries|fmripreped derivatives, [summary tables](../master/data/CATD/derivatives/summary_tables), [atlas](../master/data/CATD/derivatives/rest_processed/basc122_flat_with_ldlpfc.nii.gz)|[extracted timeseries](../master/data/CATD/derivatives/rest_processed/ts_extract_allmodels)|
|[`notebooks/rs_connectivity.ipynb`](../master/notebooks/rs_connectivity.ipynb)|Generate connectivity matrices|[dataframe](../master/references/rest_df_v21_1.csv) with paths to [extracted timeseries](../master/data/CATD/derivatives/rest_processed/ts_extract_allmodels)|[connectivity matrices](../master/data/processed/connectomes),  updated [dataframe](../master/references/rest_df_v21_1.csv)|
|[`notebooks/prep_icc_clean.ipynb`](../master/notebooks/prep_df_clean.ipynb)|Create dataframe of vectorized connectivity matrices|[connectivity matrices](../master/data/processed/connectomes)|[connectivity dataframe](../master/data/connectivity_data/connectivity_data_v21_1.csv)|

### Analysis

|Filename|Description|Input|Output|
|---------|------------------------|----------------|----------------|
|[`icc_bootstrap_main.R`](../master/markdown/icc_bootstrap_main.R)|Creates slurm commands for bootstrapped ICCs using [bootrun script](../master/markdown/scripts/icc_bootrun.R)|[subject dataframe](../master/references/rest_df_v21_1.csv), [connectivity dataframe](../master/data/processed/connectomes)|[slurm command file](../master/slurm/cmds/)|
|[`run_icc_stats.ipynb`](../master/notebooks/`run_icc_stats.ipynb)|Maps ICCs to brain atlas|[all ICC dataframe](../master/data/output/bootstrapped_all_icc_means_v21_1.csv)|[Figure 1](../master/figures)|
|[`markdown/RS_MDD_Reliability.Rmd`](../master/markdown/RS_MDD_Reliability)|Run analyses|[subject dataframe](../master/references/rest_df_v21_1.csv), [connectivity dataframe](../master/data/processed/connectomes), [all ICC dataframe](../master/data/output/bootstrapped_all_icc_means_v21_1.csv), [MDD ICC dataframe](../master/data/output/bootstrapped_mdd_icc_means_v21_1.csv), [HV ICC dataframe](../master/data/output/bootstrapped_hv_icc_means_v21_1.csv), [MDD 4mo ICCs](../master/output/mdd_iccs_v2.csv), [MDD 8mo ICCs](../master/output/mdd_iccs_v3.csv), [subject demographics](../master/references/rest_paths_no_dates_race.csv), [bootstrapped ICC data](data/output/bootstrapped_icc_data.csv), [bootstrapped FI data](../master/data/output/), [edge ranklist](../master/data/output/edge_ranklist.csv), [differential power](../master/data/processed/connectivity_data/DP_edges_all_valmat_v21_1.csv), [group consistency](../master/data/processed/connectivity_data/gr_edges_all_valmat_v21_1.csv)|[fingerprinting session matrices](../master/data/processed), [figures 2-3, SI tables 1-3, SI figures 1-3](../master/figures/)|
|[`discr_edgeranks_main.R`](../master/markdown/discr_edgeranks_main.R)|Creates slurm commands for edge discriminability using [discr_edgeranks script](../master/markdown/scripts/discr_edgeranks.R)|[connectivity dataframe](../master/data/processed/connectomes)|[slurm command file](../master/slurm/cmds/)|
|[ID_scripts](../master/ID_scripts)|Matlab scripts from [Horien et al. 2019](https://www.sciencedirect.com/science/article/pii/S1053811919300886) for fingerprinting analyses|[fingerprinting session matrices](../master/data/processed)|[bootstrapped FI data](../master/data/output/), [differential power](../master/data/processed/connectivity_data/DP_edges_all_valmat_v21_1.csv), [group consistency](../master/data/processed/connectivity_data/gr_edges_all_valmat_v21_1.csv)|

