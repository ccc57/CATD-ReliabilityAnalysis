# Test-retest reliability of functional connectivity in depressed adolescents

### important documents

*TODO: Put links here*

## Authors

*TODO: Put list of authors (possibly with links to orcid?*

## Abstract

*TODO: describe the purpose of the repo

## Where/how to get the data

*TODO: get data on openneuro and get link to data

## How to install dependencies

*TODO: how to install dependencies required to run code*

# List of steps to reproduce results

|Filename|Description|Input|Output|
|---------|------------------------|----------------|----------------|
|[`notebooks/run_fmriprepv21.0.0.ipynb`](../master/notebooks/run_fmriprepv21.0.0.ipynb)|Run fmriprep|[bids data](../master/data/CATD)|fmripreped derivatives|
|[fmriprep group report](https://github.com/transatlantic-comppsych/fmriprep-group-report)|Run fmriprep group report command line tool|fmripreped derivatives|
|[`notebooks/add_additional_confounds.ipynb`](../master/notebooks/add_additional_confounds.ipynb)|Add additional confounds to fmriprep files|fmripreped derivatives|fmripreped derivatives|
|[`notebooks/prep_tables.ipynb`](../master/notebooks/prep_tables.ipynb)|Create tables for subsequent processing|fmripreped derivatives|[summary tables](../master/data/CATD/derivatives/summary_tables)|
|[`notebooks/extract_all_timeseries.ipynb`](../master/notebooks/extract_all_timeseries.ipynb)|Extract all the timeseries|fmripreped derivatives, [summary tables](../master/data/CATD/derivatives/summary_tables), [atlas](../master/data/CATD/derivatives/rest_processed/basc122_flat_with_ldlpfc.nii.gz)|[extracted timeseries](../master/data/CATD/derivatives/rest_processed/ts_extract_allmodels)|
