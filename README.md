Scripts to replicate analyses from Durvasula and Price, 2024 "Distinct explanations underlie gene-environment interactions in the UK Biobank." (Preprint)[https://www.medrxiv.org/content/10.1101/2023.09.22.23295969.abstract]

---

# Directory structure

`data/parsed` contains the parsed output of running the analyses (see below)

`simulations/` contains code for the various simulations included in the paper. These can be run standalone in R. Some of the scripts are set up to run on a SLURM cluster for parallelization.

`scipts/0_data_management` code for binning individuals by E variables

`scripts/1_rg` code for running genetic correlation analyses. BOLT-LMM is first used to get GWAS summary statistics and then LDSC is used to get the genetic correlation

`/scripts/2_prs` code for running the PRSxE analyses. Individuals are scored using Plink and then the `reg_update.py` script is used to run the PRSxE regression with covarates

`scripts/3_h2` code for running heritability analyses with BOLT-REML

`scripts/4_parse_results` ipython notebooks for parsing the complicated output of the different programs run above

`scripts/5_plots` scripts for making plots for the paper

# Data availability

The raw, unparsed output of the scripts is available here: Durvasula, A. (2024). Data for "Distinct explanations underlie gene-environment interactions in the UK Biobank." [Data set]. Zenodo. https://doi.org/10.5281/zenodo.13621189

Feel free to email me at first.last@med.usc.edu if you have questions.  