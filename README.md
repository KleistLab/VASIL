# VA*riant-resolved* S*ars-cov-2* I*mmunological* L*andscape* 
[![GPLv3 License](https://img.shields.io/badge/License-GPL%20v3-yellow.svg)](https://opensource.org/licenses/)
![](https://img.shields.io/github/v/release/kleistlab/ginpipe)
[![DOI](https://zenodo.org/badge/681107366.svg)](https://zenodo.org/badge/latestdoi/681107366)
## System requirements 

### Operating System

This workflow was tested on macOS Monterey Version 12.5, macOS Ventura Version 13.3, CentOS Linux 7 (Core), as well as Ubuntu Version 20.04.5 LTS. 
For these OS, the conda explicit spec-files are found in [env/spec-files](https://github.com/KleistLab/VASIL/tree/main/env/spec-files) and to install VASIL, run 

`conda create --name VASIL --file <the spec-file of the OS>` ([see conda doc](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#activating-an-environment)

### Prerequisites
#### Python

version 3.11.3 

Packages:

numpy,
scipy,
openpyxl,
pandas (1.5.3),
matplotlib,
seaborn,
joblib,
regex,
pip,
pyOpenSSL,
patsy,
scikit-learn,
mpl_axes_aligner,

#### Install Conda/Miniconda
Conda will manage the dependencies of our pipeline. Instructions can be found here:
[https://docs.conda.io/projects/conda/en/latest/user-guide/install](https://docs.conda.io/projects/conda/en/latest/user-guide/install)

#### Collect required data from other pipelines

#### SARS-CoV-2 genomic data 
To extract SARS-CoV-2 genomic data from [GISAID](https://gisaid.org/) or [RKI](https://github.com/robert-koch-institut/SARS-CoV-2-Sequenzdaten_aus_Deutschland/), the covsonar tool can be used:
https://github.com/rki-mf1/covsonar

The usage of covsonar after installing and building the database is the following:

```
python3 sonar.py match --db database.db --date 2021-07-01:2023-04-16 --collection DESH_STICHPROBE RKI_STICHPROBE --tsv > covsonardata.tsv
```

#### Escape data
Please download the escape data as provided from the Bloom lab:
[escape_data.csv](https://github.com/jbloomlab/SARS2_RBD_Ab_escape_maps/blob/main/processed_data/escape_data.csv)

#### GInPipe incidence data
Please use GInPipe pipeline to generate case ascertainment data
https://github.com/KleistLab/GInPipe/tree/main

#### Antibody classification
Please download antibody class file required for epitope classification
[antibody_classes.csv](https://github.com/KleistLab/VASIL/blob/main/scripts/epitopelandscape/antibody_classes.csv)

#### Create the working environment

Create a new environment from the given environment config in [env.yml](https://github.com/KleistLab/VASIL/blob/main/env/env.yml)

```
conda env create -f env/env.yml
```

This step may take a few minutes. 

To activate the eviromnent

```
conda activate VASIL
```

#### Install Snakemake
Snakemake is the workflow management system we use. Install it in your activated environment like this:

```
conda install -c conda-forge -c bioconda snakemake
```

NOTE: In case conda is not able to find the packages for snakemake (which was the case for the Linux version), you can install mamba in your environment

```
conda install -c conda-forge mamba
```

and download snakemake with

```
mamba install -c conda-forge -c bioconda snakemake
```

Detailed Snakemake installation instruction using mamba can be found here:
[https://snakemake.readthedocs.io/en/stable/getting_started/installation.html](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)

#### Install R
version 4.2.3 (2023-03-15)

Packages:

stringr,
reshape2,
gplots,
RColorBrewer,
readr,
pheatmap.

To run R routines, R including Rscript needs to be installed for the workflow. If it is not yet, you can install it together with the needed packages in your activated environment with conda or mamba

If your environment is not yet activated, type
```
conda activate VASIL
```

```
conda install -c conda-forge -c bioconda r-base=4.2.3 r-stringr r-reshape2 r-gplots r-RColorBrewer r-readr r-pheatmap
```
Now the VASIL environment also contains all the dependencies in R

## Input
As an input, the pipeline requires the paths to the consonar data, Escape data, and GInPipe case ascertainment data.
These variables are stored in [`config.yaml`](https://github.com/KleistLab/VASIL/blob/main/config.yaml).
For more information about the YAML markup format refer to documentation: [https://yaml.org](https://yaml.org)
### File formatting requirements
`.tsv` files must be tab separated, if not first the shell command 

`sed 's/,/\t/g' oldfile_commasep.tsv > newfile_tabsep.tsv`

Covsonar data pre-processing that

1) order dates and double check for correct country name (column 'zip' in covsonar data), 

2) removes `invalid date entries`, entries that have lineages `UNASSIGNED` or `nan`,

run the following on the terminal:

`Rscript scripts/check_dataset.R <covsonar data file> <country> <output file>`

## Execution

If your environment is not yet activated, type

```
conda activate VASIL
```
Go to the pipeline directory (where the Snakefile named [`VASIL`](https://github.com/KleistLab/VASIL/blob/main/VASIL) is located) and enter the following command to execute the pipeline

```
snakemake --snakefile VASIL --configfile path/to/config.yaml -j -d path/to/workdir
```
With parameter `--configfile` you can give the configuration file, described above. The `-j` parameter determines the number of available CPU cores to use in the pipeline. Optionally you can provide the number of cores, e.g. `-j 4`. With parameter `-d` you can set the work directory, i.e. where the results of the pipeline are written to.

After generating the main results, our manuscripts figures are obtained by running

```
snakemake --snakefile VASILplots --configfile path/to/config_plots.yaml -j -d path/to/workdir

```
### Further figure improvement
The scripts in `scripts/plotting/` folder can be modified for further figure improvements.

## Output
The main pipeline (`config.yaml`) creates a folder *results*, containing all (intermediate) output, with the following structure:

```
|-- Spikegroups_membership.pck # pickle python dictionary with lineages as keys, indicating the spikegroups of each lineages were assigned to
|-- results
 	|-- Cross_react_dic_spikegroups_ALL.pck	# Pairwise cross reactivity between all spikegroups (can be provided as input, see config.yaml)
	|-- Cross_react_dic_spikegroups_present.pck	# Pairwise cross reactivity between all spikegroups imputing a previous cross-file (e.g. from shorter timeframe runs) and only computing cross for missing spikegroups (see config.yaml)
	|-- Cross_react_dic_spikegroups_*.pck   # Cross reactivity between lineage_focuss and spikegroups
	|-- Cross_with_delta_validation.pck	# Cross reactivity between Delta variant and Wild Type (used of VE fitting)
	|-- Cross_to_major_variants.pck	# Cross reactivity for major variants
	|-- Daily_Lineages_Freq.csv		# Daily frequency of all lineages present within the timeframe (in percentage)
	|-- Daily_SpikeGroups_Freq.csv 	# Daily frequency of spikegroups exceeding a  chosen filter_out threshold (in percentage)
	|-- Daily_Lineages_Freq_*_percent.csv 	# Daily frequency lineages composing the spikegroups exceeding a chosen (*) filter_out  (percentage)	
	|-- Fold_Resistance_DMS_Sites_Epitopes.csv. # Fold resistance of all the sites present in the DMS data (FR_DMS_sites: TRUE)
	|-- Mutation_Profile.pck	# Mutation profile of each spikegroups (post-processed from files in mutation_data/)
	|-- SpikeGroups.pck		# List of all spikegroups names ((post-processed from files in mutation_data/))
	|-- PK_for_all_Epitopes.csv    # Pharmacokynetics for any epitope class with all combinations of PK parameters
	|-- Immunological_Landscape # for lineage_focus (if requested)
		|-- Immunized_SpikeGroup_*_all_PK.csv # Expected number of immunized for lineage_focuss with all combinations of PK parameters
		|-- Susceptible_SpikeGroup_*_all_PK   # Expected number of Susceptible for lineage_focuss with all combinations of PK parameters
		|-- simulation_status_*.csv  	      # Writes if Immunological Landscape "done" or "Error"
		|-- P_neut_*.csv  	      # Virus neutralization probability against selected antigen present in covsonar data (if requested)
	|-- Immunological_Landscape_ALL # for all spikesgroups extracted from covsonar data (if requested)
		|-- Immunized_SpikeGroup_*_all_PK.csv # Expected number of immunized with all combinations of PK parameters
		|-- Susceptible_SpikeGroup_*_all_PK   # Expected number of Susceptible with all combinations of PK parameters
		|-- simulation_status_*.csv	      # Write if Immunological Landscape "done" or "Error"
		|-- P_neut_*.csv  	      # Virus neutralization probability against selected antigen present in covsonar data (if requested) 
	|-- Immunological_Landscape_group # for compare_groups (if requested)
		|-- Immunized_SpikeGroup_*_all_PK.csv # Expected number of immunized with all combinations of PK parameters
		|-- Susceptible_SpikeGroup_*_all_PK   # Expected number of Susceptible with combinations of PK parameters
		|-- simulation_status_*.csv	      # Write if Immunological Landscape "done" or "Error"
                 |-- P_neut_*.csv  	      # Virus neutralization probability against selected antigen present in covsonar data (if requested)
	|-- mutation_data # Cross reactivity files generated for each variants in compare_groups
		|-- cross_status_*.csv	      # Write if cross "done" 
		|-- Cross_*.csv	      # Cross react files	
	|-- mutation_data
		|-- mutationprofile_mutations_spike_lists # Full mutation profile
		|-- mutationprofile_mutations_spike	  # Mutation status for each spike
		|-- mutationprofile_mutations_spikenumber_of_genomes_per_lineage.csv # Number of spike mutations per lineages
		|-- mutationprofile_positiongroups_RBD_NTD_groups_zoom.pdf # Heatmap mutation profile (Zoomed)
		|-- mutationprofile_positiongroups_RBD_NTD_groups.pdf	   # Heatmap mutation profile 
		|-- mutationprofile_RBD_NTD_mutations.csv	# RBD-NTD mutation profile
		|-- mutationprofile_RBD_NTD_pseudogroups.csv    # Spikegroups and their members
```

The figure pipeline (`config_plots.yaml`) add new data to *results* folder and creates a folder *plots*, containing important figures, with the following structure

```
|-- results
 	|-- mean_proportion_changes_over_pseudogroups.csv  	  # Mean change in daily spikegroups proportions 
	|-- Susceptible_weighted_mean_over_spikegroups_all_PK.csv # Mean growth advantage accross  all spikegroups > 1% 
|-- plots	
	|-- absolute
		|-- absolute_estimate.pdf/svg 	# Absolute growth of the epidemics
	|-- relative           # Relative fitness for lineage focus
		|-- plot_status.csv	  	# Writes if plot is "done" or "Error"
		|-- relative_fitness_*.pdf/svg # relative growth advantage of lineage_focus
	|-- relative_groups      # Relative fitness for compare_groups
		|-- As_Spikegroups/relative_fitness_group.pdf/svg # plot relative fitness for chosen variants treating them as spikegroups (compare_groups)
		|-- As_Lineages/relative_fitness_group.pdf/svg # plot relative fitness for chosen variants treating them as Lineages (compare_groups)
        |-- relative_all # relative fitness for all spikegroups
		|-- plot_status_all.csv	  	# Writes if plot is "done" or "Error"
		|-- relative_fitness_*.pdf/svg	# relative growth advantage for all spikegroups present in data
	|-- Cross_spikegroups # Heatmap FR for each epitope classes (max of the 10 first spikegroups)
		|-- Cross_React_AB_*.pdf/svg	
	|-- Cross_Major # Heatmap FR for each epitope classes (for major variants of your choice within spikegroups or with provided mutation data file)
		|-- major_Cross_React_AB_*.svg 		
	|--- FR_sites # Heatmap Fold resistance of all the sites present in the DMS data (FR_DMS_sites: TRUE)
		|--- foldresistance_DMS_sites_epitopes.pdf 
	|--- P_neut_PK_lineage_focus
		|--- PK_Epitope_ranges.pdf/svg # Pharmacokynetics of for any epitope class
		|--- P_Neut_*.pdf/svg # Virus neutralization probability against specific antigen (if requested)
	|--- P_neut_PK_groups
		|--- PK_Epitope_ranges.pdf/svg# Pharmacokynetics of for any epitope class
		|--- P_Neut_*.pdf/svg # Virus neutralization probability against specific antigen (if requested)
```
## Demo
Demo datasets are provided in the repository folder [`demo`](https://github.com/KleistLab/VASIL/tree/main/demo)

If your environment is not yet activated, type

```
conda activate VASIL
```

To run the pipeline go into the repository where the snakefile [`VASIL`](https://github.com/KleistLab/VASIL/blob/main/VASIL) is located and run

```
snakemake --snakefile VASIL --configfile demo/demo_config.yaml -j -d demo

```

After generating the main results, our manuscripts figures are obtained by running

```
snakemake --snakefile VASILplots --configfile demo/demo_config_plots.yaml -j -d demo

```
### Manuscript Figures
The folder `MS_data` can be used to reproduce or find the non-conceptual figures our manuscript figure where only plotting should be done using the config file `ms_fig.yaml`:

```
snakemake --snakefile VASILplots --configfile MS_data/ms_fig.yaml -j -d MS_data

```
The figures are located as follows

```
|-- Figure 1: c: MS_data/results/mutation_data/mutationprofile_positiongroups_RBD_NTD_groups.pdf
              b: MS_data/results/relative_groups/Groups_proportions.pdf,svg
              c: GinPipe
|-- Figure 2: a: MS_data/plots/FR_sites/foldresistance_DMS_sites_epitopes.pdf
              b: MS_data/plots/Cross_major/major_Cross_React_AB_*.pdf,svg
|-- Figure 3: a: MS_data/plots/P_neut_*/PK_Epitopes_ranges.pdf,svg
              b: generated off pipeline
              c-d: MS_data/plots/P_neut_*/P_Neut_.pdf,svg*
|-- Figure 4: a: MS_data/plots/P_neut_groups/P_Neut_.pdf/svg*         
              b: MS_data/plots/relative or relative_all/relative_fitness_*.pdf,svg 
              c: MS_data/plots/relative_groups_germany/As_Spikegroups/relative_fitness_groups.pdf,svg (special file generated from scripts/plotting/relative_advantage_spikegroups_germany_special.py)         
|-- Figure 5: a/b: run command line: python scripts/plotting/compare_countries.py "path/to/Germany/new/path/to/USA" "Germany/USA" "GER/USA" 0.01 "2022-04-09/2022-04-26" "2022-08-21/2023-04-15" 2 path/to/folder_result "BA.2.12.1" "BE.10.ALL" "red" "orange"                          
              c: MS_data/plots/absolute/absolute_estimate.pdf,svg
              d: MS_data/plots/absolute/derivative_trends.pdf,svg
```

### Expected Runtime for demo
Less than 1 min

Deactivate the environment to exit the pipeline
```
conda deactivate
```
The result folder is created in the [`demo`](./demo) folder where you find the output files, as described above. 

The plot for the cross reactivity to chosen major variants should look like this:
![alt text](https://github.com/KleistLab/VASIL/blob/main/MS_data/plots/Cross_Major/major_Cross_React_AB_A.svg)

The plot for the estimation of change in absolute fitness should look like this:
![alt text](https://github.com/KleistLab/VASIL/blob/main/MS_data/plots/absolute/absolute_estimate.svg)

The plot for the estimation of change in relative fitness should look like this:
![alt text](https://github.com/KleistLab/VASIL/blob/main/MS_data/plots/relative_all/relative_fitness_Spike.%20BA.4.svg)

The plot for the estimation of change in relative fitness for compared groups should look like this:
![alt text](https://github.com/KleistLab/VASIL/blob/main/MS_data/plots/relative_groups/As_Spikegroups/relative_fitness_groups.svg)

## Caution
Caution must be taken for all re-parameterization of simulations made with `config.yaml`, snakemake does not execute the rules for which the result files are already present (unless an input file is updated by another rule), remove older files from the *results* folder when needed.

## Resolving some package issues
Some package-related issues might still arise during code excecution, however, most solutions to these issues can be found online. For example here are some issue we encountered

### Issue 1
Error message about parameter issues in snakemake file. This might be a snakefile formating issue, which can be solved by

First install [snakefmt](https://github.com/snakemake/snakefmt) into the `VASIL` enviroment
```
pip install snakefmt
```
Then, when needed, reformat snakefile

```
snakefmt VASIL
```
In case you had to interrupt snakemake run code (e.g. by Ctr + Z), you need to remove the folder `workdir/.snakemake/locks/`

```
rm -rf workdir/.snakemake/locks
```

### Issue 2
```
Importing the numpy C-extensions failed. This error can happen for
many reasons, often due to issues with your setup or how NumPy was
installed
```

Typing the following sequence of code solves this issue [(see stackoverflow)](https://stackoverflow.com/questions/58868528/importing-the-numpy-c-extensions-failed)
```
pip uninstall -y numpy
pip uninstall -y setuptools
pip install setuptools
pip install numpy
```

### Issue 3

```
File ... from scipy.optimize import root ... .../usr/lib/liblapack.3.dylib (no such file)   
```
This is a problem with scipy that is resolved by unistalling and re-installing scipy with `pip` [saturncloud](https://saturncloud.io/blog/pip-installation-of-scipyoptimize-or-scipy-correctly/)

```
pip unistall scipy
```
```
pip install scipy

```

### Issue 4

```
AttributeError: module 'lib' has no attribute 'OpenSSL_add_all_algorithms'
```

Solution:
```
pip uninstall pyOpenSSL
pip install pyOpenSSL
```

### Issue 5 (Apple silicon)
```
Library not loaded: @rpath/liblapack.3.dylib
```

Solution:

```
pip install --upgrade --force-reinstall scikit-learn
```

### Issue 6 
If the conda installation fails, please use the following commands to install it manually:

```
conda create --name VASIL
conda activate VASIL
conda install -c conda-forge -c bioconda -c anaconda python==3.11.3 numpy scipy openpyxl pandas==1.5.3 matplotlib seaborn joblib regex pip pyOpenSSL patsy scikit-learn
```

```
pip install mpl_axes_aligner
```

If mpl_axes_aligner still cannot be found, please install it by:

```
conda install -c d-ice mpl_axes_aligner
```




