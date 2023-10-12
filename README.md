# VA*riant-resolved* S*ars-cov-2* I*mmunological* L*andscape* 
[![GPLv3 License](https://img.shields.io/badge/License-GPL%20v3-yellow.svg)](https://opensource.org/licenses/)
![](https://img.shields.io/github/v/release/kleistlab/ginpipe)
[![DOI](https://zenodo.org/badge/681107366.svg)](https://zenodo.org/badge/latestdoi/681107366)
## System requirements 

### Operating System

This workflow was tested on macOS Ventura Version 13.3.

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
scikit-learn.

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

## Output
The main pipeline (`config.yaml`) creates a folder *results*, containing all (intermediate) output, with the following structure:

```
|-- results
 	|-- Cross_react_dic_spikegroups_ALL.pck	# Pairwise cross reactivity between spikegroups (all_il: TRUE or lineage_focus: "ALL")
	|-- Cross_react_dic_spikegroups_*.pck   # Cross reactivity between lineage_focuss and spikegroups
	|-- Cross_with_delta_validation.pck	# Cross reactivity between Delta variant and Wild Type (used of VE fitting)
	|-- Daily_Lineages_Freq.csv		# Daily frequency of specific lineages (in percentage)
	|-- Daily_SpikeGroups_Freqs.csv 	# Daily frequency of spikegroups (in percentage)
	|-- Fold_Resistance_DMS_Sites_Epitopes.csv. # Fold resistance of all the sites present in the DMS data (FR_DMS_sites: TRUE)
	|-- Mutation_Profile.pck	# Mutation profile of each spikegroups 
	|-- SpikeGroups.pck		# List of all spikegroups names
	|-- Immunological_Landscape
		|-- Immunized_SpikeGroup_*_all_PK.csv # Expected number of immunized for lineage_focuss with all combinations of PK parameters
		|-- Susceptible_SpikeGroup_*_all_PK   # Expected number of Susceptible for lineage_focuss with all combinations of PK parameters
		|-- simulation_status_*.csv  	      # Writes if Immunological Landscape "done" or "Error"
	|-- Immunological_Landscape_ALL
		|-- Immunized_SpikeGroup_*_all_PK.csv # Expected number of immunized for all spikegroups with all combinations of PK parameters
		|-- Susceptible_SpikeGroup_*_all_PK   # Expected number of Susceptible for all spikegroupswith all combinations of PK parameters
		|-- simulation_status_*.csv	      # Write if Immunological Landscape "done" or "Error"
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
		|-- absolute_estimate.pdf 	# Absolute growth of the epidemics
		|-- absolute_estimate.svg
	|-- relative
		|-- plot_status.csv	  	# Writes if plot is "done" or "Error"
		|-- relative_fitness_*.pdf	# relative growth advantage of lineage_focuss
		|-- relative_fitness_*.svg
	|-- Cross_spikegroups
		|-- Cross_React_AB_*.pdf	# Heatmap FR for each epitope classes (max of 10 spikegroups)
		|-- Corss_React_AB_*.svg		
	|--- FR_sites
		|--- foldresistance_DMS_sites_epitopes.pdf # Heatmap Fold resistance of all the sites present in the DMS data (FR_DMS_sites: TRUE)
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


Deactivate the environment to exit the pipeline
```
conda deactivate
```
It should take less than a minute to run the pipeline.
The result folder is created in the [`demo`](./demo) folder where you find the output files, as described above. The relative fitness plot of the demo sample should look like this:

![alt text](https://github.com/KleistLab/VASIL/blob/main/demo/plots/relative/relative_fitness_lineage_XXX.svg)

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







