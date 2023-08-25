

## System requirements 


### Operating systems

### Prerequisites
#### Python

version 3.11.3 

#### Install Conda/Miniconda
Conda will manage the dependencies of our pipeline. Instructions can be found here:
[https://docs.conda.io/projects/conda/en/latest/user-guide/install](https://docs.conda.io/projects/conda/en/latest/user-guide/install)


#### R 
version 4.2.3 (2023-03-15)

Packages:

stringr
reshape2
gplots
RColorBrewer
readr
pheatmap

#### Collect required data from other pipelines

##### SARS-CoV-2 genomic data 
To extract SARS-CoV-2 genomic data from GISAID or RKI, the tool covsonar can be used:
https://github.com/rki-mf1/covsonar

The usage of covsonar after installing and building the database is the following:

```
python3 sonar.py match --db database.db --date 2021-07-01:2023-04-16 --collection DESH_STICHPROBE RKI_STICHPROBE --tsv > covsonardata.tsv
```

##### Escape data
Please download the escape data as provided from the Bloom lab:
https://github.com/jbloomlab/SARS2_RBD_Ab_escape_maps/blob/main/processed_data/escape_data.csv

##### GInPipe incidence data
Please use GInPipe pipeline to generate case ascertainment data
https://github.com/KleistLab/GInPipe/tree/main

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

## Mutation Profile 
To generate a mutation profile for set of variants from covsonar use the script scripts/mutationprofile/generate_mutation_profile.R

```
Rscript generate_mutation_profile.R <covsonar data file> <output directory> <prefix> <mutation threshold>
```

It takes input data as obtained by covsonar (each lineage one tsv file) , filters for mutations with a predefined threshold (e.g. 75%) of prevalence in all samples per lineage,
greps for spike proteins and RBD sites. 


## Epitope Landscape 
To generate the epitope landscape use the script scripts/epitopelandscape/epitope_landscape.R

```
Rscript epitope_landscape.R <escape data file> <antibody mapping file> <fold resistance file> <output directory>
```

It uses escape values from escape_data.csv to compute epitope landscape for all antibodies in a class, as well as across all antibodies in a class. 
For aggregating antibodies per class at each site, their mean values are being used. Classes E1, E2.1 and E2.2 are being merged into E12.


