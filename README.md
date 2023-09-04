

## System requirements 


### Operating systems

### Prerequisites
#### Python

version 3.11.3 

#### Install Conda/Miniconda
Conda will manage the dependencies of our pipeline. Instructions can be found here:
[https://docs.conda.io/projects/conda/en/latest/user-guide/install](https://docs.conda.io/projects/conda/en/latest/user-guide/install)


#### Collect required data from other pipelines

#### SARS-CoV-2 genomic data 
To extract SARS-CoV-2 genomic data from GISAID or RKI, the tool covsonar can be used:
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

stringr
reshape2
gplots
RColorBrewer
readr
pheatmap

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
The pipeline creates a folder *results*, containing all (intermediate) output, with the following structure:

TBA

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







