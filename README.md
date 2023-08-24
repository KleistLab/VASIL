

## System requirements 


### Operating systems

### Prerequisites

...


## SARS-CoV-2 genomic data 
To extract SARS-CoV-2 genomic data from GISAID or RKI, the tool covsonar can be used:
https://github.com/rki-mf1/covsonar

The usage of covsonar after installing and building the database is the following:
python3 sonar.py match --db database.db --date 2021-07-01:2023-04-16 --collection DESH_STICHPROBE RKI_STICHPROBE --tsv > covsonardata.tsv
