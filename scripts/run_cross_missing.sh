#!/bin/bash

#SBATCH --job-name=vasil
#SBATCH --partition=main
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=5GB
#SBATCH --qos=standard
#SBATCH --time=5-00:00:00

python "path_to/scripts/Cross_neutralization/Compute_FR.py" "path_to/Country/results/SpikeGroups.pck" "path_to/Country/results/Mutation_Profiles.pck" "path_to/Country/results/epitope_data/dms_per_ab_per_site.csv" "missing" "path_to/Country/results/Cross_react_dic_spikegroups_ALL.pck" "output_path_to/Country/results/Cross_with_delta.pck" "output_path_to/Country/results/Cross_react_dic_spikegroups_present.pck" "cluster_True" 50
