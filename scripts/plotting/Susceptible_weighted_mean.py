#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
import pickle
import sys
import pandas
import numpy as np
import pdb

"""Spike groups and frequencies"""
file1 = open(sys.argv[1], "rb") 
SpikeGroups_list = pickle.load(file1)["names"]
SpikeGroups_list = np.array(SpikeGroups_list)
file1.close()
frequency_spk_df = pd.read_csv(sys.argv[2])
"""Threshold variant proportion"""
threshold = float(sys.argv[3])
"""Load total population"""
N_pop = float(sys.argv[4])
# imputing frequencies below threshold and normalization

"""Results directory and timeline of simulations"""
results_dir = sys.argv[5]

phold_df = pd.read_csv(results_dir+"/Immunological_Landscape_ALL/Immunized_SpikeGroup_%s_all_PK.csv"%SpikeGroups_list[0]) # just a place holder
try:
    phold_df.drop(columns = "Unnamed: 0", inplace = True)
except:
    pass

pk_cols = phold_df.columns
t = np.arange(1, len(phold_df['Days'])+1, 1) 

# processing of frequency data
try:
    frequency_spk_df.drop(columns = "Unnamed: 0", inplace = True)
except:
    pass

frequency_spk_df = frequency_spk_df[frequency_spk_df['date'].isin(phold_df['Days'])]
frequency_spk_df = frequency_spk_df.loc[:, frequency_spk_df.columns != 'date']
frequency_spk_df = frequency_spk_df.mask(frequency_spk_df < threshold)
frequency_spk_df = frequency_spk_df.fillna(0)
col_sums = frequency_spk_df.sum(axis = 1).values
frequency_spk_df = frequency_spk_df.divide(col_sums, axis="rows")
frequency_spk_df = frequency_spk_df.fillna(0)

pS_all = np.zeros((len(t)-1, len(SpikeGroups_list[SpikeGroups_list!="Wuhan-Hu-1"]), len(pk_cols[pk_cols!="Day since activation"])))
dprop_all = np.zeros((len(t)-1, len(SpikeGroups_list[SpikeGroups_list!="Wuhan-Hu-1"])))
for x in range(len(SpikeGroups_list)):
    if SpikeGroups_list[x] != "Wuhan-Hu-1":
        EI_df = pd.read_csv(results_dir+"/Immunological_Landscape_ALL/Immunized_SpikeGroup_%s_all_PK.csv"%SpikeGroups_list[x])
        try:
            EI_df.drop(columns = "Unnamed: 0", inplace = True) ## column corresponds to indexes
        except:
            pass

        ei_cols = EI_df.columns
        EI_ranges = EI_df.to_numpy()[:, ei_cols!="Days"].astype(float) # all t_half, t_max simulations on the columns
        
        #Compute timeline prop group (for pseudogroups, this is thee same as loading proportion data)
        Pseudo_Prop = frequency_spk_df["Spike. "+SpikeGroups_list[x]]
                                       
        dprop_all[:, x] = np.diff(np.log(Pseudo_Prop[:len(t)]))
        #Estimated susceptible
        for i in range(EI_ranges.shape[1]):
            S_i = (N_pop - EI_ranges[:, i])
            pS_all[:, x, i] = Pseudo_Prop[:len(t)-1]*S_i[:len(t)-1] ### weighted susceptible, remove the last value to be comparable the prop change timeline
         
### Mean timecourse over the pseudogroups for all PK        
pS_all_mean = np.sum(pS_all, axis = 1) ### it's already weighted
dprop_mean = np.sum(dprop_all[~np.isnan(dprop_all)])

#### save for future runs
S_mean_df = pd.DataFrame(data = pS_all_mean, index = phold_df['Days'][:len(t)-1], columns = pk_cols[pk_cols!="Day since activation"])
S_mean_df.to_csv(results_dir+"/Susceptible_weighted_mean_over_spikegroups_all_PK.csv")

dprop_df = pd.DataFrame(data = dprop_mean, index = phold_df['Days'][:len(t)-1], columns = ["Mean dProp"])
dprop_df.to_csv(results_dir+"/mean_proportion_change_over_pseudogroups.csv")
