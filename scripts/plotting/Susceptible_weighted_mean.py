#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
import pickle
import sys
import pandas
import numpy as np
import pdb
import numpy.ma as ma  
import os      

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
# processing of frequency data
try:
    frequency_spk_df.drop(columns = "Unnamed: 0", inplace = True)
except:
    pass
    
"""Remove missing spikegroups data """
i = 0
if results_dir not in ("vaccination_special_ver1", "vaccination_special_ver2"):
    results_dir_full = results_dir+"/Immunological_Landscape_ALL"
    for spk in SpikeGroups_list:
        if not os.path.exists(results_dir_full+"/Immunized_SpikeGroup_%s_all_PK.csv"%spk):
            if spk != "Wuhan-Hu-1": 
                frequency_spk_df.drop(columns = "Spike. "+spk , inplace = True)  
            else:
                frequency_spk_df.drop(columns = "Wuhan-Hu-1" , inplace = True)  
        else:
            if i == 0:
                phold_df = pd.read_csv(results_dir_full+"/Immunized_SpikeGroup_%s_all_PK.csv"%spk) # just a place holder
                try:
                    phold_df.drop(columns = "Unnamed: 0", inplace = True)
                except:
                    pass
    
                pk_cols = phold_df.columns
                t = np.arange(1, len(phold_df['Days'])+1, 1) 
            i+=1
else:
    ver = results_dir[-4:]
    results_dir_prev = results_dir
    results_dir= "vaccination" # hard-coded
    results_dir_full = results_dir+"/ImL_ALL_vs_Vacc_"+ver
    for spk in SpikeGroups_list:
        if not os.path.exists(results_dir_full+"/Immunized_SpikeGroup_%s_all_PK.csv"%spk):
            if spk != "Wuhan-Hu-1": 
                frequency_spk_df.drop(columns = "Spike. "+spk , inplace = True)  
            else:
                frequency_spk_df.drop(columns = "Wuhan-Hu-1" , inplace = True)  
        else:
            if i == 0:
                phold_df = pd.read_csv(results_dir_full+"/Immunized_SpikeGroup_%s_all_PK.csv"%spk) # just a place holder
                try:
                    phold_df.drop(columns = "Unnamed: 0", inplace = True)
                except:
                    pass
    
                pk_cols = phold_df.columns
                t = np.arange(1, len(phold_df['Days'])+1, 1) 
            i+=1
    results_dir = results_dir_prev
    
if len(frequency_spk_df.columns) > 1:  
    days_prop = frequency_spk_df['date'][frequency_spk_df['date'].isin(phold_df['Days'])]
    frequency_spk_df = frequency_spk_df[frequency_spk_df['date'].isin(phold_df['Days'])]
    frequency_spk_df = frequency_spk_df.loc[:, frequency_spk_df.columns != 'date']
    frequency_spk_df = frequency_spk_df.mask(frequency_spk_df < threshold)
    frequency_spk_df = frequency_spk_df.fillna(0)
    col_sums = frequency_spk_df.sum(axis = 1).values
    frequency_spk_df = frequency_spk_df.divide(col_sums, axis="rows")
    frequency_spk_df = frequency_spk_df.fillna(0)
    prop_mask = np.all(frequency_spk_df.loc[:, frequency_spk_df.columns != 'date'] == 0, axis = 1)
else:
    sys.exit("All E[Immunized] are missing, first compute them: i.e., set all_il = TRUE in main config.yaml")

pS_all = np.zeros((len(t)-1, len(SpikeGroups_list[SpikeGroups_list!="Wuhan-Hu-1"]), len(pk_cols[pk_cols!="Days"])))
dprop_all = np.zeros((len(t)-1, len(SpikeGroups_list[SpikeGroups_list!="Wuhan-Hu-1"])))
p_prop = np.zeros(pS_all.shape)

for x in range(len(SpikeGroups_list)):
    if SpikeGroups_list[x] != "Wuhan-Hu-1":
        if os.path.exists(results_dir_full+"/Immunized_SpikeGroup_%s_all_PK.csv"%SpikeGroups_list[x]):
            EI_df = pd.read_csv(results_dir_full+"/Immunized_SpikeGroup_%s_all_PK.csv"%SpikeGroups_list[x])
            try:
                EI_df.drop(columns = "Unnamed: 0", inplace = True) ## column corresponds to indexes
            except:
                pass
            
            ei_cols = EI_df.columns
            EI_ranges = EI_df.to_numpy()[:, ei_cols!="Days"].astype(float) # all t_half, t_max simulations on the columns
            days = list(EI_df["Days"])
            
            #Compute timeline prop group (for pseudogroups, this is thee same as loading proportion data)
            Pseudo_Prop = frequency_spk_df["Spike. "+SpikeGroups_list[x]]
            Pseudo_Prop = list(Pseudo_Prop)
            
            gamma_prop = np.zeros(len(t)-1)
            Prop_sub = np.zeros(len(t)-1)
            mask_ind = []
            for l in range(len(t)-1):
                if days[l] in list(days_prop):
                    w_l = list(days_prop).index(days[l])
                    try: # should work if proportions date where proprerly aligned with incidence data
                        if Pseudo_Prop[w_l] == 0 or Pseudo_Prop[w_l+1] == 0:
                            gamma_prop[l] = float('nan')
                        else:
                            gamma_prop[l] = Pseudo_Prop[w_l+1]/Pseudo_Prop[w_l] -1
                    except:
                        gamma_prop[l] = float("nan")
                        
                    mask_ind.append(False)
                    Prop_sub[l] = Pseudo_Prop[w_l]
                else:
                    gamma_prop[l] = float('nan')
                    mask_ind.append(True)
                    
            mask_ind = np.array(mask_ind) + prop_mask[:len(t)-1]
            dprop_all[:, x] = gamma_prop
            Pseudo_Prop_Masked = ma.masked_array(Prop_sub, mask = mask_ind)
            #dprop_all[:, x] = np.diff(np.log(Pseudo_Prop[:len(t)][0]))
            #Estimated susceptible
            for i in range(EI_ranges.shape[1]):
                S_i = (N_pop - EI_ranges[:, i])
                pS_all[:, x, i] = Pseudo_Prop_Masked*S_i[:len(t)-1] ### weighted susceptible, remove the last value to be comparable the prop change timeline
        
        """
        else:
            
            print(SpikeGroups_list[x], "File not found")
            a = np.empty(len(dprop_all[:, x]))
            b = np.empty(pS_all[:, x, :].shape)
            c = np.zeros(p_prop[:, x, :].shape)
            a[:] = np.nan
            b[:] = np.nan
            dprop_all[:, x] = a
            pS_all[:, x, :] = b
            p_prop[:, x, :] = c
        """

### Mean timecourse over the pseudogroups for all 
pS_all =  ma.masked_array(pS_all, mask=np.isnan(pS_all))
pS_all_mean = np.sum(pS_all, axis = 1) ## if already weighted
dprop_mean = np.sum(dprop_all[~np.isnan(dprop_all)])

#### save for future runs
S_mean_df = pd.DataFrame(data = pS_all_mean, index = phold_df['Days'][:len(t)-1], columns = pk_cols[pk_cols!="Days"])
if results_dir not in ("vaccination_special_ver1", "vaccination_special_ver2"):
    S_mean_df.to_csv(results_dir+"/Susceptible_weighted_mean_over_spikegroups_all_PK.csv")
else:
    ver = results_dir[-4:]
    results_dir = "vaccination"
    S_mean_df.to_csv(results_dir+"/Susceptible_weighted_mean_over_spikegroups_vs_Vacc_%s_all_PK.csv"%ver)

dprop_df = pd.DataFrame(data = dprop_mean, index = phold_df['Days'][:len(t)-1], columns = ["Mean dProp"])
dprop_df.to_csv(results_dir+"/mean_proportion_change_over_pseudogroups.csv")
