#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import numpy.ma as ma
import pickle
import re
import sys
import pdb
from datetime import datetime, timedelta

""" Loads lineage frequency data aready matched with timeline of infection """
frequency_lineage_df = pd.read_csv(sys.argv[1])
try:
	frequency_lineage_df.drop(columns = "Unnamed: 0", inplace = True)
except:
	pass

fq_cols = frequency_lineage_df.columns.astype(str)
days_prop = frequency_lineage_df["date"] ### already matched with timeline of infection -- including excess

""" Remove NONE and UNASSIGNED and nans"""
"""
if "NONE" in fq_cols:
    frequency_lineage_df.drop(columns = "NONE", inplace = True)
    fq_cols = frequency_lineage_df.columns.astype(str)
if "UNASSIGNED" in fq_cols:
    frequency_lineage_df.drop(columns = "UNASSIGNED", inplace = True)
    fq_cols = frequency_lineage_df.columns.astype(str)
if "nan" in fq_cols:
    frequency_lineage_df.drop(columns = "nan", inplace = True)
    fq_cols = frequency_lineage_df.columns.astype(str)
"""

if "week_num" in list(fq_cols):
    unique_lineage = fq_cols[(fq_cols != "date")&(fq_cols != "week_num")]
    weeks_all = frequency_lineage_df["week_num"].astype(str)
    weeks = np.unique(frequency_lineage_df["week_num"].astype(str))
    frequency_lineage = frequency_lineage_df.to_numpy()[:, (fq_cols != "date")&(fq_cols != "week_num")].T.astype(float)
else:
    unique_lineage = fq_cols[fq_cols != "date"]
    weeks = None
    frequency_lineage = frequency_lineage_df.to_numpy()[:, fq_cols != "date"].T.astype(float)


"""Load spikegroups and mutation data  """
variant_mut_data = pd.read_csv(sys.argv[2])
variant_x_name_orig = np.array(variant_mut_data["lineage"].values).astype(str)
mut_x_sites_orig = np.array(variant_mut_data["RBD_NTD_mutations"].values).astype(str)
mut_x_sites_dic = {}
AA_change_dic = {}
unique_muts = np.unique(mut_x_sites_orig)
for i in range(len(variant_x_name_orig)):
    x = variant_x_name_orig[i]
    mut_x = mut_x_sites_orig[i]
    split_mut = mut_x.split("/")
    aa_x = {}
    pos_list = []
    check = 0
    for mut in split_mut:
        pos0 = re.findall(r'\d+', mut)
        if len(pos0) == 1:
            pos = str(pos0[0])
            if pos not in list(aa_x.keys()):
                aa_x[pos] = [mut]
                pos_list.append(pos)
            else:
                aa_x[pos].append(mut)
                check += 1
        
    mut_x_sites_dic[x] = pos_list
    AA_change_dic[x] = aa_x
    

NormProp = np.sum(frequency_lineage, axis = 0)
prop_rounded = np.round(frequency_lineage,decimals = 10)
proportion_lineage = np.divide(prop_rounded, NormProp, out = np.zeros(prop_rounded.shape), where = NormProp != 0)
"""Finalizing variant proportion parameter and aligning with case ascertainement timeline if parameter is given"""
"""Load Infection Data"""
Population_Data = pd.read_csv(sys.argv[9])
days_incidence = list(Population_Data['date'])
if days_incidence[len(days_incidence) - 1] < days_prop[len(days_prop) - 1]:
    sdate = datetime.strptime(days_incidence[len(days_incidence) - 1], "%Y-%m-%d")
    edate = datetime.strptime(days_prop[len(days_prop) - 1], "%Y-%m-%d")
    date_list = list(days_incidence)[:-1]+pd.date_range(sdate,edate-timedelta(days=1),freq='d').strftime('%Y-%m-%d').tolist()
    if days_prop[len(days_prop) - 1] not in date_list:
        date_list.append(days_prop[len(days_prop) - 1])
else:
    date_list = list(days_incidence)

N_days = len(date_list)

Lineages_list = list(unique_lineage)
variant_proportion_orig = np.zeros((len(Lineages_list), N_days))
mask_missing = np.zeros((len(Lineages_list), N_days)).astype(bool)
missing_var = []
for x in range(len(Lineages_list)):
    variant = Lineages_list[x]
    x_lin = Lineages_list.index(variant)
    
    if variant not in list(variant_x_name_orig) and np.any(proportion_lineage[x_lin, :]!=0.):
        missing_var.append(variant)
        
    for k in range(len(date_list)):
        date = date_list[k]
        if date in list(days_prop):
            ik = list(days_prop).index(date)
            variant_proportion_orig[x, k] = proportion_lineage[x_lin, ik]
        else:
            mask_missing[x, k] = True

if len(missing_var) != 0:
    DiffVar = list(set(variant_x_name_orig).symmetric_difference(set(Lineages_list)))
    miss_freq = [DiffVar[i] for i in range(len(DiffVar)) if DiffVar[i] not in missing_var]
    print("Mutation_profile missing/ignored for variants (correct mutation_profile.R if they are included the data time horizon)", missing_var, "\n Frequency data missing for variants", miss_freq)

NormProp = np.sum(variant_proportion_orig, axis = 0)
prop_rounded = np.round(variant_proportion_orig, decimals = 10)
proportion_lineage = np.divide(prop_rounded, NormProp, out = np.zeros(prop_rounded.shape), where = NormProp != 0)

""" Merge proportions by Spike Groups"""
SpikeGroups_list = []
Pseudogroups = pd.read_csv(sys.argv[3])
variant_x_pseudo = np.array(Pseudogroups["lineage"].values).astype(str)
pseudo_members = np.array(Pseudogroups["group"].values).astype(str)

unique_group = np.unique(pseudo_members)
variant_proportion = []
SpikeGroups_dic  = {}
SpikeGroups_dic["Wuhan-Hu-1"] = "Wuhan-Hu-1" ### place holder for wild type

"""Filter out if indicated"""
try:
    filt_params = int(sys.argv[7])/100. ### provided in percentage
except:
    filt_params = 0

Lin_dic = {}
for x in range(len(unique_group)):
    if "/" not in unique_group[x]:
        #if unique_group[x] != "nan":
        if unique_group[x] in Lineages_list:
            where_x = list(Lineages_list).index(unique_group[x])
            prop_x = variant_proportion_orig[where_x, :]
            if np.max(prop_x) > filt_params:
                variant_proportion.append(prop_x)
                SpikeGroups_list.append(unique_group[x])
                SpikeGroups_dic[unique_group[x]] = unique_group[x]
                Lin_dic[unique_group[x]] = prop_x

    else:
        splited_var = np.array(unique_group[x].split("/"))
        "set empty lineage information as str(NaN) = nan because python treat empty entries as NaN and we transform then into strings"
        splited_var[splited_var==""] = "nan"
        where_x = []            
        gName = ""
        for var_x in splited_var:
            if var_x in Lineages_list:
                where_x.append(list(Lineages_list).index(var_x))
                if gName == "" and var_x != "nan":
                    gName = var_x
                    
        
        if len(where_x)!=0:
            prop_x = np.sum(variant_proportion_orig[np.array(where_x), :], axis = 0)
            if np.max(prop_x) > filt_params:
                variant_proportion.append(prop_x)
                SpikeGroups_list.append(gName)
                for s in range(len(where_x)):
                    Om = Lineages_list[where_x[s]]
                    Lin_dic[Om] = variant_proportion_orig[where_x[s], :]
                    SpikeGroups_dic[Om] = gName
    
variant_proportion = np.array(variant_proportion)
NormProp = np.sum(variant_proportion, axis = 0)
prop_rounded = np.round(variant_proportion,decimals = 10)
variant_proportion = np.divide(prop_rounded, NormProp, out = np.zeros(prop_rounded.shape), where = NormProp != 0)

### Re-normalize lineage frequencies
Lin_df = pd.DataFrame(Lin_dic)
col_sums = Lin_df.sum(axis = 1).values
Lin_df = 100*Lin_df.divide(col_sums, axis="rows")
Lin_df["date"] = date_list
Lin_df = Lin_df.fillna(0)
Lin_df.drop(index = Lin_df.index[np.array(Lin_df["date"]) == 0], inplace= True)
try:
   save_lin = str(sys.argv[8])
   if save_lin in ("True", "TRUE", "1"):
       Lin_df.to_csv("results/Daily_Lineages_Freq_%s_percent.csv"%str(int(filt_params*100)))
except:
    if "week_num" in list(fq_cols):
        frequency_lineage_df.drop(columns = "week_num", inplace = True)
    
### lineage frequency is always saved spikegroups_membership.pck  
SpikeGroups_dic["Frequencies"] = Lin_df

""" Save frequency pseudogroup data """
freq_dic = {}
freq_dic["date"] = date_list
for x in range(len(SpikeGroups_list)):
    freq_dic["Spike. "+SpikeGroups_list[x]] = variant_proportion[x, :]*100 

if "Wuhan-Hu-1" not in SpikeGroups_list:
    freq_dic["Wuhan-Hu-1"] = np.zeros(variant_proportion.shape[1])
    SpikeGroups_list.append("Wuhan-Hu-1")
    print("Python Output: Number of Spikegroups: %d (above %.1f %% in some calendar day) + 1 Wuhan-Hu-1"%(len(SpikeGroups_list) - 1, filt_params*100))
else:
    print("Python Output: Number of Spikegroups: %d (above %.1f %% in some calendar day) + 1 Wuhan-Hu-1"%(len(SpikeGroups_list) - 2, filt_params*100))
        

freq_df = pd.DataFrame(freq_dic, index = np.arange(0, len(date_list)))
freq_df.to_csv(sys.argv[4])


print("Python Output: Number of lineages composing Spikegroups (above %.1f %% in some calendar day): %d"%(filt_params*100, len(list(SpikeGroups_dic.keys())) - 1))
    
### Save SpikeGroups_list and Mutation_Profiles
spk_file = open(sys.argv[5], "wb")
pickle.dump({"names":SpikeGroups_list}, spk_file)
spk_file.close()
mut_file = open(sys.argv[6], "wb")
pickle.dump({"positions": mut_x_sites_dic, "AA_changes":AA_change_dic}, mut_file)
mut_file.close()

### Save SpikeGroups_dic also containing lineage frequencies for easier data handling
file = open("Spikegroups_membership.pck", "wb")
pickle.dump(SpikeGroups_dic,file)
file.close()







    
