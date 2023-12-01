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
#import pdb

""" Loads lineage frequency data aready matched with timeline of infection """
frequency_lineage_df = pd.read_csv(sys.argv[1])
try:
	frequency_lineage_df.drop(columns = "Unnamed: 0", inplace = True)
except:
	pass

fq_cols = frequency_lineage_df.columns.astype(str)
days_prop = frequency_lineage_df["date"] ### already matched with timeline of infection -- including excess

""" Remove NONE and UNASSIGNED """
if "NONE" in fq_cols:
    frequency_lineage_df.drop(columns = "NONE", inplace = True)
    fq_cols = frequency_lineage_df.columns.astype(str)
if "UNASSIGNED" in fq_cols:
    frequency_lineage_df.drop(columns = "UNASSIGNED", inplace = True)
    fq_cols = frequency_lineage_df.columns.astype(str)

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
variant_x_names = []
unique_muts = np.unique(mut_x_sites_orig)
Grouped_in_Spike = []
for i in range(len(unique_muts)):
    where_mut = mut_x_sites_orig == unique_muts[i]
    mut_x_sites_dic[variant_x_name_orig[where_mut][0]] = re.findall(r'\d+', unique_muts[i])
    variant_x_names.append(variant_x_name_orig[where_mut][0])
    Grouped_in_Spike += list(variant_x_name_orig[where_mut])
    
NormProp = np.sum(frequency_lineage, axis = 0)
prop_rounded = np.round(frequency_lineage,decimals = 10)
proportion_lineage = np.divide(prop_rounded, NormProp, out = np.zeros(prop_rounded.shape), where = NormProp != 0)
"""Finalizing variant proportion parameter and aligning with case ascertainement timeline if parameter is given"""
#try:
"""Load Infection Data"""
Population_Data = pd.read_csv(sys.argv[9])
days_incidence = list(Population_Data['date'])
if days_incidence[len(days_incidence) - 1] < days_prop[len(days_prop) - 1]:
    sdate = datetime.strptime(days_incidence[len(days_incidence) - 1], "%Y-%m-%d")
    edate = datetime.strptime(days_prop[len(days_prop) - 1], "%Y-%m-%d")
    date_list = list(days_incidence)[:-1]+pd.date_range(sdate,edate-timedelta(days=1),freq='d').strftime('%Y-%m-%d').tolist()
else:
    date_list = list(days_incidence)

N_days = len(date_list)

#except:
#    N_days = len(days_prop)

Lineages_list = list(unique_lineage)
variant_proportion_orig = np.zeros((len(Lineages_list), N_days))
mask_missing = np.zeros((len(Lineages_list), N_days)).astype(bool)
#missing_var_prop = {}
for x in range(len(Lineages_list)):
    variant = Lineages_list[x]
    x_lin = list(unique_lineage).index(variant)
    for k in range(len(date_list)):
        date = date_list[k]
        if date in list(days_prop):
            ik = list(days_prop).index(date)
            variant_proportion_orig[x, k] = proportion_lineage[x_lin, ik]
        else:
            mask_missing[x, k] = True
        
variant_proportion_orig = ma.masked_array(variant_proportion_orig, mask = mask_missing)
NormProp = np.sum(variant_proportion_orig, axis = 0)
prop_rounded = np.round(variant_proportion_orig, decimals = 10)
proportion_lineage = np.divide(prop_rounded, NormProp, out = np.zeros(prop_rounded.shape), where = NormProp != 0)

""" Merge proportions by Spike Groups"""
SpikeGroups_list = []
Pseudogroups = pd.read_csv(sys.argv[3])
variant_x_pseudo = np.array(Pseudogroups["lineage"].values).astype(str)
pseudo_members = np.array(Pseudogroups["group"].values).astype(str)
### remove nans
where_not_nans = variant_x_pseudo != "nan"
variant_x_pseudo = variant_x_pseudo[where_not_nans]
pseudo_members = pseudo_members[where_not_nans]

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
        if unique_group[x] != "nan":
            if unique_group[x] in Lineages_list:
                where_x = list(Lineages_list).index(unique_group[x])
                prop_x = variant_proportion_orig[where_x, :]
                if np.max(prop_x) > filt_params:
                    variant_proportion.append(prop_x)
                    SpikeGroups_list.append(variant_x_pseudo[pseudo_members == unique_group[x]][0])
                    SpikeGroups_dic[unique_group[x]] = variant_x_pseudo[pseudo_members == unique_group[x]][0]
                    Lin_dic[unique_group[x]] = prop_x
    else:
        splited_var = unique_group[x].split("/")
        where_x = []
        for var_x in splited_var:
            if var_x in Lineages_list:
                if var_x!="nan":
                    where_x.append(list(Lineages_list).index(var_x))
        
        if len(where_x)!=0:
            prop_x = np.sum(variant_proportion_orig[np.array(where_x), :], axis = 0)
            if np.max(prop_x) > filt_params:
                variant_proportion.append(prop_x)
                SpikeGroups_list.append(variant_x_pseudo[pseudo_members == unique_group[x]][0])
                
                for s in range(len(where_x)):
                    Om = Lineages_list[where_x[s]]
                    Lin_dic[Om] = variant_proportion_orig[where_x[s], :]
                    SpikeGroups_dic[Om] = variant_x_pseudo[pseudo_members == unique_group[x]][0]

variant_proportion = np.array(variant_proportion)
NormProp = np.sum(variant_proportion, axis = 0)
prop_rounded = np.round(variant_proportion,decimals = 10)
variant_proportion = np.divide(prop_rounded, NormProp, out = np.zeros(prop_rounded.shape), where = NormProp != 0)

try:
   save_lin = str(sys.argv[8])
   if save_lin in ("True", "TRUE", "1"):
       Lin_df = pd.DataFrame(Lin_dic)
       col_sums = Lin_df.sum(axis = 1).values
       Lin_df = 100*Lin_df.divide(col_sums, axis="rows")
       Lin_df["date"] = days_prop
       Lin_df.to_csv("results/Daily_Lineages_Freq_%s_percent.csv"%str(int(filt_params*100)))
except:
    pass
    

""" Save frequency pseudogroup data """
freq_dic = {}
freq_dic["date"] = date_list
for x in range(len(SpikeGroups_list)):
    if SpikeGroups_list[x]!= "Wuhan-Hu-1":
        freq_dic["Spike. "+SpikeGroups_list[x]] = variant_proportion[x, :]*100 


freq_df = pd.DataFrame(freq_dic, index = np.arange(0, len(date_list)))
freq_df.to_csv(sys.argv[4])
if "Wuhan-Hu-1" not in SpikeGroups_list:
    freq_dic["Wuhan-Hu-1"] = np.zeros(variant_proportion.shape[1])
    SpikeGroups_list.append("Wuhan-Hu-1")
    print("Number of Spikegroups: %d + 1 Wuhan-Hu-1"%(len(SpikeGroups_list) - 1))
else:
    print("Number of Spikegroups: %d + 1 Wuhan-Hu-1"%(len(SpikeGroups_list) - 2))

print("Number of lineages composing Spikegroups: %d"%len(list(SpikeGroups_dic.keys())))
### Save SpikeGroups_list and Mutation_Profiles
spk_file = open(sys.argv[5], "wb")
pickle.dump({"names":SpikeGroups_list}, spk_file)
spk_file.close()
mut_file = open(sys.argv[6], "wb")
pickle.dump({"positions": mut_x_sites_dic}, mut_file)
mut_file.close()

### Save SpikeGroups_dic
file = open("Spikegroups_membership.pck", "wb")
pickle.dump(SpikeGroups_dic,file)
file.close()







    
