#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import pickle
import re
import sys
import pdb
#import pdb

""" Loads lineage frequency data aready matched with timeline of infection """
frequency_lineage_df = pd.read_csv(sys.argv[1])
try:
	frequency_lineage_df.drop(columns = "Unnamed: 0", inplace = True)
except:
	pass

fq_cols = frequency_lineage_df.columns.astype(str)
days_prop = frequency_lineage_df["date"] ### already matched with timeline of infection -- including excess

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

"""Finalizing variant proportion parameter """
Lineages_list = list(unique_lineage)
variant_proportion_orig = np.zeros((len(Lineages_list), len(days_prop)))
if weeks is not None:
    weekly_prop = np.zeros((len(Lineages_list), len(weeks)))
#missing_var_prop = {}
for x in range(len(Lineages_list)):
    variant = Lineages_list[x]
    x_lin = list(unique_lineage).index(variant)
    variant_proportion_orig[x, :] = proportion_lineage[x_lin, :]
    if weeks is not None:
        for wk in range(len(weeks)):
            locs_wk = np.where(np.array(weeks_all) == weeks[wk])[0]
            w_lin = np.where(np.array(Lineages_list) == variant)[0]
            if np.sum(locs_wk & w_lin)!=0:
                weekly_prop[x, wk] = np.sum(proportion_lineage[x_lin, locs_wk])/np.sum(locs_wk & w_lin) ### average weekly proportions because they sum up to 1 for each day of that week
            else:
                weekly_prop[x, wk] = np.sum(proportion_lineage[x_lin, locs_wk])
    """
    # used in previous versions, now obsolete, spikegroups and mutation profiles are restricted a specific timeline we do not need to care for those that have misssing mutation profiles
    if (variant not in list(mut_x_sites_dic.keys())) & (variant not in Grouped_in_Spike):
        prop_miss = proportion_lineage[x_lin, :]
        missing_var_prop[variant] = prop_miss 
        print("Missing mutation profile for variant %s :\n"%variant, "proportion (min, mean, max):", (np.min(prop_miss), np.mean(prop_miss), np.max(prop_miss)))
    """

### Load filtering threshold
filt = float(sys.argv[7])
if filt !=0:
    ## Filter
    keep_inds = np.any((weekly_prop>=filt), axis = 1)
    Lineages_list = list(np.array(Lineages_list)[keep_inds])
    variant_proportion_orig = variant_proportion_orig[keep_inds, :]
    print("Number of kept lineages (appearing above %.3f in some calendar week): %d"%(filt, len(Lineages_list)))

    
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

check_var = []
for x in range(len(unique_group)):
    if "/" not in unique_group[x]:
        if unique_group[x] in Lineages_list:
            where_x = list(Lineages_list).index(unique_group[x])
            variant_proportion.append(variant_proportion_orig[where_x, :])
            SpikeGroups_list.append(variant_x_pseudo[pseudo_members == unique_group[x]][0])
            SpikeGroups_dic[unique_group[x]] = variant_x_pseudo[pseudo_members == unique_group[x]][0]
            check_var.append(unique_group[x])
    else:
        splited_var = unique_group[x].split("/")
        where_x = []
        for var_x in splited_var:
            if var_x in Lineages_list:
                where_x.append(list(Lineages_list).index(var_x))
                check_var.append(var_x)
        if len(where_x) != 0:
            variant_proportion.append(np.sum(variant_proportion_orig[np.array(where_x), :], axis = 0))
            SpikeGroups_list.append(variant_x_pseudo[pseudo_members == unique_group[x]][0])
        
            for Om in splited_var:
                SpikeGroups_dic[Om] = variant_x_pseudo[pseudo_members == unique_group[x]][0]

variant_proportion = np.array(variant_proportion)

""" Check that all lineages were grouped  (obsolete) """ 
"""
test = False
for x in check_var:
    test = x not in Lineages_list
    if test:
        print("Lineage not included in pseudogroup list", x)
"""

NormProp = np.sum(variant_proportion, axis = 0)
prop_rounded = np.round(variant_proportion,decimals = 10)
variant_proportion = np.divide(prop_rounded, NormProp, out = np.zeros(prop_rounded.shape), where = NormProp != 0)

""" Save frequency pseudogroup data """
freq_dic = {}
freq_dic["date"] = days_prop
for x in range(len(SpikeGroups_list)):
    if SpikeGroups_list[x]!= "Wuhan-Hu-1":
        freq_dic["Spike. "+SpikeGroups_list[x]] = variant_proportion[x, :]*100 


freq_df = pd.DataFrame(freq_dic, index = np.arange(0, len(days_prop)))
freq_df.to_csv(sys.argv[4])

if "Wuhan-Hu-1" not in SpikeGroups_list:
    freq_dic["Wuhan-Hu-1"] = np.zeros(variant_proportion.shape[1])
    SpikeGroups_list.append("Wuhan-Hu-1")
    print("Number of Spikegroups: %d + 1 Wuhan-Hu-1"%(len(SpikeGroups_list) - 1))
else:
    print("Number of Spikegroups: %d + 1 Wuhan-Hu-1"%(len(SpikeGroups_list) - 2))

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







    
