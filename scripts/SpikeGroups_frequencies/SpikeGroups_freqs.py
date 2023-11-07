#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import pickle
import re
import sys
#import pdb

""" Loads lineage frequency data aready matched with timeline of infection """
frequency_lineage_df = pd.read_csv(sys.argv[1])
try:
	frequency_lineage_df.drop(columns = "Unnamed: 0", inplace = True)
except:
	pass

fq_cols = frequency_lineage_df.columns.astype(str)
unique_lineage = fq_cols[fq_cols != "date"]
days_prop = frequency_lineage_df["date"] ### already matched with timeline of infection -- including excess
frequency_lineage = frequency_lineage_df.to_numpy()[:, fq_cols != "date"].T.astype(float)

"""Load lineage data  """
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
Lineages_list = ["Wuhan-Hu-1"] + list(unique_lineage)
variant_proportion_orig = np.zeros((len(Lineages_list), len(days_prop)))
missing_var_prop = {}
for x in range(len(Lineages_list)):
    variant = Lineages_list[x]
    if variant == "Wuhan-Hu-1":
        variant_proportion_orig[x, :] = np.zeros(proportion_lineage.shape[1])
    else:
        x_lin = list(unique_lineage).index(variant)
        variant_proportion_orig[x, :] = proportion_lineage[x_lin, :]
        if (variant not in list(mut_x_sites_dic.keys())) & (variant not in Grouped_in_Spike):
            prop_miss = proportion_lineage[x_lin, :]
            missing_var_prop[variant] = prop_miss 
            print("Missing mutation profile for variant %s :\n"%variant, "proportion (min, mean, max):", (np.min(prop_miss), np.mean(prop_miss), np.max(prop_miss)))
            
""" Merge proportions by Spike Groups"""
SpikeGroups_list = []
Pseudogroups = pd.read_csv(sys.argv[3])
variant_x_pseudo = np.array(Pseudogroups["lineage"].values).astype(str)
pseudo_members = np.array(Pseudogroups["group"].values).astype(str)

unique_group = np.unique(pseudo_members)
variant_proportion = np.zeros((len(unique_group), len(days_prop)))
pseudo_members_dic = {}
SpikeGroups_dic  = {}
SpikeGroups_dic ["Wuhan-Hu-1"] = "Wuhan-Hu-1" ### place holder for wild type

check_var = []
for x in range(len(unique_group)):
    if "/" not in unique_group[x]:
        where_x = list(Lineages_list).index(unique_group[x])
        variant_proportion[x, :] = variant_proportion_orig[where_x, :]
        SpikeGroups_list.append(variant_x_pseudo[pseudo_members == unique_group[x]][0])
        SpikeGroups_dic[unique_group[x]] = variant_x_pseudo[pseudo_members == unique_group[x]][0]
        check_var.append(unique_group[x])
    else:
        splited_var = unique_group[x].split("/")
        where_x = []
        for var_x in splited_var:
            where_x.append(list(Lineages_list).index(var_x))
            check_var.append(var_x)
        variant_proportion[x, :] = np.sum(variant_proportion_orig[np.array(where_x), :], axis = 0)
        SpikeGroups_list.append(variant_x_pseudo[pseudo_members == unique_group[x]][0])
        
        for Om in splited_var:
            SpikeGroups_dic[Om] = variant_x_pseudo[pseudo_members == unique_group[x]][0]
    
    pseudo_members_dic[variant_x_pseudo[pseudo_members == unique_group[x]][0]] = unique_group[x]

""" Check that all lineages were grouped """ 
test = False
for x in check_var:
    test = x not in Lineages_list
    if test:
        print("Lineage not included in pseudogroup list", x)

""" Save frequency pseudogroup data """
freq_dic = {}
freq_dic["date"] = days_prop
for x in range(len(SpikeGroups_list)):
    if SpikeGroups_list[x]!= "Wuhan-Hu-1":
        freq_dic["Spike. "+SpikeGroups_list[x]] = variant_proportion[x, :]*100 


freq_df = pd.DataFrame(freq_dic, index = np.arange(0, len(days_prop)))
freq_df.to_csv(sys.argv[4])

if "Wuhan-Hu-1" not in SpikeGroups_list:
    variant_proportion = np.row_stack((variant_proportion, np.zeros(variant_proportion.shape[1])))
    SpikeGroups_list.append("Wuhan-Hu-1")

### Save SpikeGroups_list and Mutation_Profiles
spk_file = open(sys.argv[5], "wb")
pickle.dump({"names":SpikeGroups_list}, spk_file)
spk_file.close()
mut_file = open(sys.argv[6], "wb")
pickle.dump({"positions": mut_x_sites_dic}, mut_file)
mut_file.close()






    
