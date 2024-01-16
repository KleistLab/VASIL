#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 12:29:17 2023

@author: raharinirina
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
from datetime import datetime, timedelta
import joblib as jb
from functools import partial
from scipy.interpolate import interp1d
import sys
import pdb
import re
import pickle

"""Load covsonar_data"""
try:
	covsonar_data = pd.read_csv(sys.argv[1], sep = "\t")
except:
	covsonar_data = pd.read_csv(sys.argv[1])

try:
    days_prop0 = covsonar_data["date"].values.astype(str)
except:
    sys.exit("date column not found in covsonar data")

lineages_all_0 = covsonar_data["lineage"].values.astype(str)

aa_seq = covsonar_data["aa_profile"].values
invalid_aa_index = [i for i in range(len(aa_seq)) if (str(aa_seq[i]) == "nan" and str(lineages_all_0[i]) == "nan")]
kept_aa_index = np.ones(len(aa_seq)).astype(bool)
if len(invalid_aa_index) > 0:
    kept_aa_index[np.array(invalid_aa_index)] = False
    
"""Simulation timeframe"""
date_start = sys.argv[2]
inds_keep = (days_prop0!="nan")&(kept_aa_index)
days_prop = np.array(days_prop0)[inds_keep]
lineages_all = np.array(lineages_all_0)[inds_keep]

"""Start computing Variant-proportions"""
# iniializing variant proportion for all lineages
unique_days_prop_all = list(np.unique(days_prop))
### make sure that dates are sorted and not nan

unique_days_prop_sub = []
for i in range(len(unique_days_prop_all)):
    try:
        try_this = datetime.strptime(unique_days_prop_all[i], "%Y-%m-%d") 
        unique_days_prop_sub.append(unique_days_prop_all[i])
    except:
        """Keep only well formated dates"""
        keep = days_prop != unique_days_prop_all[i]
        days_prop = days_prop[keep]
        lineages_all = lineages_all[keep]

unique_days_prop_sub.sort(key = lambda date: datetime.strptime(date, "%Y-%m-%d")) 
unique_days_prop = unique_days_prop_sub

if date_start in unique_days_prop:
    where_day = list(unique_days_prop).index(date_start)
    unique_days_prop = unique_days_prop[where_day:]
    
print("Timeline of lineage proportions: %s -- %s"%(unique_days_prop[0], unique_days_prop[-1]))
unique_lineage_timeframe = np.unique(lineages_all[np.array([i for i in range(len(days_prop)) if days_prop[i] in unique_days_prop])])
print("Python Output: Number of lineages: ", len(unique_lineage_timeframe))
try:
    seq_thres = int(sys.argv[4])
except:
    seq_thres = None
    
if seq_thres is not None:

    num_seq = []
    for u in range(len(unique_days_prop)):
        count_u = np.sum(days_prop == unique_days_prop[u])
        num_seq.append(count_u)
    
    print("Python Output: Number of genomes: ", np.sum(num_seq))
    i = 0
    n_i = 0
    days_prop_new = days_prop.copy()
    while i<len(unique_days_prop):
        num_goal = 0
        j = i
        where_days = np.zeros(len(days_prop)).astype(bool)
        while num_goal <= seq_thres and j<len(num_seq):
            num_goal += num_seq[j]
            where_days += days_prop == unique_days_prop[j]
            j +=1
        
        grouped = days_prop[where_days]
        useq = np.unique(grouped)
        chosen_mid = useq[len(useq)//2]
        days_prop_new[where_days]= chosen_mid
        
        i = j
        n_i +=len(grouped)
    
    unique_days_prop_new = list(np.unique(days_prop_new))
    if unique_days_prop_new[0]!=unique_days_prop[0]:
        unique_days_prop_new = [unique_days_prop[0]] + unique_days_prop_new
    if unique_days_prop_new[len(unique_days_prop_new) - 1]!=unique_days_prop[-1]:
        unique_days_prop_new = unique_days_prop_new + [unique_days_prop[-1]]
        
    days_prop = np.array(days_prop_new)
    unique_days_prop = unique_days_prop_new
    unique_days_prop.sort(key = lambda date: datetime.strptime(date, "%Y-%m-%d")) 
    
""" Remove NONE and UNASSIGNED """
"""
if "NONE" in unique_lineage_timeframe:
    unique_lineage_timeframe = np.array(unique_lineage_timeframe)[np.array(unique_lineage_timeframe) != "NONE"]
if "UNASSIGNED" in unique_lineage_timeframe:
    unique_lineage_timeframe = np.array(unique_lineage_timeframe)[np.array(unique_lineage_timeframe) != "UNASSIGNED"]
"""

"""Compute lineage frequencies"""
def sub_func2(x, days_prop, unique_days_prop, lineages_all, unique_lineage):
    lin_locs = lineages_all == unique_lineage_timeframe[x]
    locs = days_prop[lin_locs, np.newaxis] == np.array(unique_days_prop)[np.newaxis, :]
    res = np.sum(locs, axis = 0)
    return res

pfunc = partial(sub_func2, days_prop = days_prop, unique_days_prop = unique_days_prop, lineages_all = lineages_all, unique_lineage = unique_lineage_timeframe)
try:
    res = list(jb.Parallel(n_jobs = -1, backend = "loky")(jb.delayed(pfunc)(x) for x in range(len(unique_lineage_timeframe))))
except:
    try:
        res = list(jb.Parallel(n_jobs = -1, backend = "multiprocessing")(jb.delayed(pfunc)(x) for x in range(len(unique_lineage_timeframe))))
    except: 
        res = list(jb.Parallel(n_jobs = -1, prefer = "threads")(jb.delayed(pfunc)(x) for x in range(len(unique_lineage_timeframe))))

frequency_lineage = np.array(res)

"""
def sub_func(s, x, days_prop, unique_days_prop, lineages_all, unique_lineage):
    res = np.sum((days_prop == unique_days_prop[s]) & (lineages_all == unique_lineage[x]))
    return res

frequency_lineage = np.zeros((len(unique_lineage_timeframe), len(unique_days_prop))) # indexing of t correspondis to timeline of infection days_incidence
for x in range(len(unique_lineage_timeframe)):
    pdb.set_trace()
    lin_locs = lineages_all == unique_lineage_timeframe[x]
    locs = days_prop[lin_locs, np.newaxis] == np.array(unique_days_prop)[np.newaxis, :]
    
    frequency_lineage[x, :] = np.sum(locs, axis = 0)
    
    pfunc = partial(sub_func, x = x, days_prop = days_prop, unique_days_prop = unique_days_prop, lineages_all = lineages_all, unique_lineage = unique_lineage_timeframe)
    try:
        res = list(jb.Parallel(n_jobs = -1, backend = "loky")(jb.delayed(pfunc)(d) for d in range(len(unique_days_prop))))
    except:
        try:
            res = list(jb.Parallel(n_jobs = -1, backend = "multiprocessing")(jb.delayed(pfunc)(d) for d in range(len(unique_days_prop)))) 
        except: 
            res = list(jb.Parallel(n_jobs = -1, prefer = "threads")(jb.delayed(pfunc)(d) for d in range(len(unique_days_prop)))) 

    for s in range(len(unique_days_prop)):
        frequency_lineage[x, s] = res[s]
"""

"""Interpolate in-between dates"""
if seq_thres is not None:
    sdate = datetime.strptime(unique_days_prop[0], "%Y-%m-%d")
    edate = datetime.strptime(unique_days_prop[len(unique_days_prop) - 1], "%Y-%m-%d")
    date_list = pd.date_range(sdate,edate-timedelta(days=1),freq='d').strftime('%Y-%m-%d').tolist()
    if unique_days_prop[len(unique_days_prop) - 1] not in date_list:
        date_list.append(unique_days_prop[len(unique_days_prop) - 1])
    
    x = []
    for i in range(len(unique_days_prop)):
        da = np.array(unique_days_prop[i].split("-")).astype(str)
        da = da[~(da == "")]
        da = da[~(da== " ")]
        da = list(da)
        reformat = False
        if int(da[1]) < 10 and len(da[1]) != 2:
            da[1] = "0%d"%int(da[1])
            reformat = True
        if int(da[2]) < 10 and len(da[2]) != 2:
            da[2] = "0%d"%int(da[2])
            reformat = True
        if reformat:
            dr = "-".join(da)
            if dr in date_list:
                x.append(date_list.index(dr))
        else:
            x.append(date_list.index(unique_days_prop[i]))
    
    if len(x) != len(unique_days_prop):
        sys.exit("Some dates are not properly fomated in covsonar data: Please only use format Year-month-days \n and put 0 before single digit days/months e.g. May 3rd, 2022 = 2022-03-03 \n We recommend to first replace covsonar data with the output of (assuming that country = Brazil): Rscript scripts/check_dataset.R path/to/covsonar_data_Brazil.tsv Brazil path/to_new_filename/covsonar_data_Brazil_checked.tsv")
        
    x = np.array(x)
    y = frequency_lineage
    sub_x = np.arange(0, len(date_list))
    f = interp1d(x, y, axis = -1)
    props_fill = f(sub_x)
    unique_days_prop = np.array(date_list)
    frequency_lineage = props_fill
        
""" Re-Normalize variant proportion data to reflect daily proportion distribution """
NormProp = np.sum(frequency_lineage, axis = 0)
freq_rounded = np.round(frequency_lineage,decimals = 10)
frequency_lineage = np.divide(freq_rounded, NormProp, out = np.zeros(freq_rounded.shape), where = NormProp != 0)*100

### check
""" Save frequency data """
freq_dic = {}
freq_dic["date"] = list(unique_days_prop)

try:
    """ Save frequency calendar week """
    KW_all = list(covsonar_data["week_num"].values.astype(str))
    KW_kept = []
    for i in range(len(freq_dic["date"])):
        if freq_dic["date"][i] in days_prop0:
            KW_kept.append(KW_all[list(days_prop0).index(freq_dic["date"][i])])
        else:
            KW_kept.append(np.nan)
    freq_dic["week_num"] = KW_kept       
except:
    print("Consonar data does not have the column week_num of the calendar week")
    
for x in range(len(unique_lineage_timeframe)):
    freq_dic[unique_lineage_timeframe[x]] = frequency_lineage[x, :]

freq_df = pd.DataFrame(freq_dic, index = np.arange(0, len(freq_dic["date"])))
freq_df.to_csv(sys.argv[3]) ### output file

if len(sys.argv) > 5 :
    variant_mut_data = pd.read_csv(sys.argv[5])
    variant_x_name_orig = np.array(variant_mut_data["lineage"].values).astype(str)
    mut_x_sites_orig = np.array(variant_mut_data["RBD_NTD_mutations"].values).astype(str)
    mut_x_sites_dic = {}
    AA_change_dic = {}
    for i in range(len(variant_x_name_orig)):
        x = variant_x_name_orig[i]
        mut_x = mut_x_sites_orig[i]
        split_mut = mut_x.split("/")
        aa_x = {}
        pos_list = []
        for mut in split_mut:
            pos0 = re.findall(r'\d+', mut)
            if len(pos0) == 1:
                pos = str(pos0[0])
                if pos not in list(aa_x.keys()):
                    aa_x[pos] = [mut]
                    pos_list.append(pos)
                else:
                    aa_x[pos].append(mut)
                
        mut_x_sites_dic[x] = pos_list
        AA_change_dic[x] = aa_x
    
    mut_file = open(sys.argv[6], "wb")
    pickle.dump({"positions": mut_x_sites_dic, "AA_changes":AA_change_dic, "Group":variant_x_name_orig}, mut_file)
    mut_file.close()
    
    
