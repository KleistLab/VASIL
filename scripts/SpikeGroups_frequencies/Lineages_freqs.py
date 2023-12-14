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
"""Simulation timeframe"""
date_start = sys.argv[2]

inds_keep = np.array(days_prop0)!="nan"
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
        d = np.array(unique_days_prop_all[i].split("-"))
        d = d[~(d == "")]
        d = d[~(d == " ")]
        if int(d[1]) < 10:
            d[1] = "0%d"%int(d[1])
        if int(d[2]) < 10:
            d[2] = "0%d"%int(d[2])
        
        dr  = "-".join(d)
        unique_days_prop_sub.append(dr)
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
    

def sub_func(s, x, days_prop, unique_days_prop, lineages_all, unique_lineage):
    res = np.sum((days_prop == unique_days_prop[s]) & (lineages_all == unique_lineage[x]))
    return res

print("Timeline of lineage proportions: %s -- %s"%(unique_days_prop[0], unique_days_prop[-1]))
unique_lineage_timeframe = np.unique(lineages_all)

try:
    seq_thres = int(sys.argv[4])
except:
    seq_thres = None
    
if seq_thres is not None:

    num_seq = []
    for u in range(len(unique_days_prop)):
        count_u = np.sum(days_prop == unique_days_prop[u])
        num_seq.append(count_u)
    
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
frequency_lineage = np.zeros((len(unique_lineage_timeframe), len(unique_days_prop))) # indexing of t correspondis to timeline of infection days_incidence
for x in range(len(unique_lineage_timeframe)):
    
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

"""Interpolate in-between dates"""
if seq_thres is not None:
    sdate = datetime.strptime(unique_days_prop[0], "%Y-%m-%d")
    edate = datetime.strptime(unique_days_prop[len(unique_days_prop) - 1], "%Y-%m-%d")
    date_list = pd.date_range(sdate,edate-timedelta(days=1),freq='d').strftime('%Y-%m-%d').tolist()
    if unique_days_prop[len(unique_days_prop) - 1] not in date_list:
        date_list.append(unique_days_prop[len(unique_days_prop) - 1])

    x = np.array([date_list.index(d) for d in unique_days_prop])
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
    for i in range(len(variant_x_name_orig)):
        mut_i = re.findall(r'\d+', mut_x_sites_orig[i])
        mut_x_sites_dic[variant_x_name_orig[i]] = mut_i
    
    mut_file = open(sys.argv[6], "wb")
    pickle.dump({"positions": mut_x_sites_dic, "Group":variant_x_name_orig}, mut_file)
    mut_file.close()
    
        
