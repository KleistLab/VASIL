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
from datetime import datetime
import joblib as jb
from functools import partial
import sys
import pdb


"""Load covsonar_data"""
try:
	covsonar_data = pd.read_csv(sys.argv[1], sep = "\t")
except:
	covsonar_data = pd.read_csv(sys.argv[1])

try:
    days_prop = covsonar_data["Date"].values.astype(str)
except:
    try:
        days_prop = covsonar_data["date"].values.astype(str)
    except:
        sys.exit("Date or date column not found in covsonar data")
    
lineages_all = covsonar_data["lineage"].values.astype(str)
"""Simulation timeframe"""
date_start = sys.argv[2]


"""Start computing Variant-proportions"""
# iniializing variant proportion for all lineages
unique_days_prop_all = list(np.unique(days_prop))
### make sure that dates are sorted and not nan
unique_days_prop_all = [x for x in unique_days_prop_all if x!= "nan"]
unique_days_prop_all.sort(key = lambda date: datetime.strptime(date, "%Y-%m-%d")) 

unique_days_prop = np.array(unique_days_prop_all[unique_days_prop_all.index(date_start):])
lineages_all = np.array(lineages_all)[unique_days_prop_all.index(date_start):]
        
def sub_func(s, x, days_prop, unique_days_prop, lineages_all, unique_lineage):
    res = np.sum((days_prop == unique_days_prop[s]) & (lineages_all == unique_lineage[x]))
    return res


unique_lineage_timeframe = np.unique(lineages_all)
""" Remove NONE and UNASSIGNED """
if "NONE" in unique_lineage_timeframe:
    unique_lineage_timeframe = np.array(unique_lineage_timeframe)[np.array(unique_lineage_timeframe) != "NONE"]
if "UNASSIGNED" in unique_lineage_timeframe:
    unique_lineage_timeframe = np.array(unique_lineage_timeframe)[np.array(unique_lineage_timeframe) != "UNASSIGNED"]


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
        KW_kept.append(KW_all[list(days_prop).index(freq_dic["date"][i])])
        freq_dic["week_num"] = KW_kept       
except:
    print("Consonar data does not have the column week_num of the calendar week")
    
for x in range(len(unique_lineage_timeframe)):
    freq_dic[unique_lineage_timeframe[x]] = frequency_lineage[x, :]


freq_df = pd.DataFrame(freq_dic, index = np.arange(0, len(freq_dic["date"])))
freq_df.to_csv(sys.argv[3]) ### output file
