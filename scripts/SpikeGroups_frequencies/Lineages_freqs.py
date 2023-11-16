#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
from datetime import datetime
import joblib as jb
from functools import partial
import sys
import pdb

"""Load Infection Data: used to make sure that infection timeline matches with lineage proportion timeline"""

Population_Data = pd.read_csv(sys.argv[1])

### timeline of infection dates
days_incidence = list(Population_Data['date']) 

"""Simulation timeframe"""
date_start = sys.argv[3]
#date_end = sys.argv[4]

#### covsonar data 
try:
	covsonar_data = pd.read_csv(sys.argv[2], sep = "\t")
except:
	covsonar_data = pd.read_csv(sys.argv[2])

try:
    days_prop = covsonar_data["Date"].values.astype(str)
except:
    try:
        days_prop = covsonar_data["date"].values.astype(str)
    except:
        sys.exit("Date or date column not found in covsonar data")
    
lineages_all = covsonar_data["lineage"].values.astype(str)

"""Start computing Variant-proportions"""
where_first_day = list(days_incidence).index(date_start) 

# initializing variant proportion for all lineages
unique_days_prop_all = list(np.unique(days_prop))
### make sure that dates are sorted
unique_days_prop_all.sort(key = lambda date: datetime.strptime(date, "%Y-%m-%d")) 

unique_days_prop = np.array(unique_days_prop_all[unique_days_prop_all.index(date_start):])
if len(unique_days_prop)>len(days_incidence[where_first_day:]):
    extra = len(unique_days_prop) - len(days_incidence[where_first_day:])
else:
    extra = 0

extra_days = []
"""
### Restricted to frequency timeline to date_end was necessary in older versions (this part of code will be deleted later)
if extra != 0:
    for k in range(extra):
        if len(days_incidence[where_first_day:]) + k <= list(unique_days_prop).index(date_end):
            extra_days.append(unique_days_prop[len(days_incidence[where_first_day:]) + k])
"""    

if extra != 0:
    for k in range(extra):
        if len(days_incidence[where_first_day:]) + k <= len(unique_days_prop) - 1:
            extra_days.append(unique_days_prop[len(days_incidence[where_first_day:]) + k])
            
total_days = len(days_incidence[where_first_day:]) + len(extra_days)
        
def sub_func(s, x, days_prop, days_incidence, lineages_all, unique_lineage):
    res = np.sum((days_prop == days_incidence[where_first_day + s]) & (lineages_all == unique_lineage[x]))
    return res

def sub_func2(k, x, days_prop, days_incidence, lineages_all, unique_lineage):
    res = np.sum((days_prop == unique_days_prop[len(days_incidence[where_first_day:]) + k]) & (lineages_all == unique_lineage[x]))
    return res

""" 
### Restricting to lineages until date_end was necessary in older versions but now the feature is directly implemented in scripts/mutationprofile/mutation_profile.R (this part of code will be deleted later)
between_dates = unique_days_prop_all[unique_days_prop_all.index(date_start):unique_days_prop_all.index(date_end)+1]
to_keep = []
for i in range(len(days_prop)):
    if days_prop[i] in between_dates:
        to_keep.append(i)
unique_lineage_timeframe = np.unique(lineages_all[np.array(to_keep)])
"""
unique_lineage_timeframe = np.unique(lineages_all)

"""Compute lineage frequencies"""
frequency_lineage = np.zeros((len(unique_lineage_timeframe), total_days)) # indexing of t correspondis to timeline of infection days_incidence
for x in range(len(unique_lineage_timeframe)):
    pfunc = partial(sub_func, x = x, days_prop = days_prop, days_incidence = days_incidence, lineages_all = lineages_all, unique_lineage = unique_lineage_timeframe)
    try:
        res = list(jb.Parallel(n_jobs = -1, backend = "loky")(jb.delayed(pfunc)(d) for d in range(len(days_incidence[where_first_day:]))))
    except:
        try:
            res = list(jb.Parallel(n_jobs = -1, backend = "multiprocessing")(jb.delayed(pfunc)(d) for d in range(len(days_incidence[where_first_day:])))) 
        except: 
            res = list(jb.Parallel(n_jobs = -1, prefer = "threads")(jb.delayed(pfunc)(d) for d in range(len(days_incidence[where_first_day:])))) 

    
    for s in range(len(days_incidence[where_first_day:])):
        frequency_lineage[x, s] = res[s]
        
    if len(unique_days_prop)>len(days_incidence[where_first_day:]):
        pfunc2 = partial(sub_func2, x = x, days_prop = days_prop, days_incidence = days_incidence, lineages_all = lineages_all, unique_lineage = unique_lineage_timeframe)
        try:
            res2 = list(jb.Parallel(n_jobs = -1, backend = "loky")(jb.delayed(pfunc2)(d) for d in range(extra))) 
        except:
            try:
                res2 = list(jb.Parallel(n_jobs = -1, backend = "multiprocessing")(jb.delayed(pfunc2)(d) for d in range(extra))) 
            except:
                res2 = list(jb.Parallel(n_jobs = -1, prefer = "threads")(jb.delayed(pfunc2)(d) for d in range(extra))) 
                
        for k in range(len(extra_days)):
            frequency_lineage[x, len(days_incidence[where_first_day:]) + k] = res2[k]
 
""" Re-Normalize variant proportion data to reflect daily proportion distribution """
NormProp = np.sum(frequency_lineage, axis = 0)
freq_rounded = np.round(frequency_lineage,decimals = 10)
frequency_lineage = np.divide(freq_rounded, NormProp, out = np.zeros(freq_rounded.shape), where = NormProp != 0)*100

### check
""" Save frequency data """
freq_dic = {}
freq_dic["date"] = list(days_incidence[where_first_day:]) + list(extra_days)

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
freq_df.to_csv(sys.argv[5]) ### output file
