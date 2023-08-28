#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
from datetime import datetime
import joblib as jb
from functools import partial
import sys

"""Load Infection Data: used to make sure that infection timeline matches with lineage proportion timeline"""

Population_Data = pd.read_csv(sys.argv[1])

### timeline of infection dates
days_incidence = list(Population_Data['date']) 

#### covsonar data 
try:
	covsonar_data = pd.read_csv(sys.argv[2], sep = "\t")
except:
	covsonar_data = pd.read_csv(sys.argv[2])
	
days_prop = covsonar_data["date"].values.astype(str)
lineages_all = covsonar_data["lineage"].values.astype(str)
unique_lineage = np.unique(lineages_all)

"""Start computing Variant-propotions from the first day July 2021"""
date_start = sys.argv[3]
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
if extra != 0:
    for k in range(extra):
        extra_days.append(unique_days_prop[len(days_incidence[where_first_day:]) + k])
        
total_days = len(days_incidence[where_first_day:])+extra
frequency_lineage = np.zeros((len(unique_lineage), total_days)) # indexing of t correspondis to timeline of infection days_incidence
        
def sub_func(s, x, days_prop, days_incidence, lineages_all, unique_lineage):
    res = np.sum((days_prop == days_incidence[where_first_day + s]) & (lineages_all == unique_lineage[x]))
    return res

def sub_func2(k, x, days_prop, days_incidence, lineages_all, unique_lineage):
    res = np.sum((days_prop == unique_days_prop[len(days_incidence[where_first_day:]) + k]) & (lineages_all == unique_lineage[x]))
    return res
    
for x in range(len(unique_lineage)):
    pfunc = partial(sub_func, x = x, days_prop = days_prop, days_incidence = days_incidence, lineages_all = lineages_all, unique_lineage = unique_lineage)
    res = list(jb.Parallel(n_jobs = -1, prefer = "threads")(jb.delayed(pfunc)(d) for d in range(len(days_incidence[where_first_day:])))) 
    
    for s in range(len(days_incidence[where_first_day:])):
        frequency_lineage[x, s] = res[s]
        
    if len(unique_days_prop)>len(days_incidence[where_first_day:]):
        pfunc2 = partial(sub_func2, x = x, days_prop = days_prop, days_incidence = days_incidence, lineages_all = lineages_all, unique_lineage = unique_lineage)
        res2 = list(jb.Parallel(n_jobs = -1, prefer = "threads")(jb.delayed(pfunc2)(d) for d in range(extra))) 
        for k in range(extra):
            frequency_lineage[x, len(days_incidence[where_first_day:]) + k] = res2[k]
 
""" Re-Normalize variant proportion data to reflect daily proportion distribution """
NormProp = np.sum(frequency_lineage, axis = 0)
freq_rounded = np.round(frequency_lineage,decimals = 10)
frequency_lineage = np.divide(freq_rounded, NormProp, out = np.zeros(frequency_lineage.shape), where = NormProp != 0)*100

### check
""" Save frequency data """
freq_dic = {}
freq_dic["date"] = list(days_incidence[where_first_day:]) + list(extra_days)                       
for x in range(len(unique_lineage)):
    freq_dic[unique_lineage[x]] = frequency_lineage[x, :]

freq_df = pd.DataFrame(freq_dic, index = np.arange(0, frequency_lineage.shape[1]))
freq_df.to_csv(sys.argv[4]) ### output file
