#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 30 12:02:05 2023

@author: raharinirina
"""

import sys
import pandas as pd
import pdb
import numpy as np
import pickle

try:
    df_Vaccines = pd.read_csv(sys.argv[1])
    print(df_Vaccines["date"])
except:
    df_Vaccines = pd.read_csv(sys.argv[1], sep="\t")

date_start = str(sys.argv[2])
date_end = str(sys.argv[3])
switch_vacc = str(sys.argv[4])
vacc_considered = str(sys.argv[5])

save_mut_to = str(sys.argv[-2])
save_to = str(sys.argv[-1])

### Drop dates between 
if date_start in df_Vaccines["date"].tolist():
    df_Vaccines.drop(index = df_Vaccines.index[:list(df_Vaccines["date"]).index(date_start)], inplace = True)

if date_end in df_Vaccines["date"].to_list():
    df_Vaccines.drop(index = df_Vaccines.index[list(df_Vaccines["date"]).index(date_end)+1:], inplace = True)


df_Vaccines.drop(columns = "impfungen_boost1", inplace = True)
df_Vaccines.drop(columns = "impfungen_boost2", inplace = True)
df_Vaccines.drop(columns = "impfungen_boost3", inplace = True)
df_Vaccines.drop(columns = "impfungen_boost4", inplace = True)


dates_vacc = df_Vaccines["date"].to_list()
### Extract relevant vaccine columns
Columns = df_Vaccines.columns

Timeline = {}
Timeline["date"] = dates_vacc
col_kept = []
for col in Columns:
    if ("boost" in col):
    #if ("gi" in col) or ("boost" in col):
        if ("kumulativ" not in col) & ("impfungen" in col) :
            for v_n in ("biontech", "moderna"): #"astra", "novavax", "valneva", "johnson"
                if v_n in col:
                    if switch_vacc not in ("none", "None"): ### Not used (might be Deprecated), Still need to be updates later if is needed
                        print("kept", col)
                        col_kept.append(col)
                        col_new = col
                        vacc_wt = np.zeros(len(dates_vacc))
                        vacc_wt[:dates_vacc.index(switch_vacc)] = np.array(df_Vaccines[col])[:dates_vacc.index(switch_vacc)]
                        Timeline["Wuhan-Hu-1*_as_"+col_new] = vacc_wt
                        
                        vacc_var1 = np.zeros(len(dates_vacc))
                        vacc_var1[dates_vacc.index(switch_vacc):] = 0.5*np.array(df_Vaccines[col])[dates_vacc.index(switch_vacc):]
                        Timeline["BA.4*_as_"+col_new] = vacc_var1
                        
                        vacc_var2 = np.zeros(len(dates_vacc))
                        vacc_var2[dates_vacc.index(switch_vacc):] = 0.5*np.array(df_Vaccines[col])[dates_vacc.index(switch_vacc):]
                        Timeline["BA.5*_as_"+col_new] = vacc_var2
                    else:
                        if vacc_considered == "all_boosters":
                            print("kept", col)
                            col_kept.append(col)
                            col_new = col
                            if "bivalent" in col:
                                Timeline["BA.5*_as_"+col_new] = np.array(df_Vaccines[col])
                            else:
                                Timeline["Wuhan-Hu-1*_as_"+col_new] = np.array(df_Vaccines[col])
                        elif vacc_considered == "bivalent_boosters":
                            if "bivalent" in col:
                                print("kept", col)
                                col_kept.append(col)
                                col_new = col
                                Timeline["BA.5*_as_"+col_new] = np.array(df_Vaccines[col])
                                

print("Vaccines condidered:", col_kept)
df_Timeline = pd.DataFrame(Timeline)
df_Timeline.to_csv(save_to+"/Vaccination_Timeline.csv")

Total_vacc = {}
Total_vacc["date"] = dates_vacc
Total_vacc["vaccinated"] = np.sum(df_Timeline.to_numpy()[:, np.array(df_Timeline.columns) != "date"], axis = 1)

df_total = pd.DataFrame(Total_vacc)    
df_total.to_csv(save_to+"/Vaccination_Total.csv")