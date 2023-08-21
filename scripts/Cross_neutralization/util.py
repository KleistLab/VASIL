#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import numpy.ma as ma
import joblib as jb
from functools import partial

def sub_Bind(d, tiled_esc, Where_Mut, Where_Cond):
    Inter_Cond_Mut = Where_Mut & Where_Cond[np.newaxis, d, :]
    Bind_list_d = np.prod(1 - tiled_esc[:, np.newaxis, :]*Inter_Cond_Mut, axis = (1,2))
    Missing_cond_data_d = ~np.any(np.any(Where_Mut & Where_Cond[np.newaxis, d, :], axis = 2), axis = 1)
    return [Bind_list_d, Missing_cond_data_d]

def FR_xy(i, mut_sites, mut_bool, escape_ab_dic, ab, variant_name, mut_sites_per_variant, EF_func = "MEAN", GM = False, quiet = True):
    vars_num = mut_bool[i+1:, :].shape[0]
    
    test_i = np.tile(mut_bool[i, :], (vars_num, 1))
    diff_sites = (test_i ^ mut_bool[i+1:, :])

    escape_data_ab = escape_ab_dic["escape_data_ab"]
    conditions = escape_ab_dic["conditions"]
    ab_sub_list = escape_ab_dic["ab_sub_list"]   
    escape_sites = escape_ab_dic["escape_sites"]
    IC50_list = escape_ab_dic["IC50_list"]
    
    tiled_esc = np.tile(escape_data_ab, (vars_num, 1))
    Bind_list  = np.ones((vars_num, len(conditions)))
    Missing_cond_data = np.zeros((vars_num, len(conditions)), dtype = bool)
    Where_Cond = conditions[:, np.newaxis] == ab_sub_list[np.newaxis, :]
    tiled_mut = ma.array(np.tile(mut_sites, (vars_num, 1)), mask = ~diff_sites)
    Where_Mut = tiled_mut[:, :, np.newaxis] == escape_sites[np.newaxis, np.newaxis,  :]
    
    """Parallel codes"""
    pfunc = partial(sub_Bind, tiled_esc = tiled_esc, Where_Mut = Where_Mut, Where_Cond = Where_Cond)
    try:
        jb_res = list(jb.Parallel(n_jobs = 10)(jb.delayed(pfunc)(d) for d in range(len(conditions))))
    except:
        jb_res = list(jb.Parallel(n_jobs = 10, prefer = "threads")(jb.delayed(pfunc)(d) for d in range(len(conditions))))
    
    for d in range(len(conditions)):
        Bind_list[:, d]   = np.array(jb_res[d][0])
        Missing_cond_data[:, d] = np.array(jb_res[d][1])
    
    """
    for d in range(len(conditions)):
        print(d, len(conditions))
        Inter_Cond_Mut = Where_Mut & Where_Cond[np.newaxis, d, :]
        Bind_list[:, d] = np.prod(1 - tiled_esc[:, np.newaxis, :]*Inter_Cond_Mut, axis = (1,2))  
        Missing_cond_data[:, d] = ~np.any(np.any(Where_Mut & Where_Cond[np.newaxis, d, :], axis = 2), axis = 1)
    """
    
    """
    ### Crashes: Unable to allocate 38.1 GB for an array with shape (400, 418, 109, 2247) and data type bool
    Inter_Cond_Mut = Where_Mut[:, :, np.newaxis] & Where_Cond[np.newaxis, :, :]
    Bind_list = np.prod(1 - tiled_esc[:, np.newaxis, np.newaxis, :]*Inter_Cond_Mut, axis = (1,3))  
    Missing_cond_data = ~np.any(np.any(Inter_Cond_Mut, axis = 3), axis = 1)
    """
    
    retained_binding = Bind_list
    FR_list = np.maximum((0.4/IC50_list[np.newaxis, :])*((1/retained_binding) - 1), 1)# DMS data was performed at fixed antibody concentration of 0.4 micro-gram/L (IC50 values are in micro-gram/L)
    #FR could be less than 1 happen if the IC50 is close to 0.4 and thus the binding retained is not reliable, we set those to 1
    FR_list = ma.array(FR_list, mask = Missing_cond_data) # do not take into account the missing data
    FR_ab = FR_list
    
    Missed = Missing_cond_data
    Greater_one = []
    
    return FR_ab, Missed, Greater_one 

def cross_reactivity(variant_name, escape_per_sites, Ab_classes, mut_sites_per_variant, EF_func = "MEAN", GM = False, quiet = True):
    FRxy = {}
    Missed = []
    Greater_one = []
    
    ### extract all sites
    mut_sites = []
    for var in list(mut_sites_per_variant.keys()):
        if var in variant_name:
            mut_sites += list(np.array(mut_sites_per_variant[var]))
    
    mut_sites = np.unique(mut_sites).astype(str)
    ### Constriuct a boolean array for location of mutation for all variants    
    mut_bool = np.zeros((len(variant_name), len(mut_sites)), dtype = bool)
    
    #print(variant_name)
    for i in range(len(variant_name)):
        for j in range(len(mut_sites)): 
            mut_bool[i, j] = mut_sites[j] in list(np.array(mut_sites_per_variant[variant_name[i]]).astype(str))
        
    for ab in Ab_classes:
        if not quiet:
            print("AB", ab)
            
        FRxy_ab = np.ones((len(variant_name), len(variant_name)))
        Missing_Cond = np.zeros((len(variant_name), len(variant_name)))
        
        escape_ab_dic = {}
        escape_data = escape_per_sites["mut_escape"].values
        where_ab_group = (escape_per_sites["group"].values).astype(str) == str(ab)
        escape_ab_dic["ab_sub_list"]   = (escape_per_sites["condition"].values).astype(str)[where_ab_group]

        sub_table = escape_per_sites.loc[(escape_per_sites['condition'][where_ab_group].drop_duplicates().index), ['condition', 'IC50']]
        escape_ab_dic["conditions"] = sub_table["condition"].values
        escape_ab_dic["IC50_list"] = sub_table["IC50"].values
        
        escape_ab_dic["escape_sites"] = ((escape_per_sites["site"].values).astype(str))[where_ab_group]
        
        escape_ab_dic["escape_data_ab"] = escape_data[where_ab_group]
        
        for i in range(len(variant_name)):
            FR, missed, gOne = FR_xy(i, mut_sites, mut_bool, escape_ab_dic, ab, variant_name, mut_sites_per_variant, GM = GM, quiet = quiet)
            if EF_func == "MEAN":
                FRxy_ab[i, i+1:] = np.mean(FR, axis = 1)
                FRxy_ab[i+1:, i] = FRxy_ab[i, i+1:]
                
            elif EF_func == "PROD":
                FRxy_ab[i, i+1:] = np.prod(FR, axis = 1)
                FRxy_ab[i+1:, i] = FRxy_ab[i, i+1:]
            
            Missing_Cond[i, i+1:] = np.all(missed, axis = 1)
            Missing_Cond[i+1:, i] = Missing_Cond[i, i+1:]
            
            Missed = []
            Greater_one += gOne
        
        #FRxy[ab] = ma.array(FRxy_ab, mask = Missing_Cond.astype(bool), fill_value = np.nan)
        FRxy_ab[Missing_Cond.astype(bool)] = 1
        FRxy[ab] = FRxy_ab
        
   
    Missed = np.unique(np.array(Missed))
    Greater_one = np.unique(np.array(Greater_one))           
    return FRxy, Missed, Greater_one
