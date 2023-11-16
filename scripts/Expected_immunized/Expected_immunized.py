#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import pickle
import re
import joblib as jb
from functools import partial
import sys
import warnings
import pdb
"""Load Infection Data"""
Population_Data = pd.read_csv(sys.argv[1])

"""Load delta relevant cross_neutralization files"""
file1 = open(sys.argv[2], "rb") # premade simulations
Cross_with_delta_validation = pickle.load(file1)
variants_x_names_show = Cross_with_delta_validation["variant_list"]
Cross_with_delta_validation.pop("variant_list")
file1.close()
Ab_classes = list(Cross_with_delta_validation.keys()) ### the same for all cross reactivity files computed in the pipeline
try:
    file1 = open(sys.argv[3], "rb") # load cross reactivity if parameter is not a directory
    Cross_react_dic = pickle.load(file1)
    variants_in_cross = Cross_react_dic["variant_list"]
    Cross_react_dic.pop("variant_list")
    file1.close()
    run_group=False
except:
    # result directory is instead provided
    run_group = True

"""Spike groups and frequencies"""
file1 = open(sys.argv[4], "rb") 
SpikeGroups_list = pickle.load(file1)["names"]
file1.close()
frequency_spk_df = pd.read_csv(sys.argv[5])

"""Escape fraction data """
Escape_Fraction = pd.read_csv(sys.argv[6])

"""Load vaccing efficacy data for fitting """
Load_Delta = pd.read_excel(sys.argv[7], engine='openpyxl')

"""Load total population and re-adjust infection timeline"""
total_population = float(sys.argv[8])

date_start = str(sys.argv[9])
Population_Data = Population_Data.drop(index = Population_Data.index[:list(Population_Data['date']).index(date_start)])

date_end = sys.argv[10]
if date_end in list(Population_Data['date']):
    where_last_day = list(Population_Data['date']).index(date_end) + 1
else:
    where_last_day = len(list(Population_Data['date']))

t = np.arange(1, len(Population_Data['date'][:where_last_day]) + 1, 1)

cases_col = str(sys.argv[len(sys.argv) - 1])
try:
    infection_data_corrected = Population_Data[cases_col].values[:where_last_day]
except:
    sys.exit("The column %s was not found in case ascertainment data"%cases_col)
    
t_dates = Population_Data['date'].values[:where_last_day]
days_incidence = list(Population_Data['date'][:where_last_day]) 
        
def Antibody_ranges(thalf_vec, tmax_vec, t, Ab_classes):
    N = len(Ab_classes) # number of antibody classes
    is_log = False # if True, it returns the log of the antibody concentration
    dataname = "Ab_%d"%N
    solver = "lm" # root solver method for finding absorption rate ka (see scipy.optimize.root)
    c_t_vec = np.zeros((len(thalf_vec), len(tmax_vec), len(Ab_classes), len(t)))
    c_dframe_dic = {}
    for m in range(len(thalf_vec)):
        for n in range(len(tmax_vec)):
            t_half = thalf_vec[m]*np.ones(N) # antibody half-life for all antibody classes, respecively
            t_max = tmax_vec[n]*np.ones(N) # time to the peak of antibody concentration for each antibody class
            params_dic = {"t_max":t_max, "t_half":t_half}
            c_t, c_dframe, ka, ke, c_max_dic = Antibody(t = t, params_dic = params_dic, is_log = is_log, Ab_names = Ab_classes, ka_solver = solver)
            c_t_vec[m, n, :, :] = c_t
            c_dframe_dic["(%d, %d)"%(m, n)] = c_dframe
            
    return c_t_vec, c_dframe_dic, dataname

def Antibody(t, params_dic, is_log = False, Ab_names = None, ka_solver = "lm"):
    """
    @brief: Compute Antibody Concentration as a function of time for N antibody classes
    
    Parameters
    ----------
    t : T time points ndarray (T) 
    params_dic : dictionary of parameters
                params_dic["t_max"] = t_{max}, ndarray (N, )
                params_dic["t_half"] = t_{1/2}, ndarray (N, ) 
    is_log : bool, optional
             True if return the log of concentration. The default is False.
    ka_solver : root solver method for finding absorption rate k_a (see scipy.optimize.root), optional. The default is lm

    Returns
    -------
    Antibody concentration at time t.

    """
    t_max = params_dic["t_max"]
    t_half = params_dic["t_half"]
    
    # Antibody elimination rate
    ke = np.log(2)/t_half
    
    # Antibody absorption rate
    guess = np.ones(len(t_max))
    warnings.filterwarnings("ignore")
    ka = root(ka_solve, guess, args = (ke, t_max), method = ka_solver).x
    
    if not np.all(np.isclose(ka_solve(ka, ke, t_max), np.zeros(len(t_max)))):
        print("\n k_a was found correctly:", np.all(np.isclose(ka_solve(ka, ke, t_max), np.zeros(len(t_max)))), "\n", params_dic)

    # Compute Normalized Concentration
    c_max = (np.exp(- ke*t_max) - np.exp(- ka*t_max))
    c_t = (np.exp(- ke[:, np.newaxis]*t) - np.exp(- ka[:, np.newaxis]*t))/c_max[:, np.newaxis]
    
    if is_log:
        c_t = np.log(c_t)

    # Build pandas dataframe‚
    df = {}
    c_max_dic = {}
    df["Days"] = t
    for i in range(len(t_max)):
        if Ab_names is None:
            df["Ab class %d"%(i+1)] = c_t[i, :]
            c_max_dic["Ab class %d"%(i+1)] = c_max[i]
        else:
            df[Ab_names[i]] = c_t[i, :]
            c_max_dic[Ab_names[i]] = c_max[i]
        
        
    df = pd.DataFrame(df)
             
    return c_t, df, ka, ke, c_max_dic

def ka_solve(ka, ke, t_max):
    if np.all(ka)>0:
        res = np.divide(t_max*(ka - ke) - (np.log(ka) - np.log(ke)), (ka - ke), out = np.ones(len(ke)), where = (ka - ke)!=0)
    else:
        res = 1
    return res
        

"""Fit IC50 parameter"""
from scipy.optimize import root

def vaccine_efficacy(x, ic50):
    return(x/(x + ic50))

def efficacy_n_antibodies(x, ic50):
    res = 1
    for i in range(len(x)):
        ve = vaccine_efficacy(x[i], ic50[i])
        res *= (1 - ve)
    return(1 - res)

def sqrt_diff_FR(ic50, days_list, FR, ve_data, n, c_dframe):
    res = 0
    for d in range(len(ve_data)):
        data = ve_data[d]
        days = days_list[d]
        ve_estimate = np.zeros(len(days))
        for i in range(len(data)):
            antibody_level = c_dframe.loc[int(days[i]) - 1][1:n+1]
            ve_estimate[i] = efficacy_n_antibodies(antibody_level, np.array(FR)*ic50)
        res += np.linalg.norm(data-ve_estimate[0:len(data)])
    return(res)


"""Compute fold change deviation for mean IC50"""
IC50_group = Escape_Fraction.groupby('condition', as_index=False).first()[['condition', 'IC50', 'group']]
mean_IC50_per_group = IC50_group.groupby('group')['IC50'].mean().reset_index()
total_mean = mean_IC50_per_group['IC50'].mean()
mean_IC50_per_group['fold_change'] = mean_IC50_per_group['IC50']/total_mean
ntd_row = {'group': 'NTD','IC50': 1, 'fold_change': 1}
mean_IC50_per_group = pd.concat([mean_IC50_per_group, pd.DataFrame(ntd_row, index=[10])])

"""Load fold change IC50 in data"""
FC_ic50_dic = {Ab_classes[i]:(mean_IC50_per_group["fold_change"].values[mean_IC50_per_group["group"] == Ab_classes[i]])[0] for i in range(len(Ab_classes))}

"""Extract information from Clinical data"""
def transform(x):
    x[x<0] = 0
    return x

def extract_yerr(x, CI):
    lower_diff = np.minimum(transform(x), np.abs(x - CI[:, 0]))
    upper_diff = np.minimum(transform(x), np.abs(CI[:, 1] - x))
    return np.array([lower_diff, upper_diff]), CI[:, 1], CI[:, 0]
    

"""Fit Delta Vaccine Data """
days_fitting = []
ve_fitting = []

Delta_Sources = Load_Delta["Source"].values.astype(str)
Delta_Vaccine = Load_Delta["Vaccine"].values.astype(str)

All_Days_Delta = np.array([])
All_Delta_Data = np.array([])
All_Days_xerr_Delta = np.array([])
All_Delta_yerr = np.array([])

""" Set up Disease status to restrict the studies on any infected people"""
keep_status_d = Load_Delta["Disease Status"].values.astype(str) == "Any infection"
keep_method_d = (Load_Delta["Method"].values.astype(str) != "(1 - Adjusted OR)") & (Load_Delta["Method"].values.astype(str) != "(1 - OR)")

u_Delta_Sources = np.unique(Delta_Sources[keep_status_d & keep_method_d])
u_Delta_Vaccine = np.unique(Delta_Vaccine[keep_status_d & keep_method_d])

u_Delta_Sources = u_Delta_Sources[~(u_Delta_Sources  == "nan")] 
u_Delta_Vaccine = u_Delta_Vaccine[~(u_Delta_Vaccine  == "nan")]


Delta_done = []
for source in u_Delta_Sources:
    for vacc in u_Delta_Vaccine:
        where_source = (Delta_Sources == source)&(Delta_Vaccine == vacc)
        
        if (np.sum(where_source)!=0) and ("%s (%s)"%(vacc, source) not in Delta_done):
            
            #if (not re.search("Feikin", source)):
            days_fitting.append(Load_Delta["Days (Mid)"].values[where_source])
            ve_fitting.append(transform(Load_Delta["VE (value)"].values[where_source]))
            
            All_Days_Delta = np.concatenate((All_Days_Delta, Load_Delta["Days (Mid)"].values[where_source]))
            All_Delta_Data = np.concatenate((All_Delta_Data, Load_Delta["VE (value)"].values[where_source]))
            All_Days_xerr_Delta = np.concatenate((All_Days_xerr_Delta, Load_Delta["Days Err (+/-)"].values[where_source]))
            
            x = Load_Delta["VE (value)"].values[where_source]
            CI = np.array([Load_Delta["VE (Lower CI)"].values[where_source], Load_Delta["VE (Upper CI)"].values[where_source]]).T
            yerr, upper_CI, lower_CI = extract_yerr(x = x, CI = CI)
            
            if len(All_Delta_yerr.flatten()) !=0:
                All_Delta_yerr = np.concatenate((All_Delta_yerr , yerr), axis = 1)
            else:
                All_Delta_yerr = yerr
        
        Delta_done.append("%s (%s)"%(vacc, source))
	
def Fitting_IC50(thalf, tmax, t, Ab_classes, Cross_dic, quiet = False):  
   
    
    N = len(Ab_classes)
    t_max = tmax*np.ones(N) # time to the peak of antibody concentration for each antibody class
    t_half = thalf*np.ones(N) # antibody half-life for all antibody classes, respecively
    params_dic = {"t_max":t_max, "t_half":t_half}
    ### Compute  PK 
    is_log = False # if True, it returns the log of the antibody concentration
    solver = "lm" # root solver method for finding absorption rate ka (see scipy.optimize.root)
    c_t, c_dframe, ka, ke, c_max_dic = Antibody(t = t, params_dic = params_dic, is_log = is_log, Ab_names = Ab_classes, ka_solver = solver)
    
    if not quiet:
        print("t_max = %.3f, t_half = %.3f"%(t_max[0], t_half[0]),"\n k_a:", ka, "\n k_e:", ke, "\n c_max", c_max_dic)
    
    ### select FR data for delta computed in the cross_reac_dic_show
    where_wt = list(variants_x_names_show).index("Wuhan-Hu-1")
    where_delta = list(variants_x_names_show).index("Delta: B.1.617.2")
    FR_delta = [Cross_dic[Ab_classes[i]][where_wt, where_delta] for i in range(len(Ab_classes))]
    ### Estimate one IC50 for all ABs
    guess = 0.3
    FC_ic50_list = [FC_ic50_dic[Ab_classes[i]] for i in range(len(Ab_classes))]
    IC50_data = root(sqrt_diff_FR, guess, args = (days_fitting, np.array(FC_ic50_list)*np.array(FR_delta), ve_fitting, len(Ab_classes), c_dframe), method = "lm").x
    #print("fitted", IC50_data)
    
    """Extract the IC50xx -- conserving the fold changes deviation from the mean of each epitope classes"""
    IC50xx = {Ab_classes[i]:IC50_data[0]*FC_ic50_dic[Ab_classes[i]] for i in range(len(Ab_classes))} 
    #print("IC50xx", IC50xx)
    
    return IC50xx, IC50_data[0], FC_ic50_list, FR_delta, c_dframe

def Find_IC50_ranges(thalf_vec, tmax_vec, t, Ab_classes, Cross_dic):
    IC50xx_dic = {}
    for m in range(len(thalf_vec)):
        for n in range(len(tmax_vec)):
            thalf = thalf_vec[m]
            tmax = tmax_vec[n]
            IC50xx, IC50_data, FC_ic50_list, FR_delta, c_dframe = Fitting_IC50(thalf, tmax, t, Ab_classes, Cross_dic, quiet = True)           
            IC50xx_dic["(%d, %d)"%(m, n)] = IC50_data

    IC50xx = np.mean(list(IC50xx_dic.values()))
    mean_IC50xx_dic = {Ab_classes[i]:IC50xx*FC_ic50_dic[Ab_classes[i]] for i in range(len(Ab_classes))}
    return IC50xx_dic, mean_IC50xx_dic


"""Expected Immunity Efficacy as a function of COVI19 variant proportions -- vectorization 2"""
def P_Neut(t, present_variant_index, PK_dframe, tested_variant_list, variant_name, Ab_classes, IC50xx, Cross_react_dic, escape_per_sites = None, mut_sites_per_variant = None):
    x0 = present_variant_index 
    y = list(variant_name).index(tested_variant_list[0])# knowing that for now there is always one variant to be tested at the time
    
    P_neut_ab = [PK_dframe[ab].values[:, np.newaxis]/(PK_dframe[ab].values[:, np.newaxis] + Cross_react_dic[ab][x0, y][np.newaxis, :]*IC50xx[ab]) for ab in Ab_classes]        
    
    P_neut = 1 - np.prod(1 - np.array(P_neut_ab), axis = 0)
        
    return P_neut.T

"""Immunity dynamics using Fast Frourier Transform: scipy.signal.fftconvolve"""
from scipy import signal
def Immunity_dynamics_fftconvolve(t, PK_dframe, infection_data, present_variant_index, tested_variant_list, variant_name, variant_proportion, Ab_classes, 
                                  IC50xx, Cross_react_dic, escape_per_sites = None, mut_sites_per_variant = None, parallel = False, mode_func = None):
    
    stop = min(len(infection_data), variant_proportion.shape[0])
    Infected_l_vect = infection_data[np.newaxis, :stop]*variant_proportion[:, :stop]    
    
    Prob_Neut = P_Neut(t, present_variant_index, PK_dframe, tested_variant_list, variant_name, Ab_classes, IC50xx, Cross_react_dic)
    
    Conv_Mat = signal.fftconvolve(Infected_l_vect, Prob_Neut, axes = 1)[:, :len(t)]
    
    # No normalization
    Expected_Immuned = np.sum(Conv_Mat, axis = 0)
    """
    tested that this gives the as Immunity_dynamics and is 200x faster
    """
    return Expected_Immuned

### Load spikegroups membership file
file = open("Spikegroups_membership.pck", "rb")
Pseudogroup_dic = pickle.load(file)
file.close()
def PNeut_Envelope(s, t, variants, variant_x_names, Cross_react_dic, c_dframe_dic, IC50xx_dic, antigen_list = ["Wuhan-Hu-1"],mean_IC50xx = True):
    
    if mean_IC50xx:
        IC50xx = np.mean(list(IC50xx_dic.values()))
        to_print = "Computing P_Neut, used mean fitted IC50 %.5f"%IC50xx
    else:
        to_print = "Compute P_Neut, used fitted IC5 for each PK paramteters"
    
    to_print = to_print + " for %s vs. %s antigen"%("/".join(variants), antigen_list[s])
    
    num = "%d/%d"%(s+1, len(antigen_list))
    to_print = to_print + " (%s)"%num
    print(to_print)
    
    splited_var = np.array(antigen_list[s].split("/"))
    splited_var = splited_var[~(splited_var == "")]
    splited_var = splited_var[~(splited_var == " ")]
    
    res = np.zeros((len(variants), len(list(c_dframe_dic.keys())), len(t)))
    
    for j in range(len(splited_var)):
        spl_sub = np.array(splited_var[j].split("="))
        spl_sub = spl_sub[~(spl_sub == "")]
        spl_sub = spl_sub[~(spl_sub == " ")]

        lin = spl_sub[0]
        try:
            where_x = list(variant_x_names).index(lin)
        except:
            try:
                where_x = list(variant_x_names).index(Pseudogroup_dic[lin])
            except:
                print("Error in antigen parameter: %s is not present in covsonar data"%lin)

        if len(spl_sub)==1:
            prop_lin = 1/len(splited_var)
        else:
            prop_0 = re.findall(r"[-+]?(?:\d*\.*\d+)", spl_sub[1])[0] ### anything else is error
            prop_lin = float(prop_0)
        
        for i in range(len(variants)):
            where_y = list(variant_x_names).index(variants[i])
            
            for j in range(len(list(c_dframe_dic.keys()))):
                if not mean_IC50xx:
                    IC50xx = IC50xx_dic[list(c_dframe_dic.keys())[j]]
                
                IC50xy = [Cross_react_dic[Ab_classes[i]][where_x, where_y]*IC50xx*FC_ic50_dic[Ab_classes[i]] for i in range(len(Ab_classes))]
                    
                c_dframe = c_dframe_dic[list(c_dframe_dic.keys())[j]]
                for l in range(len(t)): 
                    antibody_level = c_dframe.loc[l][1:]
                    res[i, j, l] += prop_lin*efficacy_n_antibodies(antibody_level, IC50xy) 
    
    return np.min(res, axis = (0, 1)), np.max(res, axis = (0, 1))

"""Compute Antibody concentration over time for a range of t_half and t_max"""
thalf_vec = np.linspace(25, 69, 15) 
tmax_vec = np.linspace(14, 28, 5)
t_conc = np.arange(1, 656, 1) ### used in older codes
c_t_vec, c_dframe_dic, dataname = Antibody_ranges(thalf_vec, tmax_vec, t_conc, Ab_classes)
IC50xx_dic, mean_IC50xx_dic = Find_IC50_ranges(thalf_vec, tmax_vec, t_conc, Ab_classes,  Cross_with_delta_validation)

"""Save Dynamics: The Antibody PK is assumed to be the same for all Epitope Classes"""
PK = {}
PK["Day since activation"] = t_conc
for key in c_dframe_dic.keys():
    key_num = np.array(re.findall(r"\d+", key)).astype(int)
    PK["t_half = %.3f \nt_max = %.3f"%(thalf_vec[key_num[0]], tmax_vec[key_num[1]])] = c_dframe_dic[key]["A"]


PK_df = pd.DataFrame(PK)
PK_df.to_csv("results/PK_for_all_Epitopes.csv")

### spikegroup frequency
print("Fitted mean IC50=%.3f \n IC50 per epitope class is "%mean_IC50xx_dic["NTD"], mean_IC50xx_dic)
try:
	frequency_spk_df.drop(columns = "Unnamed: 0", inplace = True)
except:
	pass

### Make sure proportions axis 0 is aligned with Spikegroup_list
prop_start = list(frequency_spk_df["date"]).index(date_start)
spikegroups_freq = np.zeros((len(SpikeGroups_list), len(t)))
for i in range(len(SpikeGroups_list)):
    if SpikeGroups_list[i]!="Wuhan-Hu-1":
        spikegroups_freq[i, :] = frequency_spk_df["Spike. "+SpikeGroups_list[i]][prop_start:prop_start+len(t)]

NormProp = np.sum(spikegroups_freq, axis = 0)
prop_rounded = np.round(spikegroups_freq,decimals = 10)
spikegroups_proportion = np.divide(prop_rounded, NormProp, out = np.zeros(prop_rounded.shape), where = NormProp != 0)


### end of simulation
def ei_util(Lin_name, variants_in_cross, antigen_list, 
            Cross_react_dic = None, save_pneut=None, w_save=len(sys.argv)-2, 
            var_list_index = None, spikegroups_proportion_adjust=None):
    
    variant_to_sim = [Lin_name]
    EI = {}
    EI["Days"] = days_incidence
    
    Susc = {}
    Susc["Days"] = days_incidence
    
    if antigen_list == ["ALL"]:
        antigen_list = np.array(SpikeGroups_list[var_list_index])
    
    if len(antigen_list) != 0:
        if save_pneut in ("TRUE", "True"):
            VE = {}
            VE["Day since infection"] = t_conc
            pfunc = partial(PNeut_Envelope, t=t_conc, 
                                variants=[Lin_name], 
                                variant_x_names = variants_in_cross, 
                                Cross_react_dic = Cross_react_dic,
                                c_dframe_dic = c_dframe_dic, 
                                IC50xx_dic = IC50xx_dic, 
                                antigen_list = antigen_list, 
                                mean_IC50xx = True, 
                                ) 
            status = False
            try:
                jb_res = list(jb.Parallel(n_jobs = -1, backend = "loky")(jb.delayed(pfunc)(d) for d in range(len(antigen_list))))
                status = True
                print("run joblib.Parallel")
            except:
                try:
                    jb_res = list(jb.Parallel(n_jobs = -1, backend = "multiprocessing")(jb.delayed(pfunc)(d) for d in range(len(antigen_list))))
                    status=True
                    print("run joblib.Parallel")
                except:
                    jb_res = list(jb.Parallel(n_jobs = -1, prefer = "threads")(jb.delayed(pfunc)(d) for d in range(len(antigen_list))))
                    status=True
                    print("run joblib.Parallel")
            
            if status:
                 for i in range(len(antigen_list)): ## "Wuhan-Hu-1" is always in each cross reactivity files produced by our pipeline
                     antigen = antigen_list[i]
                     EnvD_Min,EnvD_Max = jb_res[i]
                     VE["Proba Neut Min\n vs. %s antigen"%antigen] = EnvD_Min
                     VE["Proba Neut Max\n vs. %s antigen"%antigen] = EnvD_Max
            
            else:
               
                for i in range(len(antigen_list)): ## "Wuhan-Hu-1" is always in each cross reactivity files produced by our pipeline
                    antigen = antigen_list[i]
                    EnvD_Min,EnvD_Max = PNeut_Envelope(1, t_conc, [Lin_name], variants_in_cross, Cross_react_dic, c_dframe_dic, IC50xx_dic, antigen_list = antigen_list, mean_IC50xx = True)
                    VE["Proba Neut Min\n vs. %s antigen"%antigen] = EnvD_Min
                    VE["Proba Neut Max\n vs. %s antigen"%antigen] = EnvD_Max
                
            """ Save P_Neut ranges"""
            VE_df = pd.DataFrame(VE)
            VE_df.to_csv(sys.argv[w_save]+"/P_neut_"+Lin_name+".csv")
        else:
            pass
    
    
    try:
        if var_list_index is None:
            SpikeGroups_list_index = []
            for j in range(len(SpikeGroups_list)):
                if SpikeGroups_list[j] in variants_in_cross:
                    SpikeGroups_list_index.append(list(variants_in_cross).index(SpikeGroups_list[j]))
            SpikeGroups_list_index = np.array(SpikeGroups_list_index)
        else:
            SpikeGroups_list_index = var_list_index
        
        if spikegroups_proportion_adjust is None:
            if len(SpikeGroups_list_index)!=len(SpikeGroups_list):
                spikegroups_proportion_adjust = np.zeros((len(SpikeGroups_list_index), spikegroups_proportion.shape[1]))
                for j in range(len(SpikeGroups_list_index)):
                    w_j = list(SpikeGroups_list).index(variants_in_cross[SpikeGroups_list_index[j]])
                    spikegroups_proportion_adjust[j, :] = spikegroups_proportion[w_j, :]
                
                # renormalization
                spikegroups_proportion_adjust = spikegroups_proportion_adjust/np.sum(spikegroups_proportion_adjust, axis = 0)
            else:
                spikegroups_proportion_adjust = spikegroups_proportion.copy()

        else:
            spikegroups_proportion_adjust = spikegroups_proportion_adjust
            
        for key in c_dframe_dic.keys():
            PK_dframe = c_dframe_dic[key]
            key_num = np.array(re.findall(r"\d+", key)).astype(int)
            Res_sub_0 = Immunity_dynamics_fftconvolve(t, PK_dframe, infection_data = infection_data_corrected, 
                                                         present_variant_index = SpikeGroups_list_index, ### indexes of variant in variant_in_cross
                                                         tested_variant_list =  variant_to_sim, 
                                                         variant_name = variants_in_cross, ### Aligned with Cross_react_dic["variant_list"]
                                                         variant_proportion =  spikegroups_proportion_adjust, ### rows are aligned with present_variant_list that are present in cross
                                                         Ab_classes = Ab_classes, 
                                                         IC50xx= mean_IC50xx_dic,
                                                         Cross_react_dic = Cross_react_dic, 
                                                         )
            EI["t_half = %.3f \nt_max = %.3f"%(thalf_vec[key_num[0]], tmax_vec[key_num[1]])] = Res_sub_0
            Susc["t_half = %.3f \nt_max = %.3f"%(thalf_vec[key_num[0]], tmax_vec[key_num[1]])] = total_population - Res_sub_0
        """ Save Dynamics Without Vaccination """
        EI_df = pd.DataFrame(EI)
        EI_df.to_csv(sys.argv[w_save]+"/Immunized_SpikeGroup_%s_all_PK.csv"%variant_to_sim[0])
        
        Susc_df = pd.DataFrame(Susc)
        Susc_df.to_csv(sys.argv[w_save]+"/Susceptible_SpikeGroup_%s_all_PK.csv"%variant_to_sim[0])
        return "Done"
        
    except:
        
        return "Error"
    
w_save = 11 # index of resdir
save_pneut = str(sys.argv[12])
"""Load Lineage to assess """
Lin_name = str(sys.argv[13])
"""Update name if necessary -- this is the same update as in Compute_FR"""

if Lin_name != "ALL":
    ### Load spikegroups membership file
    file = open("Spikegroups_membership.pck", "rb")
    Pseudogroup_dic = pickle.load(file)
    file.close()
    if not run_group:
        num_antigen = int(sys.argv[14])
        k=15
        antigen_list = []
        while k<15+num_antigen:
            antigen = str(sys.argv[k])
            antigen_list.append(antigen)
            k+=1
        
        if Lin_name not in Pseudogroup_dic.keys():
            status_var = ei_util(Lin_name, variants_in_cross, antigen_list, Cross_react_dic, save_pneut=save_pneut, w_save = w_save) 
        else:
            try:
                file1 = open("results/Cross_react_dic_spikegroups_ALL.pck", "rb") # Check that this is the file you want to load
                Cross_react_dic = pickle.load(file1)
                variants_cross = list(Cross_react_dic["variant_list"])
                ### insert lin_sim in the position of it's Pseudogroup
                Cross_react_dic["variant_list"][variants_cross.index(Pseudogroup_dic[Lin_name])] = Lin_name
                file1.close()
                status_var = ei_util(Lin_name, variants_in_cross, antigen_list, Cross_react_dic, save_pneut=save_pneut, w_save = w_save) 
            except:
                try:
                    file1 = open("results/Cross_react_dic_spikegroups_present.pck", "rb") # Check that this is the file you want to load
                    Cross_react_dic = pickle.load(file1)
                    variants_cross = list(Cross_react_dic["variant_list"])
                    ### insert lin_sim in the position of it's Pseudogroup
                    Cross_react_dic["variant_list"][variants_cross.index(Pseudogroup_dic[Lin_name])] = Lin_name
                    file1.close()
                    status_var = ei_util(Lin_name, variants_in_cross, antigen_list, Cross_react_dic, save_pneut=save_pneut, w_save = w_save) 
                except:
                    pass
        # Save file as a placeholder for exectuted codes, required for snakemake
        sim_df = pd.DataFrame({"Lineage":[Lin_name], "Simulation status":[status_var]})
        sim_df.to_csv(sys.argv[w_save]+"/simulation_status_%s.csv"%Lin_name)
    else:
        Lin_list = []
        cross_list = []
        k = 13
        run_k = True
        not_pres = []
        while k<len(sys.argv) and run_k:
            lin_sim = str(sys.argv[k])                
            try:
                # will not give any error at num_antigen because num_antigen is an int
                num_antigen = int(sys.argv[k])
                loc_num_anti = k
                run_k = False
            except:
                run_k = True
            
            if run_k:
                try:
                    file1 = open(sys.argv[3]+"/Cross_%s.pck"%lin_sim, "rb") # Check that this is the file you want to load
                    cross_list.append(pickle.load(file1))
                    file1.close()
                    Lin_list.append(lin_sim)
                    k +=1
                except:
                    try:
                        try:
                            file1 = open("results/Cross_react_dic_spikegroups_%s.pck"%lin_sim, "rb") # Check that this is the file you want to load
                            Cross_dic = pickle.load(file1)
                            cross_list.append(Cross_dic)
                            file1.close()
                        except:
                            try:
                                file1 = open("results/Cross_react_dic_spikegroups_ALL.pck", "rb") # Check that this is the file you want to load
                                Cross_dic = pickle.load(file1)
                                variants_cross = list(Cross_dic["variant_list"])
                                file1.close()
                                if lin_sim in Pseudogroup_dic.keys():
                                    ### insert lin_sim in the position of it's Pseudogroup
                                    Cross_dic["variant_list"][variants_cross.index(Pseudogroup_dic[lin_sim])] = lin_sim
                                    cross_list.append(Cross_dic)
                                    file1.close()
                                    Lin_list.append(lin_sim)
                                else:
                                    not_pres.append(lin_sim)
                            except:
                                try:
                                    file1 = open("results/Cross_react_dic_spikegroups_present.pck", "rb") # Check that this is the file you want to load
                                    Cross_dic = pickle.load(file1)
                                    variants_cross = list(Cross_dic["variant_list"])
                                    file1.close()
                                    if lin_sim in Pseudogroup_dic.keys():
                                        ### insert lin_sim in the position of it's Pseudogroup
                                        Cross_dic["variant_list"][variants_cross.index(Pseudogroup_dic[lin_sim])] = lin_sim
                                        cross_list.append(Cross_dic)
                                        Lin_list.append(lin_sim)
                                    else:
                                        not_pres.append(lin_sim)   
                                except:
                                    pass
                        k +=1
                    except:
                        k +=1
        
        
        if len(not_pres)!=0:
            to_print = "Some Lineage in group do not have a cross reactivity file and are not present in covsonar data nor in results/ folder:" + " ".join(["%s"%var for var in not_pres])
            sys.exit(to_print)
            
        antigen_list = []
        k = loc_num_anti + 1
        while k<14+len(Lin_list)+num_antigen:
            antigen = str(sys.argv[k])
            antigen_list.append(antigen)
            k +=1
            
        for i in range(len(Lin_list)):
            lin_sim = Lin_list[i]
            Cross_react_dic = cross_list[i]
            variants_in_cross = Cross_react_dic["variant_list"]
            Cross_react_dic.pop("variant_list")
            status_var = ei_util(lin_sim, variants_in_cross, antigen_list, Cross_react_dic, save_pneut = save_pneut, w_save = w_save) 
            # Save file as a placeholder for exectuted codes, required for snakemake
            if i != 0:
                sim_df = pd.read_csv(sys.argv[w_save]+"/simulation_status_group.csv")
                Lin_List = sim_df["Lineage"].tolist() + [lin_sim]
                status_var_list = sim_df["Simulation status"].tolist() + [status_var]
                sim_df = pd.DataFrame({"Lineage":Lin_List, "Simulation status":status_var_list})
            else:    
                sim_df = pd.DataFrame({"Lineage":[lin_sim], "Simulation status":[status_var]})
                
            sim_df.to_csv(sys.argv[w_save]+"/simulation_status_group.csv")    

else:
    status_var = []
    SpikeGroups_list_index = []
    
    num_antigen = int(sys.argv[14])
    k=15
    antigen_list = []
    while k<15+num_antigen:
        antigen = str(sys.argv[k])
        antigen_list.append(antigen)
        k+=1
    
    for j in range(len(SpikeGroups_list)):
        if SpikeGroups_list[j] in variants_in_cross:
            SpikeGroups_list_index.append(list(variants_in_cross).index(SpikeGroups_list[j]))       
    SpikeGroups_list_index = np.array(SpikeGroups_list_index)
    add_print = False
    if len(SpikeGroups_list_index)!=len(SpikeGroups_list):
        try: ### use updated cross ALL including missing if it was computed
            file1 = open("results/Cross_react_dic_spikegroups_present.pck", "rb") 
            Cross_react_dic = pickle.load(file1)
            variants_in_cross = Cross_react_dic["variant_list"]
            Cross_react_dic.pop("variant_list")
            file1.close()
            spikegroups_proportion_adjust = spikegroups_proportion.copy()
            # regenerage the indexes
            SpikeGroups_list_index = []
            for j in range(len(SpikeGroups_list)):
                if SpikeGroups_list[j] in variants_in_cross:
                    SpikeGroups_list_index.append(list(variants_in_cross).index(SpikeGroups_list[j]))       
            SpikeGroups_list_index = np.array(SpikeGroups_list_index)
            
        except:
            # readjust variant proportions to spikegroups available in cross react
            add_print = True
            spikegroups_proportion_adjust = np.zeros((len(SpikeGroups_list_index), spikegroups_proportion.shape[1]))
            for j in range(len(SpikeGroups_list_index)):
                w_j = list(SpikeGroups_list).index(variants_in_cross[SpikeGroups_list_index[j]])
                spikegroups_proportion_adjust[j, :] = spikegroups_proportion[w_j, :]
            
            # renormalization
            spikegroups_proportion_adjust = spikegroups_proportion_adjust/np.sum(spikegroups_proportion_adjust, axis = 0)
    else:
        spikegroups_proportion_adjust = spikegroups_proportion.copy()
    
    for i in range(len(SpikeGroups_list)):
        if SpikeGroups_list[i] in variants_in_cross:
            if add_print:
                print("A smaller set of spikesgroups are being simulated for all_il = TRUE \n Make sure this is what you want otherwise first set the parameter cross_missing to TRUE")
                
            print("Compute E[immunized] for %s (%d out of %d spikegroups + Wuhan-Hu-1)"%(SpikeGroups_list[i], i, len(SpikeGroups_list)-1))
            status_var.append(ei_util(SpikeGroups_list[i], 
                                      variants_in_cross = variants_in_cross,
                                      antigen_list = antigen_list,
                                      Cross_react_dic = Cross_react_dic,
                                      save_pneut = save_pneut, 
                                      var_list_index=SpikeGroups_list_index, 
                                      spikegroups_proportion_adjust=spikegroups_proportion_adjust, 
                                      w_save = w_save))
        else:
            status_var.append("Not in cross")
    # Save file as a placeholder for exectuted codes, required for snakemake
    sim_df = pd.DataFrame({"SpikeGroups":SpikeGroups_list, "Simulation status":status_var})
    sim_df.to_csv(sys.argv[w_save]+"/simulation_status_ALL.csv")
    