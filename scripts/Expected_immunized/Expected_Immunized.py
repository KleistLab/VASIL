#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import pickle
import re
import joblib as jb
import sys
import warnings


"""Load Infection Data"""
Population_Data = pd.read_csv(sys.argv[1])

"""Load relevant cross_neutralization files"""
file1 = open(sys.argv[2], "rb") # premade simulations
Cross_with_delta_validation = pickle.load(file1)
variants_x_names_show = Cross_with_delta_validation["variant_list"]
Cross_with_delta_validation.pop("variant_list")
file1.close()

file1 = open(sys.argv[3], "rb") # Check that this is the file you want to load
Cross_react_dic = pickle.load(file1)
variants_in_cross = Cross_react_dic["variant_list"]
file1.close()

Ab_classes = list(Cross_react_dic.keys())

"""Spike groups and frequencies"""
file1 = open(sys.argv[4], "rb") 
SpikeGroups_list = pickle.load(file1)["names"]
file1.close()
frequency_spk_df = pd.read_csv(sys.argv[5])

"""Escape fraction data """
Escape_Fraction = pd.read_csv(sys.argv[6])

"""Load vaccing efficacy data for fitting """
VE_Delta_df = pd.read_csv(sys.argv[7])

"""Load total population and re-adjust infection timeline"""
total_population = float(sys.argv[8])

date_start = str(sys.argv[9])
Population_Data = Population_Data.drop(index = Population_Data.index[:list(Population_Data['date']).index(date_start)])

date_end = sys.argv[10]
where_last_day = list(Population_Data['date']).index(date_end) + 1

t = np.arange(1, len(Population_Data['date'][:where_last_day])+1, 1) # array of timepoints at which to compute the antibody concentration
infection_data_corrected = Population_Data['minNTrue'].values[:where_last_day]
t_dates = Population_Data['date'].values[:where_last_day]
days_incidence = list(Population_Data['date'][:where_last_day]) 

"""Load Lineage to assess """
Lin_name = sys.argv[11]
"""Update name if necessary -- this is the same update as in Compute_FR"""
if Lin_name in SpikeGroups_list: 
    Lin_name = Lin_name + "_requested" ### renamed to avoid ambiguities
        
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

Ab_classes = list(Cross_with_delta_validation.keys())
"""Load fold change IC50 in data"""
FC_ic50_dic = {Ab_classes[i]:(mean_IC50_per_group["fold_change"].values[mean_IC50_per_group["group"] == Ab_classes[i]])[0] for i in range(len(Ab_classes))}

try:
	VE_Delta_df.drop(columns = "Unnamed: 0", inplace = True)
except:
	pass
	
ve_fitting = VE_Delta_df["Processed VE (used in IC50 fitting)"].values
days_fitting = VE_Delta_df["Day since vacc (used in IC50 fitting)"].values

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
    IC50_data = root(sqrt_diff_FR, guess, args = ([days_fitting], np.array(FC_ic50_list)*np.array(FR_delta), [ve_fitting], len(Ab_classes), c_dframe), method = "lm").x
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
def Immunity_dynamics_fftconvolve(t, PK_dframe, infection_data, present_variant_list, tested_variant_list, variant_name, variant_proportion, Ab_classes, 
                                  IC50xx, Cross_react_dic, escape_per_sites = None, mut_sites_per_variant = None, parallel = False, mode_func = None):
    
    Infected_l_vect = infection_data[np.newaxis, :]*variant_proportion[:, :len(infection_data)]
    present_variant_index = np.array([list(variant_name).index(present_variant_list[j]) for j in range(len(present_variant_list))])

    Prob_Neut = P_Neut(t, present_variant_index, PK_dframe, tested_variant_list, variant_name, Ab_classes, IC50xx, Cross_react_dic)
    
    """
    Conv_Mat = np.zeros(Prob_Neut.shape)
    for i in range(Conv_Mat.shape[0]):
        Conv_Mat[i, :] = signal.fftconvolve(Infected_l_vect[i, :], Prob_Neut[i])[:len(t)]
    """   
    
    Conv_Mat = signal.fftconvolve(Infected_l_vect, Prob_Neut, axes = 1)[:, :len(t)]
    
    # No normalization
    Expected_Immuned = np.sum(Conv_Mat, axis = 0)
    """
    tested that this gives the as Immunity_dynamics and is 200x faster
    """
    return Expected_Immuned


"""Compute Antibody concentration over time for a range of t_half and t_max"""
thalf_vec = np.linspace(25, 69, 15) 
tmax_vec = np.linspace(14, 28, 5)
t_conc = np.arange(1, 700, 1)
c_t_vec, c_dframe_dic, dataname = Antibody_ranges(thalf_vec, tmax_vec, t_conc, Ab_classes)
IC50xx_dic, mean_IC50xx_dic = Find_IC50_ranges(thalf_vec, tmax_vec, t_conc, Ab_classes,  Cross_with_delta_validation)
### spikegroup frequency
print("Fitted mean IC50=%.3f \n IC50 per epitope class is "%mean_IC50xx_dic["NTD"], mean_IC50xx_dic)
try:
	frequency_spk_df.drop(columns = "Unnamed: 0", inplace = True)
except:
	pass

### Make sure proportions axis 0 is aligned with Spikegroup_list
spikegroups_freq = np.zeros((len(SpikeGroups_list), len(t)))
for i in range(len(SpikeGroups_list)):
    if SpikeGroups_list[i]!="Wuhan-Hu-1":
        spikegroups_freq[i, :] = frequency_spk_df["Spike. "+SpikeGroups_list[i]][:len(t)]

NormProp = np.sum(spikegroups_freq, axis = 0)
prop_rounded = np.round(spikegroups_freq,decimals = 10)
spikegroups_proportion = prop_rounded/NormProp

### end of simulation
def ei_util(Lin_name):
    variant_to_sim = [Lin_name]
    EI = {}
    EI["Days"] = days_incidence
    
    Susc = {}
    Susc["Days"] = days_incidence
    try:
        for key in c_dframe_dic.keys():
            PK_dframe = c_dframe_dic[key]
            key_num = np.array(re.findall(r"\d+", key)).astype(int)
    
            Res_sub_0 = Immunity_dynamics_fftconvolve(t, PK_dframe, infection_data = infection_data_corrected, 
                                                         present_variant_list = SpikeGroups_list, ### Aligned with rows-indexes of variant_proportion
                                                         tested_variant_list =  variant_to_sim, 
                                                         variant_name = variants_in_cross, ### Aligned with Cross_react_dic["variant_list"]
                                                         variant_proportion =  spikegroups_proportion, ### rows are aligned with present_variant_list
                                                         Ab_classes = Ab_classes, 
                                                         IC50xx= mean_IC50xx_dic,
                                                         Cross_react_dic = Cross_react_dic, 
                                                         )
            
            EI["t_half = %.3f \nt_max = %.3f"%(thalf_vec[key_num[0]], tmax_vec[key_num[1]])] = Res_sub_0
            Susc["t_half = %.3f \nt_max = %.3f"%(thalf_vec[key_num[0]], tmax_vec[key_num[1]])] = total_population - Res_sub_0
            
        """ Save Dynamics Without Vaccination """
        EI_df = pd.DataFrame(EI)
        EI_df.to_csv(sys.argv[12]+"/Immunized_SpikeGroup_%s_all_PK.csv"%variant_to_sim[0])
        
        Susc_df = pd.DataFrame(Susc)
        Susc_df.to_csv(sys.argv[12]+"/Susceptible_SpikeGroup_%s_all_PK.csv"%variant_to_sim[0])
        return "Done"
    
    except:
        
        return "Error"
 
       
simulated_var = ei_util(Lin_name)   
# Save file as a placeholder for exectuted codes, require for snakemake
sim_df = pd.DataFrame({"SpikeGroups":[Lin_name], "Simulation status":simulated_var})
sim_df.to_csv(sys.argv[12]+"/simulation_status.csv")