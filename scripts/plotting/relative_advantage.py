import numpy as np
import pandas as pd
import numpy.ma as ma
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def moving_average(X, window = 7):
    u = np.zeros(len(X))
    u[:window] = X[:window]
    for i in range(window, len(X)):
        u[i] = np.mean(X[i-window:i+1])
    
    return u

ES_df = pd.read_csv(sys.argv[1])
lineage_freq = pd.read_csv(sys.argv[2])
threshold = sys.argv[3]
 variant = sys.arg[4]

#ES_df = pd.read_csv("demo/results/Immunological_Landscape/Susceptible_SpikeGroup_lineage_XXX_all_PK.csv")
#lineage_freq = pd.read_csv("demo/results/Daily_Lineages_Freq.csv")
#threshold = 0 #(percent)
#variant = "BA.5.1"

# processing of susceptibles 
ES_df.drop(columns = "Unnamed: 0", inplace = True)
es_cols = ES_df.columns
ES_ranges = ES_df.to_numpy()[:, es_cols!="Days"].astype(float)
    
t_dates = ES_df["Days"]

# processing of frequency data
lineage_freq.drop(columns = "Unnamed: 0", inplace = True)
lineage_freq = lineage_freq[lineage_freq['date'].isin(t_dates)]
freqs = lineage_freq.loc[:, lineage_freq.columns != 'date']

# imputing frequencies below threshold and normalization
freqs = freqs.mask(freqs < threshold)
freqs = freqs.fillna(0)
col_sums = freqs.sum(axis = 1).values
freqs = freqs.divide(col_sums, axis="rows")
lineage_freq.loc[:, lineage_freq.columns != 'date'] = freqs

# processing mean
# needs to be updated to allow individual weighting 
S_mean_df = pd.read_csv("/Users/nilsgubela/Desktop/VASIL NO GIT/demo/results/Immunological_Landscape/Susceptible_weighted_mean_over_pseudogroups_all_PK.csv")
S_mean_df = S_mean_df[S_mean_df['date'].isin(t_dates)]
S_all_mean = S_mean_df.to_numpy()[:, S_mean_df.columns != "date"].astype(float)

# calculation of change in relative frequency from model
gamma_SI = np.zeros((len(t_dates), ES_ranges.shape[1]))

for i in range(ES_ranges.shape[1]):
    S_x = ES_ranges[:, i]
    S_mean = S_all_mean[:, i]

    gamma_SI[:, i] = np.divide(S_x - S_mean, S_mean, out = S_x, where = S_mean != 0)

   
# get min max gamma over PK at each timepoints
gamma_SI_min, gamma_SI_max = np.min(gamma_SI, axis = 1), np.max(gamma_SI, axis = 1)

# change in relative frequency from genomic surveillance data 
gamma_prop = np.diff(np.log(lineage_freq[variant]))

# Smooth gamma_prop with moving_average
window = 14
gamma_prop = moving_average(gamma_prop, window = window)
    

  
# mask absent variant prop data
Prop = moving_average(lineage_freq[variant], window = window)
gamma_SI_min = ma.array(gamma_SI_min, mask = Prop<0) ## filter out by prop threshold
gamma_SI_max = ma.array(gamma_SI_max, mask = Prop<0) ## filter out by prop threshold

gSI_min = ma.array([gamma_SI_min], mask = [gamma_SI_min.mask])
gSI_max = ma.array([gamma_SI_max], mask = [gamma_SI_max.mask])


# plotting
fig, ax = plt.subplots()

plt.fill_between(t_dates, gamma_SI_min, gamma_SI_max, color = "green")
plt.plot(t_dates[:-1], gamma_prop, color = "orange")
#ax.axhline(xmin = 0, xmax = t_dates[-1], ls = "--", linewidth = 2, color = "black")
ax.axhline(xmin = 0, xmax = len(t_dates), ls = "--", linewidth = 2, color = "black")

perday = range(0,len(t_dates), 28)

ax.set_xticks(perday)
ax.set_xticklabels(t_dates[perday].tolist(),
    rotation = 45, horizontalalignment = "right")
    
pdf = PdfPages("realtive_fitness_%s.pdf"%variant)
pdf.savefig(fig, bbox_inches = "tight")
pdf.close()
fig.savefig("realtive_fitness_%s.svg"%variant, bbox_inches = "tight")

