import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from patsy import cr
from sklearn.linear_model import LinearRegression
from scipy.interpolate import CubicSpline
import sys
import pdb

# inputs
infection_data_corrected = pd.read_csv(sys.argv[1])
S_mean_file = sys.argv[2]
S_mean_df = pd.read_csv(S_mean_file)
color = sys.argv[3]

#infection_data_corrected = pd.read_csv("demo/caseAscertainmentTable_reportedCasesRatio.csv")
#S_mean_df = pd.read_csv("demo/results/Immunological_Landscape/Susceptible_weighted_mean_over_pseudogroups_all_PK.csv")
#color = "red"
t_dates = S_mean_df["Days"]
infection_data_corrected  = infection_data_corrected[infection_data_corrected['date'].isin(t_dates)]

try:
    S_mean_df.drop(columns = "Unnamed: 0", inplace = True)
except:
    pass

# fit cubic splines to model predictions 
r_ABS_min, r_ABS_max = np.min(S_mean_df, axis = 1), np.max(S_mean_df, axis = 1)

df = 20
x = range(0,len(t_dates))
y = [np.median([r_ABS_min[i],r_ABS_max[i]]) for i in range(len(r_ABS_min))]
x_basis = cr(x, df=df, constraints="center")

# Fit model to the data
model = LinearRegression().fit(x_basis, y)

# Get estimates
y_hat = model.predict(x_basis)

# fit spline to prediction
spl = CubicSpline(x, y_hat)
dspl = spl.derivative()
ddspl = dspl.derivative()

# create intervals between change points of model prediction curve
r = ddspl.roots()
intervals = []
interval = []
r_index = 0
for i in range(len(x)):
    interval.append(x[i])
    if len(r)>r_index:
        if (x[i] > r[r_index]):
            intervals.append(interval)
            interval = [x[i]]
            r_index += 1

intervals.append(interval)
intervals.pop(0)


# Fig presettings
import matplotlib
def PreFig(xsize = 12, ysize = 12):
    '''
    @brief: customize figure parameters
    '''
    matplotlib.rc('xtick', labelsize=xsize) 
    matplotlib.rc('ytick', labelsize=ysize)
PreFig(xsize = 20, ysize = 20)
fig = plt.figure(figsize = (15, 7))
ax = fig.add_subplot(1, 1, 1)

# I plot
cases_col = str(sys.argv[len(sys.argv) - 1])
plt.plot(x, infection_data_corrected[cases_col], color = "black", label = "Infections", linewidth = 3)
plt.ylabel("Infections", fontsize = 25)

# S plot 
ax2 = ax.twinx()
for i in range(len(intervals)):
    t_int = intervals[i]
    y_int = [spl(i) for i in range(t_int[0]-1, t_int[-1])]
    plt.plot(t_int, y_int, linewidth = 3, color = color)
    if color == "green":
        color = "red"
    else:
        color = "green"
ax2.fill_between(x, r_ABS_min, r_ABS_max, alpha = 0.3, color = "grey", label = "$S$ ranges")


plt.ylabel("Susceptibles", fontsize = 25)

ax.legend(loc='upper center', bbox_to_anchor=(0.6, -0.2),
          fancybox=True, shadow=True, ncol=5)
ax2.legend(loc='upper center', bbox_to_anchor=(0.3, -0.2),
          fancybox=True, shadow=True, ncol=5)

try:
    x_min = list(t_dates).index(str(sys.argv[5]))
    if str(sys.argv[6]) in t_dates:
        x_max = list(t_dates).index(str(sys.argv[6]))
    else:
        x_max = len(t_dates) - 1
except:
    x_min = x[0]
    x_max = x[-1]

if (x_min is not None):
    ax.set_xlim((x_min, x_max))
    ax2.set_xlim((x_min, x_max))    
    t_dates_show = np.array(t_dates)[x_min:x_max+1]
else:
    t_dates_show = t_dates

if len(t_dates_show)>200:
    pp = 7*4
else:
    pp = min(len(t_dates_show), 14)

perday = np.arange(0,len(t_dates_show), pp)
date_ticks = t_dates_show[perday].tolist()
if t_dates[len(t_dates) - 1] not in date_ticks:
    n=list(t_dates).index(date_ticks[-1])+pp
    while n<len(t_dates)-1:
        date_ticks.append(t_dates[n])
        perday = np.append(perday, n)
        n += pp
    date_ticks.append(t_dates[len(t_dates) - 1])
    perday = np.append(perday, len(t_dates) - 1)

if x_min is not None:
    perday_orig = []
    for i in range(len(date_ticks)):
        perday_orig.append(list(t_dates).index(date_ticks[i]))
else:
    perday_orig = perday
    
ax.set_xticks(perday_orig)
ax.set_xticklabels(date_ticks,
    rotation = 45, horizontalalignment = "right")

ax2.set_xticks(perday_orig)
ax2.set_xticklabels(date_ticks,
    rotation = 45, horizontalalignment = "right")

# add vlines
plt.vlines(x=r.real[abs(r.imag)<1e-5], ymin=min(r_ABS_min), ymax=max(r_ABS_max), colors="black", ls='--', lw=0.5, label='Delta infection')

plt.suptitle("Absolute growth Simulations", fontsize = 20)
plt.subplots_adjust(hspace=0.75, wspace=0.25)
pdf = PdfPages(sys.argv[4]+"/absolute_estimate.pdf")
pdf.savefig(fig, bbox_inches = "tight")
pdf.close()
plt.savefig(sys.argv[4]+"/absolute_estimate.svg", bbox_inches='tight')


# Fig presettings
PreFig(xsize = 20, ysize = 20)
fig4d = plt.figure(figsize = (14, 7))
ax = fig4d.add_subplot(1, 1, 1)
x = np.arange(0,len(t_dates)).astype(int)

# f' plot
plt.plot(x, dspl(x), color = "turquoise", label = r'$\dfrac{d}{dt} S(t)$', linewidth = 3)
ax.axhline(xmin = 0, xmax = x[-1], ls = "--", linewidth = 4, color = "black") 
plt.legend(loc = (0.46, 0.8), fontsize = 25)
ax.set_xticks([x[x_min]]+list(x[x_min+pp:x_max][::pp])+[x[x_max]])
ax.set_xticklabels([t_dates[x_min]]+list(t_dates[x_min+pp:x_max][::pp])+[t_dates[x_max]], 
                   rotation = 45, horizontalalignment = "right")

# f'' plot
ax2 = ax.twinx()
ax2.set_ylim(-abs(max(ddspl(x))),abs(max(ddspl(x))))
plt.plot(x, ddspl(x), linewidth = 3, color = "orange", label = r'$\left(\dfrac{d}{dt}\right)^2 S(t)$')
ax2.set_xlim((x[x_min], x[x_max]))

ax.legend(loc = (1.2, 0.5) ,fontsize = 25)
ax2.legend(loc = (1.2, 0.) ,fontsize = 25)
# add vlines
plt.vlines(x=r.real[abs(r.imag)<1e-5], ymin=min(ddspl(x)), ymax=max(ddspl(x)), colors="black", ls='--', lw=0.5, label='Delta infection')

#plt.suptitle("First and Second Derivative of estimated Susceptibles", fontsize = 20)
#plt.subplots_adjust(hspace=0.75, wspace=0.25)
pdf = PdfPages(sys.argv[4]+"/derivative_trends.pdf")
pdf.savefig(fig4d, bbox_inches = "tight")
pdf.close()
plt.savefig(sys.argv[4]+"/derivative_trends.svg")
