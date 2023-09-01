import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from patsy import cr
from sklearn.linear_model import LinearRegression
from scipy.interpolate import CubicSpline

# inputs

infection_data_corrected = sys.argv[1]
S_mean_df = sys.argv[2]
color = sys.argv[3]

#infection_data_corrected = pd.read_csv("demo/caseAscertainmentTable_reportedCasesRatio.csv")
#S_mean_df = pd.read_csv("demo/results/Immunological_Landscape/Susceptible_weighted_mean_over_pseudogroups_all_PK.csv")
#color = "red"


t_dates = S_mean_df["date"]
infection_data_corrected  = infection_data_corrected [infection_data_corrected ['date'].isin(t_dates)]


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
    if(x[i] > r[r_index]):
        intervals.append(interval)
        interval = [x[i]]
        r_index += 1

intervals.append(interval)
intervals.pop(0)


# Fig presettings
fig, ax = plt.subplots()


# I plot
plt.plot(x, infection_data_corrected["minNTrue"], color = "black", label = "Infections", linewidth = 3)
plt.ylabel("Infections", fontsize = 25)

plt.legend(loc = (0.4, 0.8), fontsize = 25)
 
perday = range(0,len(t_dates), 28)

ax.set_xticks(perday)
ax.set_xticklabels(t_dates[perday].tolist(),
    rotation = 45, horizontalalignment = "right")
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
plt.legend(loc = (0.7, 0.8) ,fontsize = 25)

# add vlines
plt.vlines(x=r.real[abs(r.imag)<1e-5], ymin=min(r_ABS_min), ymax=max(r_ABS_max), colors="black", ls='--', lw=0.5, label='Delta infection')




plt.suptitle("Absolute growth Simulations Without Vaccination Timeline", fontsize = 20)
plt.subplots_adjust(hspace=0.75, wspace=0.25)
pdf = PdfPages("absolute_estimate.pdf")
pdf.savefig(fig, bbox_inches = "tight")
pdf.close()
plt.savefig("absolute_estimate.svg")
