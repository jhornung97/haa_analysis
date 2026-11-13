import numpy as np
import matplotlib.pyplot as plt
import os 
import sys
import mplhep
from scipy.optimize import curve_fit

mplhep.style.use("CMS")

year = sys.argv[1] 
scope = sys.argv[2]

if year == "2016preVFP":
    lumi = 19.5
elif year == "2016postVFP":
    lumi = 16.8
elif year == "2017":
    lumi = 41.5
elif year == "2018":
    lumi = 59.8
elif year == "2024":
    lumi = 109.95

data = np.loadtxt(f"/work/jhornung/Haa/alphas_{year}_{scope}.csv", delimiter=" ")
sorted_idx = np.argsort(data[:, 0])
iso_cuts = data[sorted_idx, 0]
alphas = data[sorted_idx, 1]
alphas_err = data[sorted_idx, 2]

mean_alpha = np.mean(alphas)
mean_alpha_err = np.std(alphas) 

fig, axs = plt.subplots()
axs.axhline(mean_alpha, color='red', linestyle='--', label=f'Mean')
iso_cuts_fill = np.insert(iso_cuts, 0, 0)
iso_cuts_fill = np.append(iso_cuts_fill, np.max(iso_cuts) + 0.5)
iso_cuts_fill = np.sort(iso_cuts_fill)
axs.fill_between(iso_cuts_fill,
                mean_alpha - mean_alpha_err, 
                mean_alpha + mean_alpha_err, 
                color='gray', alpha=0.5, linewidth=0, label=r'1$\sigma$ Band')
axs.fill_between(iso_cuts_fill,
                mean_alpha - 2 * mean_alpha_err, 
                mean_alpha + 2 * mean_alpha_err,
                color='gray', alpha=0.2, linewidth=0, label=r'2$\sigma$ Band')
axs.errorbar(iso_cuts, alphas, yerr=alphas_err, fmt='o', label=r'Fitted $\alpha$', color='black')

axs.set_xlabel('Isolation Cut')
axs.set_xlim(0, np.max(iso_cuts) + 0.5)
axs.set_ylabel(r'Exponential fit $\alpha$')
plt.legend()

axs.set_title(r"$\it{" + "Private \: work" + "}$" +" "+ r"$\bf{" + "(CMS \: simulation)" + "}$",loc="left", fontsize=20)
if int(year[:4]) <= 2018:
    axs.set_title(f"{lumi} fb$^{{-1}}$, {year} (13 $\mathrm{{{{TeV}}}}$)",loc="right", fontsize=20)
elif int(year[:4]) > 2018:
    axs.set_title(f"{lumi:.1f} fb$^{{-1}}$, {year} (13.6 $\mathrm{{{{TeV}}}}$)",loc="right", fontsize=20)
plt.tight_layout()

plt.savefig(f"/web/jhornung/public_html/analysis_plots/{year}/alphas_{scope}.png")
plt.savefig(f"/web/jhornung/public_html/analysis_plots/{year}/alphas_{scope}.pdf")

plt.show()