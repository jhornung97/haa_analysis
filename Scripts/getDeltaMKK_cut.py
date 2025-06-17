import uproot
import numpy as np
import awkward as ak
import matplotlib.pyplot as plt
import mplhep
from scipy.optimize import curve_fit
import sys
import utility as util

def gauss(x, a, x0, sigma):
    return a * np.exp(-1.0 * (x - x0) ** 2 / (2 * sigma ** 2))   

def gauss_err(x, a, x0, sigma, a_err, x0_err, sigma_err):
    exp_term = np.exp(-1.0 * (x - x0) ** 2 / (2 * sigma ** 2))
    d_a = exp_term
    d_x0 = a * exp_term * (x - x0) / sigma ** 2
    d_sigma = a * exp_term * (x - x0) ** 2 / sigma ** 3
    return np.sqrt((d_a * a_err) ** 2 + (d_x0 * x0_err) ** 2 + (d_sigma * sigma_err) ** 2) 

ipath = sys.argv[1]

with uproot.open(ipath) as f:
    tree = f["ntuple"]
    branches = tree.arrays()

H_mass_mask = (branches["H_mass"] > 70) & (branches["H_mass"] < 200)
H_eta_mask = np.abs(branches["H_eta"]) < 2.4
leading_k_mask = (branches["d1_pt"] > 10) & (branches["d2_pt"] > 10)
z_pt_mask = branches["pt_vis"] > 30
z_m_mask = (branches["m_vis"] > 75) & (branches["m_vis"] < 105)
#iso_mask = (branches["d1_iso"] < 3) & (branches["d2_iso"] < 3) & (branches["d3_iso"] < 3) & (branches["d4_iso"] < 3)
ps_mask = (branches["ps_1_mass"] < 3) & (branches["ps_2_mass"] < 3)
mask = H_mass_mask & H_eta_mask & leading_k_mask & z_pt_mask & z_m_mask & ps_mask

sfs = 0
scope = None

if "/mm/" in ipath:
    scope = "mm"
    sfs = branches["id_wgt_mu_1"][mask] * branches["id_wgt_mu_2"][mask] * branches["iso_wgt_mu_1"][mask] * branches["iso_wgt_mu_2"][mask] * branches["trigger_wgt_mu_1"][mask] * branches["trigger_wgt_mu_2"][mask]
elif "/ee/" in ipath:
    scope = "ee"
    sfs = branches["id_wgt_ele_1"][mask] * branches["id_wgt_ele_2"][mask] * branches["trigger_wgt_ele_1"][mask] * branches["trigger_wgt_ele_2"][mask]

signs = np.array(list(map(util.weights, branches["genWeight"][mask])))
weights = 59.8e3 * sfs * signs * branches['evtweight'][mask] * branches['puweight'][mask] 

plt.style.use(mplhep.style.CMS)
fig, ax = plt.subplots()

plt.hist(branches['ps_1_mass'][mask] - branches['ps_2_mass'][mask], bins=60, range=(-0.15, 0.15), weights=weights, histtype='step', label=r'$Z(\ell\ell)H(aa\rightarrow 4K)$', color='red')

dmkk = branches['ps_1_mass'][mask].to_numpy() - branches['ps_2_mass'][mask].to_numpy()
w = weights.to_numpy()

hist, bins = np.histogram(dmkk, bins=60, range=(-0.15, 0.15), weights=w)
bin_centers = 0.5 * (bins[1:] + bins[:-1])

hist_err = np.sqrt(hist)
plt.errorbar(bin_centers, hist, yerr=hist_err, fmt='none', color='red')

fit_xmin, xmax = -0.05, 0.05
fit_mask = (bin_centers > fit_xmin) & (bin_centers < xmax)
fit_x = bin_centers[fit_mask]
fit_y = hist[fit_mask]

popt, pcov = curve_fit(gauss, xdata=fit_x, ydata=fit_y, p0=[np.sum(fit_y), 0, 0.02])

fit_y_exp = gauss(fit_x, *popt)

fit_err = hist_err[fit_mask]

chi2 = np.sum(((fit_y - fit_y_exp) / fit_err) ** 2)
ndof = len(fit_x) - len(popt)

#perr = np.sqrt(np.diag(pcov))

x = np.linspace(-0.05, 0.05, 1000)
plt.plot(x, gauss(x, *popt), color='black', label='Gaussian fit', linewidth=2)
#plt.fill_between(x, gauss(x, *popt) - gauss_err(x, *popt, *perr), gauss(x, *popt) + gauss_err(x, *popt, *perr), alpha=0.3, edgecolor='white', facecolor='grey')
plt.axvline(x=-1, color='white', linestyle='--', linewidth=1, label=fr'norm = ${popt[0]:.0f}\pm {np.sqrt(pcov[0][0]):.0f}$')
plt.axvline(x=-1, color='white', linestyle='--', linewidth=1, label=fr'$\mu$ = ${popt[1]:.4f}\pm {np.sqrt(pcov[1][1]):.4f}$')
plt.axvline(x=-1, color='white', linestyle='--', linewidth=1, label=fr'$\sigma$ = ${popt[2]:.4f}\pm {np.sqrt(pcov[2][2]):.4f}$')
plt.axvline(x=-1, color='white', linestyle='--', linewidth=1, label=f'$\chi^2/dof$ = {chi2/ndof:.2f}')
plt.legend(loc='upper right', fontsize=15)

plt.xlabel(r'$\Delta m_{K^+K^-}$ [GeV]')
plt.ylabel('Events')

plt.xlim(-0.15, 0.15)

ax.set_title(r"$\it{" + "Private \: work" + "}$" +" "+ r"$\bf{" + "(CMS \: simulation)" + "}$",loc="left", fontsize=20)
ax.set_title(f"59.8 fb$^{{-1}}$, 2018 (13 $\mathrm{{{{TeV}}}}$)",loc="right", fontsize=20)

#plt.savefig(f'/web/jhornung/public_html/analysis_plots/2018/deltaMKK_cut_{scope}.png')
#plt.savefig(f'/web/jhornung/public_html/analysis_plots/2018/deltaMKK_cut_{scope}.pdf')

plt.show()