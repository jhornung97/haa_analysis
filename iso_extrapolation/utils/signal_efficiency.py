import uproot
from concurrent.futures import ThreadPoolExecutor
import awkward as ak
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import crystalball, norm
import mplhep as hep
import sys
import csv
import os

sys.path.append('/work/jhornung/Haa/Scripts/')
import HistLib as hl
import utility as util

def gaussian(x, a, x0, sigma):
    return a * norm.pdf(x, loc=x0, scale=sigma)

def crystalball_func(x, a, beta, m, x0, scale):
    return a * crystalball.pdf(x, beta, m, loc=x0, scale=scale)

def combined_model(x, a, frac_gauss, x0, sigma, beta, m, scale):
    frac_cb = 1 - frac_gauss
    return a * (frac_gauss * norm.pdf(x, x0, sigma) + frac_cb * crystalball.pdf(x, beta, m, x0, scale))

def moving_average(data, window_size):
    return np.convolve(data, np.ones(window_size)/window_size, mode='same')

plt.style.use(hep.style.CMS)

scope = sys.argv[1]

signal_loc = "/ceph/jhornung/MC_2018/2018/"
signal = "HaaKKKK_MH125_ma2.0_13TeV_powheg_pythia8-RunIISummer20UL18MiniAODv2-106X"

data = util.read_file(signal_loc + signal + f"/{scope}/{signal}.root", data_type="mc")

H_mass_mask = (data["H_mass"] > 70) & (data["H_mass"] < 200)
H_eta_mask = (np.abs(data["H_eta"]) < 2.4)
leading_k_mask = (data["d1_pt"] > 10) & (data["d2_pt"] > 10)
ps_mass_mask = (data["ps_1_mass"] < 3) & (data["ps_2_mass"] < 3) & (np.abs(data["ps_1_mass"] - data["ps_2_mass"]) < 0.06)
pt_vis_mask = (data["pt_vis"] > 30)
m_vis_mask = (data["m_vis"] > 75) & (data["m_vis"] < 105)
mask = H_mass_mask & H_eta_mask & leading_k_mask & ps_mass_mask & pt_vis_mask & m_vis_mask

if scope == "mm":
    sfs = data["id_wgt_mu_1"] * data["id_wgt_mu_2"] * data["iso_wgt_mu_1"] * data["iso_wgt_mu_2"] * data["trigger_wgt_mu_1"] * data["trigger_wgt_mu_2"]
elif scope == "ee":
    sfs = data["id_wgt_ele_1"] * data["id_wgt_ele_2"] * data["trigger_wgt_ele_1"] * data["trigger_wgt_ele_2"]
signs = np.array(list(map(util.weights, data["genWeight"])))
mc_weights = 59.8e3 * sfs * signs * data["evtweight"] * data["puweight"]

reference_hist, _ = np.histogram(data["H_mass"][mask], bins=65, range=(70, 200), weights=mc_weights[mask])
reference_integral = np.sum(reference_hist)

reference_fit_hist, reference_bin_edges = np.histogram(data["H_mass"][mask & ((data["d1_iso"] < 0.5) & (data["d2_iso"] < 0.5))], bins=40, range=(115, 135), weights=mc_weights[mask & ((data["d1_iso"] < 0.5) & (data["d2_iso"] < 0.5))])
reference_sumw2, _ = np.histogram(data["H_mass"][mask & ((data["d1_iso"] < 0.5) & (data["d2_iso"] < 0.5))], bins=40, range=(115, 135), weights=mc_weights[mask & ((data["d1_iso"] < 0.5) & (data["d2_iso"] < 0.5))] ** 2)
reference_errors = np.sqrt(reference_sumw2)
reference_bin_centers = 0.5 * (reference_bin_edges[:-1] + reference_bin_edges[1:])  

reference_fit_hist = np.asarray(reference_fit_hist)
reference_bin_centers = np.asarray(reference_bin_centers)
reference_errors = np.asarray(reference_errors)

reference_mean = np.average(data["H_mass"][mask], weights=mc_weights[mask])
weighted_variance = np.average((data["H_mass"][mask] - reference_mean)**2, weights=mc_weights[mask])
reference_std = np.sqrt(weighted_variance)

print("Reference mean:", reference_mean)
print("Reference std:", reference_std)

gauss_popt, gauss_pcov = curve_fit(
    gaussian, 
    reference_bin_centers[reference_fit_hist > 0], 
    reference_fit_hist[reference_fit_hist > 0],
    p0=[np.max(reference_fit_hist[reference_fit_hist > 0]), reference_mean, reference_std],
    sigma=reference_errors[reference_fit_hist > 0], 
    absolute_sigma=True
)

crystalball_popt, crystalball_pcov = curve_fit(
    crystalball_func, 
    reference_bin_centers[reference_fit_hist > 0], 
    reference_fit_hist[reference_fit_hist > 0],
    p0=[np.max(reference_fit_hist[reference_fit_hist > 0]), 2, 1.5, reference_mean, reference_std],
    sigma=reference_errors[reference_fit_hist > 0], 
    absolute_sigma=True
)

popt, pcov = curve_fit(
    combined_model, 
    reference_bin_centers[reference_fit_hist > 0], 
    reference_fit_hist[reference_fit_hist > 0],
    p0=[np.max(reference_fit_hist[reference_fit_hist > 0]), 0.5, reference_mean, 2, 1.5, 2, 10],
    sigma=reference_errors[reference_fit_hist > 0], 
    absolute_sigma=True
)
'''
fig, axs = plt.subplots()
plt.hist(data["H_mass"][mask & ((data["d1_iso"] < 0.5) & (data["d2_iso"] < 0.5))], bins=40, range=(115, 135), weights=mc_weights[mask & ((data["d1_iso"] < 0.5) & (data["d2_iso"] < 0.5))], histtype='step', label='Signal', color='blue')

xfit = np.linspace(reference_bin_centers.min(), reference_bin_centers.max(), 200)
plt.plot(xfit, combined_model(xfit, *popt), label='Combined')
plt.plot(xfit, popt[1]*gaussian(xfit, popt[0], popt[2], popt[3]), label='Gaussian')
plt.plot(xfit, (1 - popt[1])*crystalball_func(xfit, popt[0], popt[4], popt[5], popt[2], popt[6]), label='Crystal Ball')
#plt.plot(xfit, gaussian(xfit, *gauss_popt), label='Gaussian', color='green')
#plt.plot(xfit, crystalball_func(xfit, *crystalball_popt), label='Crystal Ball', color='orange')

#plt.show()
'''
combined_chi_2 = np.sum(((reference_fit_hist[reference_fit_hist > 0] - combined_model(reference_bin_centers[reference_fit_hist > 0], *popt)) ** 2) / reference_errors[reference_fit_hist > 0] ** 2)
gauss_chi_2 = np.sum(((reference_fit_hist[reference_fit_hist > 0] - gaussian(reference_bin_centers[reference_fit_hist > 0], *gauss_popt)) ** 2) / reference_errors[reference_fit_hist > 0] ** 2)
crystalball_chi_2 = np.sum(((reference_fit_hist[reference_fit_hist > 0] - crystalball_func(reference_bin_centers[reference_fit_hist > 0], *crystalball_popt)) ** 2) / reference_errors[reference_fit_hist > 0] ** 2)

combined_dof = len(reference_fit_hist[reference_fit_hist > 0]) - len(popt)
gauss_dof = len(reference_fit_hist[reference_fit_hist > 0]) - len(gauss_popt)
crystalball_dof = len(reference_fit_hist[reference_fit_hist > 0]) - len(crystalball_popt)

combined_chi2_dof = combined_chi_2 / combined_dof
gauss_chi2_dof = gauss_chi_2 / gauss_dof
crystalball_chi2_dof = crystalball_chi_2 / crystalball_dof

print(combined_chi2_dof, gauss_chi2_dof, crystalball_chi2_dof)

iso_cut_val = 0
step = 0.01

iso_cuts = []
efficiencies = []
combined_chi2_dof = []
gauss_chi2_dof = []
crystalball_chi2_dof = []

eff_1_counter = 0

while True:
    iso_cut = (data["d1_iso"] < iso_cut_val) & (data["d2_iso"] < iso_cut_val)
    hist, _ = np.histogram(data["H_mass"][mask & iso_cut], bins=65, range=(70, 200), weights=mc_weights[mask & iso_cut])
    efficiency = np.sum(hist) / reference_integral

    iso_cuts.append(iso_cut_val)
    efficiencies.append(efficiency)

    if ((efficiency > 0) & (np.isclose(iso_cut_val % 0.1, 0))):
        fit_hist, bin_edges = np.histogram(data["H_mass"][mask & iso_cut], bins=40, range=(115, 135), weights=mc_weights[mask & iso_cut])
        sumw2, _ = np.histogram(data["H_mass"][mask & iso_cut], bins=40, range=(115, 135), weights=mc_weights[mask & iso_cut] ** 2)
        errors = np.sqrt(sumw2)
        bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

        fit_hist = np.asarray(fit_hist)
        bin_centers = np.asarray(bin_centers)
        errors = np.asarray(errors)

        gauss_popt, _ = curve_fit(
            gaussian, 
            bin_centers[fit_hist > 0], 
            fit_hist[fit_hist > 0], 
            p0=[np.max(fit_hist[fit_hist > 0]), 125, 2], 
            sigma=errors[fit_hist > 0], 
            absolute_sigma=True
        )

        crystalball_popt, _ = curve_fit(
            crystalball_func, 
            bin_centers[fit_hist > 0], 
            fit_hist[fit_hist > 0], 
            p0=[np.max(fit_hist[fit_hist > 0]), 2, 1.5, 125, 2], 
            sigma=errors[fit_hist > 0], 
            absolute_sigma=True
        )

        popt, pcov = curve_fit(
            combined_model, 
            bin_centers[fit_hist > 0], 
            fit_hist[fit_hist > 0], 
            p0=[np.max(fit_hist[fit_hist > 0]), 0.75, gauss_popt[1], gauss_popt[2], crystalball_popt[1], crystalball_popt[2], crystalball_popt[4]],
            sigma=errors[fit_hist > 0], 
            absolute_sigma=True,
            maxfev=10000
        )
        '''
        if np.isclose(iso_cut_val, 0.5):
            print(popt[1])
            fig, axs = plt.subplots()
            plt.hist(data["H_mass"][mask & iso_cut], bins=40, range=(115, 135), weights=mc_weights[mask & iso_cut], histtype='step', label='Signal', color='blue')

            xfit = np.linspace(bin_centers.min(), bin_centers.max(), 200)
            plt.plot(xfit, combined_model(xfit, *popt), label='Combined')
            plt.plot(xfit, popt[1]*gaussian(xfit, popt[0], popt[2], popt[3]), label='Gaussian')
            plt.plot(xfit, (1 - popt[1])*crystalball_func(xfit, popt[0], popt[4], popt[5], popt[2], popt[6]), label='Crystal Ball')
            #plt.plot(xfit, gaussian(xfit, *gauss_popt), label='Gaussian')
            #plt.plot(xfit, crystalball_func(xfit, *crystalball_popt), label='Crystal Ball')
            plt.legend()
            plt.show()
        '''
        chi2 = np.sum(((fit_hist[fit_hist > 0] - combined_model(bin_centers[fit_hist > 0], *popt)) ** 2) / errors[fit_hist > 0]**2)
        gauss_chi2 = np.sum(((fit_hist[fit_hist > 0] - gaussian(bin_centers[fit_hist > 0], *gauss_popt)) ** 2) / errors[fit_hist > 0]**2)
        crystalball_chi2 = np.sum(((fit_hist[fit_hist > 0] - crystalball_func(bin_centers[fit_hist > 0], *crystalball_popt)) ** 2) / errors[fit_hist > 0]**2)

        dof = len(fit_hist[fit_hist > 0]) - len(popt)
        gauss_dof = len(fit_hist[fit_hist > 0]) - len(gauss_popt)
        crystalball_dof = len(fit_hist[fit_hist > 0]) - len(crystalball_popt)

        combined_chi2_dof.append(chi2 / dof)
        gauss_chi2_dof.append(gauss_chi2 / gauss_dof)
        crystalball_chi2_dof.append(crystalball_chi2 / crystalball_dof)

    iso_cut_val += step

    if efficiency >= 0.995:
        eff_1_counter += 1

    if eff_1_counter >= 10:
        break

iso_cuts = np.array(iso_cuts)
efficiencies = np.array(efficiencies)
combined_chi2_dof = np.array(combined_chi2_dof)
gauss_chi2_dof = np.array(gauss_chi2_dof)
crystalball_chi2_dof = np.array(crystalball_chi2_dof)

fig, axs = plt.subplots()

#smoothed_efficiencies = moving_average(efficiencies, window_size=10)

#plt.plot(iso_cuts, smoothed_efficiencies, linestyle='-', color='red')
plt.plot(iso_cuts, efficiencies, linestyle='-', color='blue')

eff_at_0_5 = efficiencies[np.isclose(iso_cuts, 0.5)][0]
print("Efficiency at 0.5:", eff_at_0_5)
x_vals = np.linspace(0, 0.5, 1000)
y_vals = np.linspace(0, eff_at_0_5, 1000)

plt.plot(0.5 * np.ones_like(y_vals), y_vals, linestyle='--', color='black')
plt.plot(x_vals, eff_at_0_5 * np.ones_like(x_vals), linestyle='--', color='black')

axs.set_xlabel("Isolation Cut", fontsize=20)
axs.set_ylabel("Signal Efficiency", fontsize=20)

axs.set_xlim(0, np.max(iso_cuts))
axs.set_ylim(0, 1.05)

axs.set_title(r"$\it{" + "Private \: work" + "}$" +" "+ r"$\bf{" + "(CMS \: simulation)" + "}$",loc="left", fontsize=20)
axs.set_title(f"59.8 fb$^{{-1}}$, 2018 (13 $\mathrm{{{{TeV}}}}$)",loc="right", fontsize=20)

plt.savefig(f"/web/jhornung/public_html/analysis_plots/2018/signal_efficiency_{scope}.png")
plt.savefig(f"/web/jhornung/public_html/analysis_plots/2018/signal_efficiency_{scope}.pdf") 

fig, axs = plt.subplots()
plt.plot(iso_cuts[(efficiencies > 0) & np.isclose(iso_cuts % 0.1, 0)], gauss_chi2_dof, linestyle='-', color='green', label=r"$\chi^2 / dof$ (Gaussian)")
plt.plot(iso_cuts[(efficiencies > 0) & np.isclose(iso_cuts % 0.1, 0)], crystalball_chi2_dof, linestyle='-', color='orange', label=r"$\chi^2 / dof$ (Crystal Ball)")
plt.plot(iso_cuts[(efficiencies > 0) & np.isclose(iso_cuts % 0.1, 0)], combined_chi2_dof, linestyle='-', color='blue', label=r"$\chi^2 / dof$ (Combined)")
axs.set_xlabel("Isolation Cut", fontsize=20)
axs.set_ylabel(r"$\chi^2 / dof$", fontsize=20)

axs.set_xlim(left=np.min(iso_cuts[(efficiencies > 0) & (np.isclose(iso_cuts % 0.1, 0))]), right=np.max(iso_cuts[(efficiencies > 0) & (np.isclose(iso_cuts % 0.1, 0))]))
axs.set_ylim(top=15, bottom=0)

plt.legend(loc='upper right', fontsize=20)

axs.set_title(r"$\it{" + "Private \: work" + "}$" +" "+ r"$\bf{" + "(CMS \: simulation)" + "}$",loc="left", fontsize=20)
axs.set_title(f"59.8 fb$^{{-1}}$, 2018 (13 $\mathrm{{{{TeV}}}}$)",loc="right", fontsize=20)

