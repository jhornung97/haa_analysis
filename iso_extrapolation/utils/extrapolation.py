import uproot
from concurrent.futures import ThreadPoolExecutor
import awkward as ak
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import chi2
import mplhep as hep
import sys
import csv
import os
import json

sys.path.append('/work/jhornung/Haa/Scripts/')
import HistLib as hl
import utility as util

plt.style.use(hep.style.CMS)

def read_file(filename, data_type="data"):
    if data_type == "data":
        with uproot.open(filename) as file:
            tree = file["ntuple"]
            data = tree.arrays([
                "H_mass",
                "pt_vis",
                "m_vis",
                "H_eta",
                "d1_pt",
                "d2_pt",
                "ps_1_mass",
                "ps_2_mass",
                "d1_iso",
                "d2_iso",
                "trg_single_mu" if scope == 'mm' else "trg_single_el",
            ])
            return data
    elif data_type == "mc":
        with uproot.open(filename) as file:
            tree = file["ntuple"]
            data = tree.arrays([
                "H_mass",
                "pt_vis",
                "m_vis",
                "H_eta",
                "d1_pt",
                "d2_pt",
                "ps_1_mass",
                "ps_2_mass",
                "d1_iso",
                "d2_iso",
                "id_wgt_mu_1" if scope == 'mm' else "id_wgt_ele_1",
                "id_wgt_mu_2" if scope == 'mm' else "id_wgt_ele_2",
                "iso_wgt_mu_1" if scope == 'mm' else None,
                "iso_wgt_mu_2" if scope == 'mm' else None,
                "trigger_wgt_mu_1" if scope == 'mm' else "trigger_wgt_ele_1",
                "trigger_wgt_mu_2" if scope == 'mm' else "trigger_wgt_ele_2",
                "evtweight",
                "puweight",
                "genWeight",
                "trg_single_mu" if scope == 'mm' else "trg_single_el",
            ])
            return data
    
def model(x, a, x0, b):
    return a * (1 - np.exp(-x/x0)) ** b

def model_errors(x, a, x0, b, cov, withCorr=True):
    da = (1 - np.exp(-x/x0)) ** b
    dx0 = - a * b * x * (1 - np.exp(-x/x0)) ** (b - 1) * np.exp(-x/x0) / x0**2
    db = a * (1 - np.exp(-x/x0)) ** (b) * np.log(1 - np.exp(-x/x0) + 1e-12)

    diag_term = da**2 * cov[0,0] + dx0**2 * cov[1,1] + db**2 * cov[2,2]
    non_diag_terms = 2 * (da * cov[0,1] * dx0 + da * cov[0,2] * db + dx0 * cov[1,2] * db)
    if withCorr:
        return np.sqrt(diag_term + non_diag_terms)
    else:
        return np.sqrt(diag_term) 
    
year = sys.argv[1]
scope = sys.argv[2]

if scope == 'mm' and int(year) <= int("2018"):
    path = f"/ceph/jhornung/Data_2018/{year}/single_muon_data/{scope}/"
    exclude = "single_muon_data.root"
elif scope == 'ee' and year == '2018':
    path = f"/ceph/jhornung/Data_2018/{year}/egamma_data/{scope}/"
    exclude = "egamma_data.root"
elif scope == 'ee' and (year == '2017' or year == '2016preVFP' or year == '2016postVFP'):
    path = f"/ceph/jhornung/Data_2018/{year}/single_electron_data/{scope}/"
    exclude = "single_electron_data.root"
elif scope == 'mm' and year == '2024':
    path = f"/ceph/jhornung/Data_2018/{year}/muon_data/{scope}/"
    exclude = "muon_data.root"
elif scope == 'ee' and year == '2024':
    path = f"/ceph/jhornung/Data_2018/{year}/egamma_data/{scope}/"
    exclude = "egamma_data.root"

filelist = [path + f for f in os.listdir(path) if f.endswith(".root") and f != exclude]

with ThreadPoolExecutor(max_workers=2) as executor:
    results = list(executor.map(read_file, filelist))

data = ak.concatenate(results)

H_mass_mask = ((data["H_mass"] > 70) & (data["H_mass"] < 110)) | ((data["H_mass"] > 140) & (data["H_mass"] < 200))
H_eta_mask = (np.abs(data["H_eta"]) < 2.4)
leading_k_mask = (data["d1_pt"] > 10) & (data["d2_pt"] > 10)
ps_mass_mask = (data["ps_1_mass"] < 3) & (data["ps_2_mass"] < 3) & (np.abs(data["ps_1_mass"] - data["ps_2_mass"]) < 0.06)
pt_vis_mask = (data["pt_vis"] > 30)
m_vis_mask = (data["m_vis"] > 75) & (data["m_vis"] < 105)
trigger_branch = "trg_single_mu" if scope == 'mm' else "trg_single_el"
trigger_mask = data[trigger_branch] == 1
mask = H_mass_mask & H_eta_mask & leading_k_mask & ps_mass_mask & pt_vis_mask & m_vis_mask & trigger_mask

# fixed fit range in iso for all eras; the signal region iso < x_extrap is never part of the fit
x_extrap = 0.5
iso_max = 10
iso_min = 1
step = 0.5

ref_hist, edges = np.histogram(data["H_mass"][mask], bins=65, range=(70, 200))
ref_bin_count = np.sum(ref_hist)

bin_centers = (edges[:-1] + edges[1:]) / 2

counts = []
iso_cuts = []
for cut in np.arange(iso_max, iso_min - step / 2, -step):
    iso_mask = (data["d1_iso"] < cut) & (data["d2_iso"] < cut)
    hist, _ = np.histogram(data["H_mass"][mask & iso_mask], bins=65, range=(70, 200))
    iso_cuts.append(cut)
    counts.append(np.sum(hist))

iso_cuts = np.array(iso_cuts)
counts = np.array(counts)
error = np.sqrt(counts)
# the points are cumulative (iso < x contains iso < x - step), so for Poisson counts cov_ij = min(N_i, N_j)
cov = np.minimum.outer(counts, counts).astype(float)

fit_bounds = ([0, 1e-3, 1e-3], [np.inf, np.inf, np.inf])

popt, pcov = curve_fit(model, iso_cuts, counts, sigma=cov, p0=[ref_bin_count, 1, 1], absolute_sigma=True, bounds=fit_bounds)

x_fit = np.linspace(0, 10.5, 2001)
y_fit = model(x_fit, *popt)
error_fit = model_errors(x_fit, *popt, pcov, withCorr=True)

extrapolated_value = model(x_extrap, *popt)
extrapolated_error = model_errors(x_extrap, *popt, pcov, withCorr=True)

def fit_extrapolation(f, x, y, p0, bounds):
    """Correlated fit of f to the cumulative counts y(x); returns parameters, fit quality and the value at x_extrap."""
    cov_xy = np.minimum.outer(y, y).astype(float)
    try:
        p, _ = curve_fit(f, x, y, sigma=cov_xy, p0=p0, absolute_sigma=True, bounds=bounds, maxfev=20000)
    except RuntimeError:
        return None
    res = y - f(x, *p)
    chi2_val = res @ np.linalg.solve(cov_xy, res)
    ndf = len(y) - len(p)
    return {"func": f, "popt": p, "chi2": chi2_val, "ndf": ndf, "pvalue": chi2.sf(chi2_val, ndf), "extrapolation": f(x_extrap, *p)}

def print_fit(name, fit):
    if fit is None:
        print(f"{name:34s}: fit failed")
        return
    print(f"{name:34s}: N({x_extrap}) = {fit['extrapolation']:.1f}  chi2/ndf = {fit['chi2']:.1f}/{fit['ndf']}  p = {fit['pvalue']:.3f}  params = {np.array2string(fit['popt'], precision=3)}")

# 1) shape systematic: nested 4-parameter extension of the nominal model (reduces to it for c = 0)
def model_ext(x, a, x0, b, c):
    return model(x, a, x0, b) * (1 + c * x)

# c > -1/max(x) keeps the correction factor positive over the whole range
ext_bounds = (fit_bounds[0] + [-1 / iso_cuts.max()], fit_bounds[1] + [np.inf])
nominal_fit = fit_extrapolation(model, iso_cuts, counts, list(popt), fit_bounds)
ext_fit = fit_extrapolation(model_ext, iso_cuts, counts, list(popt) + [0], ext_bounds)
print("nominal and nested extension:")
print_fit("nominal", nominal_fit)
print_fit("nominal x (1 + c x)", ext_fit)
shape_syst = abs(ext_fit["extrapolation"] - extrapolated_value) if ext_fit is not None else 0.0
if ext_fit is not None:
    print(f"  delta chi2 (1 dof) = {nominal_fit['chi2'] - ext_fit['chi2']:.2f}")

# 2) fit range systematic: refit the nominal model with the highest / lowest iso points removed
# iso_cuts is ordered from high to low, so [n:] drops the n highest and [:-n] the n lowest points
range_variations = {
    f"iso $\\leq$ {iso_cuts[2]:g}": slice(2, None),
    f"iso $\\leq$ {iso_cuts[4]:g}": slice(4, None),
    f"iso $\\geq$ {iso_cuts[-2]:g}": slice(None, -1),
    f"iso $\\geq$ {iso_cuts[-3]:g}": slice(None, -2),
}
range_fits = {}
print("fit range variations (nominal model):")
for name, sl in range_variations.items():
    if len(iso_cuts[sl]) <= len(popt) + 1:
        continue
    range_fits[name] = fit_extrapolation(model, iso_cuts[sl], counts[sl], list(popt), fit_bounds)
    print_fit(name, range_fits[name])
range_fits = {name: fit for name, fit in range_fits.items() if fit is not None}
range_syst = max([abs(fit["extrapolation"] - extrapolated_value) for fit in range_fits.values()], default=0.0)

# 3) cross-check only: alternative 3-parameter forms (not used in the systematic, the
# exp/tanh forms can degenerate into a pure power law a (x/x0)^b for x0 -> inf)
alt_models = {
    r"$a\,(1-e^{-(x/x_0)^b})$": lambda x, a, x0, b: a * (1 - np.exp(-(x / x0)**b)),
    r"$a\,x^b/(x_0^b+x^b)$": lambda x, a, x0, b: a * x**b / (x0**b + x**b),
    r"$a\,\tanh(x/x_0)^b$": lambda x, a, x0, b: a * np.tanh(x / x0)**b,
}
print("alternative functions (cross-check):")
alt_fits = {}
for name, f in alt_models.items():
    alt_fits[name] = fit_extrapolation(f, iso_cuts, counts, [ref_bin_count, 1, 1], fit_bounds)
    print_fit(name, alt_fits[name])
alt_fits = {name: fit for name, fit in alt_fits.items() if fit is not None}

func_syst = np.sqrt(shape_syst**2 + range_syst**2)
total_error = np.sqrt(extrapolated_error**2 + func_syst**2)
print(f"N({x_extrap}) = {extrapolated_value:.1f} +- {extrapolated_error:.1f} (stat) +- {shape_syst:.1f} (shape) +- {range_syst:.1f} (fit range)")

# curves shown in the plot for the systematic variations
syst_fits = dict(range_fits)
if ext_fit is not None:
    syst_fits[r"nominal $\times\,(1+c\,x)$"] = ext_fit

def fit_summary(fit):
    return {"extrapolation": fit["extrapolation"], "chi2": fit["chi2"], "ndf": int(fit["ndf"]), "pvalue": fit["pvalue"], "params": fit["popt"].tolist()}

result = {
    "normalization": extrapolated_value,
    "normalization_error": extrapolated_error,
    "normalization_shape_syst": shape_syst,
    "normalization_range_syst": range_syst,
    "normalization_func_syst": func_syst,
    "normalization_total_error": total_error,
    "nested_extension": fit_summary(ext_fit) if ext_fit is not None else None,
    "fit_range_variations": {name: fit_summary(fit) for name, fit in range_fits.items()},
    "alternative_functions": {name: fit_summary(fit) for name, fit in alt_fits.items()},
}

output_file = f"/work/jhornung/Haa/iso_extrapolation/bkg_normalization.json"

if os.path.exists(output_file):
    with open(output_file, 'r') as f:
        data = json.load(f)
else:
    data = {}

if year not in data:
    data[year] = {}
if scope not in data[year]:
    data[year][scope] = {}

data[year][scope] = result

with open(output_file, 'w') as f:
    json.dump(data, f, indent=2)

ratio = (counts - model(iso_cuts, *popt)) / model_errors(iso_cuts, *popt, pcov, withCorr=True)
ratio_errors = np.sqrt(error**2/model_errors(iso_cuts, *popt, pcov, withCorr=True)**2)

mc_loc = f"/ceph/jhornung/MC_2018/{year}"

if year == "2016preVFP":
        bkgs = [
                "DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8+RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1+MINIAODSIM",
                "ttbar",
                "diboson",
                "wh",
                "single_top"
                #"qcd"
                ]
elif year == "2016postVFP":
        bkgs = [
                "DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8+RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1+MINIAODSIM",
                "ttbar",
                "diboson",
                "wh",
                "single_top"
                #"qcd"
                ]
elif year == "2017":
        bkgs = [
                "DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8+RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v2+MINIAODSIM",
                "ttbar",
                "diboson",
                "wh",
                "single_top"
                #"qcd"
                ]
elif year == "2018":
        bkgs = [
                "DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8+RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2+MINIAODSIM",
                "ttbar",
                "diboson",
                "wh",
                "single_top"
                #"qcd"
                ]
elif year == "2024":
        bkgs = [
                #"DYto2E-2Jets_Bin-MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8-RunIII2024Summer24MiniAODv6-150X",
                #"DYto2Mu-2Jets_Bin-MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8-RunIII2024Summer24MiniAODv6-150X",
                #"DYto2Tau-2Jets_Bin-MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8-RunIII2024Summer24MiniAODv6-150X",
                "dy",
                "TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8-RunIII2024Summer24MiniAODv6-150X",
                "diboson",
                "wh",
                #"single_top"
                #"qcd"
                ]

if int(year[:4]) <= 2018:
        bkg_labels = [
                'Drell-Yan', 
                r'$t\,\bar{t}$', 
                r'Diboson', 
                r'$W\,H\rightarrow \ell\,\nu\, b\,\bar{b}$', 
                r'Single Top'
                #r'QCD'
                ]
else:
        bkg_labels = [
                #'DYto2E-2Jets',
                #'DYto2Mu-2Jets',
                #'DYto2Tau-2Jets',
                'Drell-Yan', 
                r'$t\,\bar{t}$', 
                'Diboson', 
                r'$W\,H\rightarrow \ell\,\nu\, b\,\bar{b}$'
                #r'Single Top',
                #r'QCD'
                ]

cache_dir = f"/ceph/jhornung/cache/{year}/{scope}"

mc_bkg = hl.data(mc_loc, bkgs, bkg_labels, "zh", scope, data_type="mc", cache_dir=cache_dir, use_cache=True, write_cache=False, max_workers=4)
bkg_branches = mc_bkg.branches_dict
z_pt_weights = mc_bkg.get_Z_pt_weights()

mc_masks = {}
mc_weights = {}
mc_counts = {}
mc_errors = {}

lumi = None
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

for key, branch in bkg_branches.items():
    H_mass_mask = ((branch["H_mass"] > 70) & (branch["H_mass"] < 110)) | ((branch["H_mass"] > 140) & (branch["H_mass"] < 200))
    H_eta_mask = (np.abs(branch["H_eta"]) < 2.4)
    leading_k_mask = (branch["d1_pt"] > 10) & (branch["d2_pt"] > 10)
    ps_mass_mask = (branch["ps_1_mass"] < 3) & (branch["ps_2_mass"] < 3) & (np.abs(branch["ps_1_mass"] - branch["ps_2_mass"]) < 0.06)
    pt_vis_mask = (branch["pt_vis"] > 30)
    m_vis_mask = (branch["m_vis"] > 75) & (branch["m_vis"] < 105)
    trigger_mask = branch[trigger_branch] == 1
    mask = H_mass_mask & H_eta_mask & leading_k_mask & ps_mass_mask & pt_vis_mask & m_vis_mask & trigger_mask
    mc_masks[key] = mask

    if scope == 'mm':
        trigger_sf = branch["trigger_wgt_mu_1"]
        sfs = branch["id_wgt_mu_1"] * branch["id_wgt_mu_2"] * branch["iso_wgt_mu_1"] * branch["iso_wgt_mu_2"] * trigger_sf
    elif scope == 'ee':
        trigger_sf = branch["trigger_wgt_ele_1"]
        sfs = branch["id_wgt_ele_1"] * branch["id_wgt_ele_2"] #* trigger_sf
    signs = np.array(list(map(util.weights, branch["genWeight"])))
    if key == 'Drell-Yan':
        mc_weights[key] = lumi * 1e3 * sfs * signs * branch["evtweight"] * branch["puweight"] #* z_pt_weights
    else:
        mc_weights[key] = lumi * 1e3 * sfs * signs * branch["evtweight"] * branch["puweight"]

    mc_count = []
    mc_error = []

    for i, iso_cut in enumerate(iso_cuts):
        iso_mask = (branch["d1_iso"] < iso_cut) & (branch["d2_iso"] < iso_cut)
        hist, _ = np.histogram(branch["H_mass"][mask & iso_mask], weights=mc_weights[key][mask & iso_mask], bins=65, range=(70, 200))
        stat_errs, _ = np.histogram(branch["H_mass"][mask & iso_mask], weights=mc_weights[key][mask & iso_mask]**2, bins=65, range=(70, 200)) 
        mc_count.append(np.sum(hist))
        mc_error.append(np.sqrt(np.sum(stat_errs)))

    mc_counts[key] = np.array(mc_count)
    mc_errors[key] = np.array(mc_error)

# normalise total MC to data at each iso cut (one common factor for all backgrounds)
# c_k = m_k * D / M with M = sum_j m_j; propagate MC stat of all m_j and Poisson error of D
raw_counts = {key: mc_count.copy() for key, mc_count in mc_counts.items()}
raw_var = {key: mc_error**2 for key, mc_error in mc_errors.items()}
M = np.sum(list(raw_counts.values()), axis=0)
var_M = np.sum(list(raw_var.values()), axis=0)
norm = counts / M
for key in mc_counts:
    f = raw_counts[key] / M
    mc_counts[key] = raw_counts[key] * norm
    mc_errors[key] = np.sqrt(norm**2 * ((1 - f)**2 * raw_var[key] + f**2 * (var_M - raw_var[key])) + f**2 * counts)

fig, axs = plt.subplots(2, 1, sharex=True, gridspec_kw={"height_ratios": [3, 1]})

iso_cuts_fill = np.insert(iso_cuts, 0, 10.5)
iso_cuts_fill = np.append(iso_cuts_fill, np.min(iso_cuts) - 0.5)

previous_counts = np.zeros_like(iso_cuts_fill) 

colors = ['seagreen', 'gold', 'darkorange', 'orange', 'mediumspringgreen']

for i, (key, mc_count) in enumerate(mc_counts.items()):
    mc_count = np.insert(mc_count, 0, 0)
    mc_count = np.append(mc_count, 0)
    axs[0].fill_between(iso_cuts_fill[::-1], previous_counts[::-1], previous_counts[::-1] + mc_count[::-1], step='mid', label=key, alpha=1, color=colors[i], linewidth=0)
    previous_counts += mc_count

total_mc_count = np.sum(list(mc_counts.values()), axis=0)
# components are correlated through the normalisation: sum_k c_k = D, so the MC stat terms cancel
total_mc_error = np.sqrt(counts)

iso_cuts_edges = np.insert(iso_cuts, 0, iso_cuts[0] + step / 2)
iso_cuts_edges = np.append(iso_cuts_edges, iso_cuts[-1] - step / 2)

total_mc_count_fill = np.insert(total_mc_count, 0, total_mc_count[0])
total_mc_count_fill = np.append(total_mc_count_fill, total_mc_count[-1])
total_mc_error_fill = np.insert(total_mc_error, 0, total_mc_error[0])
total_mc_error_fill = np.append(total_mc_error_fill, total_mc_error[-1])

axs[0].fill_between(iso_cuts_edges, total_mc_count_fill - total_mc_error_fill, total_mc_count_fill + total_mc_error_fill, step='mid', facecolor='none', edgecolor='gray', hatch='///', linewidth=0, label='MC stat.')

axs[0].errorbar(iso_cuts, counts, yerr=error, fmt='o', color='black')
axs[0].plot(x_fit, y_fit, color='red', label='Fit')
axs[0].fill_between(x_fit, y_fit - error_fit, y_fit + error_fit, facecolor='red', alpha=0.25, label='Fit Error', linewidth=0)

alt_colors = ['royalblue', 'darkviolet', 'teal', 'saddlebrown', 'deeppink']
for color, (name, fit) in zip(alt_colors, syst_fits.items()):
    axs[0].plot(x_fit, fit["func"](x_fit, *fit["popt"]), color=color, linestyle='--', linewidth=1.5, label=name)

# outer bar: stat (+) function choice, inner bar: stat only
axs[0].errorbar(x_extrap, extrapolated_value, yerr=total_error, fmt='none', color='red', capsize=6, elinewidth=1.5)
axs[0].errorbar(x_extrap, extrapolated_value, yerr=extrapolated_error, fmt='o', color='red', capsize=3, elinewidth=3, label='Extrapolation')

# zoom on the extrapolation region, where the function choice matters
inset = axs[0].inset_axes([0.62, 0.08, 0.35, 0.4])
zoom = x_fit <= 2.5
inset.fill_between(x_fit[zoom], (y_fit - error_fit)[zoom], (y_fit + error_fit)[zoom], facecolor='red', alpha=0.25, linewidth=0)
inset.plot(x_fit[zoom], y_fit[zoom], color='red')
for color, fit in zip(alt_colors, syst_fits.values()):
    inset.plot(x_fit[zoom], fit["func"](x_fit[zoom], *fit["popt"]), color=color, linestyle='--', linewidth=1.5)
inset.errorbar(x_extrap, extrapolated_value, yerr=total_error, fmt='none', color='red', capsize=6, elinewidth=1.5)
inset.errorbar(x_extrap, extrapolated_value, yerr=extrapolated_error, fmt='o', color='red', capsize=3, elinewidth=3)
in_zoom = iso_cuts <= 2.5
inset.errorbar(iso_cuts[in_zoom], counts[in_zoom], yerr=error[in_zoom], fmt='o', color='black', markersize=4)
inset.set_xlim(0, 2.5)
inset.tick_params(labelsize=12)
axs[0].set_ylabel("Event yield", fontsize=20)
axs[0].set_xlim(0, 10.5)
axs[0].set_ylim(bottom=0, top=1.25 * max(counts))
# three legend columns: samples/fit | fit range variations | nested extension + extrapolation
handles, labels = axs[0].get_legend_handles_labels()
legend_columns = [[], [], []]
for handle, label in zip(handles, labels):
    if label in range_fits:
        legend_columns[1].append((handle, label))
    elif label in syst_fits or label.startswith('Extrapolation'):
        legend_columns[2].append((handle, label))
    else:
        legend_columns[0].append((handle, label))
# legend fills column by column, so pad every column to the same length with invisible entries
n_rows = max(len(column) for column in legend_columns)
blank = (plt.Line2D([], [], alpha=0), "")
legend_entries = [entry for column in legend_columns for entry in column + [blank] * (n_rows - len(column))]
axs[0].legend(*zip(*legend_entries), loc='upper left', fontsize=17, ncol=3)
axs[0].set_title(r"$\it{" + "Private \: work" + "}$" +" "+ r"$\bf{" + "(CMS \: data/simulation)" + "}$",loc="left", fontsize=20)
if "2016" in year:
    axs[0].set_title(f"{lumi} fb$^{{-1}}$, 2016 (13 $\mathrm{{{{TeV}}}}$)",loc="right", fontsize=20)
elif int(year[:4]) <= 2018:
    axs[0].set_title(f"{lumi} fb$^{{-1}}$, {year} (13 $\mathrm{{{{TeV}}}}$)",loc="right", fontsize=20)
else:
    axs[0].set_title(f"{lumi:.1f} fb$^{{-1}}$, {year} (13.6 $\mathrm{{{{TeV}}}}$)",loc="right", fontsize=20)

axs[1].errorbar(iso_cuts, ratio, yerr=ratio_errors, fmt='o', color='black')
axs[1].axhline(0, color='red', linestyle='--')
axs[1].fill_between(x_fit, -1, 1, color='red', alpha=0.25, linewidth=0)
for color, fit in zip(alt_colors, syst_fits.values()):
    axs[1].plot(x_fit[1:], (fit["func"](x_fit[1:], *fit["popt"]) - y_fit[1:]) / error_fit[1:], color=color, linestyle='--', linewidth=1.5)
axs[1].axvline(x_extrap, color='gray', linestyle=':', linewidth=1)
axs[1].set_ylabel(r"$\frac{Data - Fit}{Fit\ Error}$", fontsize=20)
axs[1].set_xlabel("Isolation Cut", fontsize=20)
axs[1].set_xlim(0, 10.25)
axs[1].set_ylim(-3, 3)

plt.tight_layout()
fig.subplots_adjust(hspace=0)

plt.savefig(f"/web/jhornung/public_html/analysis_plots/{year}/{scope}_iso_cut_extrapolation_test.png")
plt.savefig(f"/web/jhornung/public_html/analysis_plots/{year}/{scope}_iso_cut_extrapolation_test.pdf")
