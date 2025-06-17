import sys
import uproot
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.patches as mpatches
import mplhep as hep
import utility as util
import sys
#import pickle, lz4.frame
import HistLib as hl

plt.style.use(hep.style.CMS)
plt.rcParams.update({'font.size': 24})

scope = sys.argv[1]
year = sys.argv[2]
selection = sys.argv[3]

mc_loc = f"/ceph/jhornung/MC_2018/{year}"
signal_loc = "/ceph/jhornung/MC_2018/2018"
data_loc = f"/ceph/jhornung/Data_2018/{year}"

if year == "2016preVFP":
        bkgs = [
                "DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8+RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1+MINIAODSIM",
                "ttbar",
                "diboson",
                "wh",
                "single_top"
                #"qcd"
                ]
if year == "2016postVFP":
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


bkg_labels = [
        'Drell-Yan', 
        r'$t\,\bar{t}$', 
        r'Diboson', 
        r'$W\,H\rightarrow \ell\,\nu\, b\,\bar{b}$', 
        r'Single Top'
        #r'QCD'
        ]
mc_bkg = hl.data(mc_loc, bkgs, bkg_labels, scope)
bkg_branches = mc_bkg.branches_dict
z_pt_weights = mc_bkg.get_Z_pt_weights()
#mc_min_combs = mc_bkg.get_min_ps_mass_diff()

signal = ["HaaKKKK_M125_13TeV_powheg_pythia8-RunIISummer20UL18MiniAODv2-106X"]
signal_labels = [r'$Z(\ell\ell)H(aa\rightarrow 4K)/10$']
mc_signal = hl.data(signal_loc, signal, signal_labels, scope)
signal_branches = mc_signal.branches_dict
#signal_min_combs = mc_signal.get_min_ps_mass_diff()

if scope == "mm" or scope == "em":
    datasamples = ["single_muon_data"]
elif scope == "ee" and year == "2018":
    datasamples = ["egamma_data"]
elif scope == "ee" and ((year == "2017") or (year == "2016preVFP") or (year == "2016postVFP")):
    datasamples = ["single_electron_data"]

data_labels = [r"Data"]
data = hl.data(data_loc, datasamples, data_labels, scope)
data_branches = data.branches_dict
#data_min_combs = data.get_min_ps_mass_diff()

vars = ["H_mass", "H_eta", "m_vis", "pt_vis", "d1_pt", "d2_pt", "d1_iso", "d2_iso", "d3_iso", "d4_iso", "ps_1_mass", "ps_2_mass"]
lbs = [60, -3, 70, 20, 5, 5, 0, 0, 0, 0, 0, 0]
ubs = [220, 3, 110, 300, 100, 100, 1000, 1000, 1000, 1000, 5, 5]
nbs = [80, 60, 40, 56, 38, 38, 1000, 1000, 1000, 1000, 50, 50]

#var = ['PV_npvsGood']
#lbs = [0]
#ubs = [60]
#nbs = [60]

for var, lb, ub, nb in zip(vars, lbs, ubs, nbs):
        config = {
                "scope": scope,
                "mode": "data to mc",
                "era": f"{year}",
                "selection": selection,
                "mc": bkg_branches,
                "data": data_branches,
                "signal": signal_branches,
                "z_pt_weights": z_pt_weights,
                "var": var,
                "heatmap": False,
                "lower": lb,
                "upper": ub,    
                "bins": nb
            }
        hist = hl.hist(config)
        hist.plot_data_to_mc(save_as=f"/web/jhornung/public_html/analysis_plots/{config['era']}/{var}_{scope}_AN_{config['selection']}_{nb}_bins_tight_Delta_mKK_cut")
'''

config = {
        "scope": scope,
        "mode": "diff",
        "era": f"{year}",
        "selection": selection,
        "mc": bkg_branches,
        "data": data_branches,
        "signal": signal_branches,
        "z_pt_weights": z_pt_weights,
        "var_1": "ps_1_mass",
        "var_2": "ps_2_mass",
        "heatmap": False,
        "lower": -0.15,
        "upper": 0.15,    
        "bins": 60
    }

hist = hl.hist(config)
hist.plot_diff(save_as=f"/web/jhornung/public_html/analysis_plots/{config['era']}/ps_mass_diff_{scope}_AN_{config['selection']}")

config = {
        "scope": scope,
        "mode": "append",
        "era": f"{year}",
        "selection": selection,
        "mc": bkg_branches,
        "data": data_branches,
        "signal": signal_branches,
        "z_pt_weights": z_pt_weights,
        "var_1": "ps_1_mass",
        "var_2": "ps_2_mass",
        "heatmap": False,
        "lower": 0,
        "upper": 5,    
        "bins": 50
    }

hist = hl.hist(config)
hist.plot_both()#save_as=f"/web/jhornung/public_html/analysis_plots/{config['era']}/ps_mass_append_{scope}_AN_{config['selection']}")
'''
plt.show()
