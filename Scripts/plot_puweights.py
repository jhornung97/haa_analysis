import uproot
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.patches as mpatches
import mplhep as hep
import utility as util
import sys
import pickle, lz4.frame
import HistLib as hl

plt.style.use(hep.style.CMS)
plt.rcParams.update({'font.size': 24})

mc_loc = "/ceph/jhornung/MC_2018/2018"
scope = "mm"
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

signal = ["HaaKKKK_M125_13TeV_powheg_pythia8-RunIISummer20UL18MiniAODv2-106X"]
signal_labels = [r'$H\rightarrow aa \rightarrow 4\,K$, $m_a = 1.5\, GeV$']
mc_signal = hl.data(mc_loc, signal, signal_labels, scope)
signal_branches = mc_signal.branches_dict

config = {
        "scope": scope,
        "mode": "mc",
        "mc": bkg_branches,
        "signal": signal_branches,
        "z_pt_weights": z_pt_weights,
        "var": "puweight",
        "diff": False,
        "heatmap": False,
        "lower": 0,
        "upper": 2,
        "bins": 20,
        }

hist = hl.hist(config)
hist.plot_mc(save_as=f"/web/jhornung/public_html/ps_mass_cut/puweights_{scope}.png")
plt.show()