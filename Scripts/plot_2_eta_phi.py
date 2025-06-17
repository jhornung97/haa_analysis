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

signal = ["HaaKKKK_M125_13TeV_powheg_pythia8-RunIISummer20UL18MiniAODv2-106X"]
signal_labels = [r'sig']
mc_signal = hl.data(mc_loc, signal, signal_labels, scope)
signal_branches = mc_signal.branches_dict

print(signal_branches)

config = {
        "scope": scope,
        "mode": "reco to truth 2D",
        "var": "d1_eta",
        "var2": "d1_phi",
        "xlower": -2.5,
        "xupper": 2.5,
        "ylower": -3.2,
        "yupper": 3.2,
        "bins": 0,
        "diff": False,
        "heatmap": False,
        "signal": signal_branches
        }

plot = hl.hist(config)

plot.plot_2d()
plt.show()
