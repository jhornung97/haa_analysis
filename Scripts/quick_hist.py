import uproot
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.patches as mpatches
from matplotlib.legend_handler import HandlerPatch
import mplhep as hep
import utility as util
import sys
import pickle, lz4.frame
from itertools import combinations

dfile = uproot.open('/ceph/jhornung/Data_2018/2018/single_muon_data/mm/single_muon_data.root')
dtree = dfile['ntuple'].arrays()

dnpvsGood = dtree['PV_npvsGood'].to_numpy()

dnpvsGood_clipped = np.clip(dnpvsGood, 0, 60.5)

mcfile = uproot.open('/ceph/jhornung/MC_2018/2018/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8+RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2+MINIAODSIM/mm/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8+RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2+MINIAODSIM.root')
mctree = mcfile['ntuple'].arrays()

mcpvsGood = mctree['PV_npvsGood'].to_numpy()
mcpvsGood_clipped = np.clip(mcpvsGood, 0, 60.5)
print(mcpvsGood_clipped.shape)

puweights_nom = mctree['puweight'].to_numpy()
print(puweights_nom.shape)
puweights_up = mctree['puweight_up'].to_numpy()
puweights_down = mctree['puweight_down'].to_numpy()

plt.figure()

dhist = np.histogram(dnpvsGood_clipped, bins=51, range=(0, 51))
dhist_err = np.sqrt(dhist[0])

mchist = np.histogram(mcpvsGood_clipped, weights=puweights_nom, bins=51, range=(0, 51))

plt.hist(mcpvsGood_clipped, weights=puweights_nom/np.max(mchist[0]), bins=51, range=(0, 51), histtype='step', color='blue', label='MC')

bin_centers = 0.5*(dhist[1][1:] + dhist[1][:-1])

plt.errorbar(bin_centers, dhist[0]/np.max(dhist[0]), yerr=0, marker='.', linestyle='none', label='Data', color='black')

plt.show()
