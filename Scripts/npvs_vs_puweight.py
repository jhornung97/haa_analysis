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
data_loc = "/ceph/jhornung/Data_2018/2018"
scope = "mm"

bkgs = [
        "DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8+RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2+MINIAODSIM"
        ]
bkg_labels = [
        'Drell-Yan'
        ]
mc_bkg = hl.data(mc_loc, bkgs, bkg_labels, scope)
bkg_branches = mc_bkg.branches_dict
z_pt_weights = mc_bkg.get_Z_pt_weights()

fig, ax = plt.subplots()

plt.scatter(bkg_branches['Drell-Yan']['PV_npvs'], bkg_branches['Drell-Yan']['puweight'])
plt.xlabel('PV_npvs')
plt.ylabel('puweight')
hep.cms.label(data=False, year='2018', ax=ax)

#plt.savefig("/web/jhornung/public_html/analysis_plots/2018/npu_vs_puweight.png")
#plt.savefig("/web/jhornung/public_html/analysis_plots/2018/npu_vs_puweight.pdf")

fig, ax = plt.subplots()

mean = np.mean(bkg_branches['Drell-Yan']['puweight'])

normalized = bkg_branches['Drell-Yan']['puweight'] / mean

plt.scatter(bkg_branches['Drell-Yan']['PV_npvs'], normalized)
plt.xlabel('PV_npvs')
plt.ylabel('puweight/mean(puweight)')
hep.cms.label(data=False, year='2018', ax=ax)
#plt.savefig("/web/jhornung/public_html/analysis_plots/2018/npu_vs_puweight_normalized.png")
#plt.savefig("/web/jhornung/public_html/analysis_plots/2018/npu_vs_puweight_normalized.pdf")
plt.show()