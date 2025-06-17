import uproot
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.patches as mpatches
import mplhep as hep
from tqdm import tqdm
import seaborn as sns
import awkward as ak
import scipy as sp
import utility
import sys
import pickle, lz4.frame

plt.style.use(hep.style.CMS)   

scope = "mm"

files = []
labels = ['Drell-Yan', r'$t\,\bar{t}$', r'Diboson', r'$W\,H\rightarrow \ell\,\nu\, b\,\bar{b}$', r'Single Top']

idy = uproot.open("/ceph/jhornung/MC_2018/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8+RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2+MINIAODSIM/%s/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8+RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2+MINIAODSIM_final.root"%scope)

idy_ntuples = idy['ntuple']
idy_branches = idy_ntuples.arrays()

ttbar = uproot.open("/ceph/jhornung/MC_2018/ttbar/%s/ttbar.root"%scope)

ttbar_ntuples = ttbar['ntuple']
ttbar_branches = ttbar_ntuples.arrays()

diboson = uproot.open("/ceph/jhornung/MC_2018/diboson/%s/diboson.root"%scope)

diboson_ntuples = diboson['ntuple']
diboson_branches = diboson_ntuples.arrays()

wh = uproot.open("/ceph/jhornung/MC_2018/wh/%s/wh.root"%scope)

wh_ntuples = wh['ntuple']
wh_branches = wh_ntuples.arrays()

singletop = uproot.open("/ceph/jhornung/MC_2018/single_top/%s/single_top.root"%scope)

singletop_ntuples = singletop['ntuple']
singletop_branches = singletop_ntuples.arrays()

signal = uproot.open("/ceph/jhornung/MC_2018/HaaKKKK_M125_13TeV_powheg_pythia8-RunIISummer20UL18MiniAODv2-106X/%s/HaaKKKK_M125_13TeV_powheg_pythia8-RunIISummer20UL18MiniAODv2-106X_final.root"%scope)

signal_ntuples = signal['ntuple']
signal_branches = signal_ntuples.arrays()

files.append(idy_branches)
files.append(ttbar_branches)
files.append(diboson_branches)
files.append(wh_branches)
files.append(singletop_branches)

idy_zpt = idy_branches["gen_pt_vis"].to_numpy()
idy_idx = np.arange(idy_zpt.shape[0])
idy_zpt_weights = np.ones(idy_zpt.shape[0])

eejfile = uproot.open("/work/jhornung/Haa/lindert_qcd_nnlo_sf.root")
eej_histo = eejfile["eej"]

eej_weights = eej_histo.values()
eej_weights_axes = eej_histo.axes
eej_edges = eej_weights_axes[0].edges() 

for i in range(len(eej_edges)-1):
    mask = (idy_zpt >= eej_edges[i]) & (idy_zpt < eej_edges[i+1])
    tmp_idx = idy_idx[mask]
    idy_zpt_weights[tmp_idx] = eej_weights[i]*idy_zpt_weights[tmp_idx]

masks = []
weights = [] 

for i,f in enumerate(files):
    mask = (((f['H_mass'] > 70) & (f['H_mass'] < 110)) | ((f['H_mass'] > 140) & (f['H_mass'] < 200))) & (f["pt_vis"] > 30) & ((f["m_vis"] > 75) & (f["m_vis"] < 105))
    mask = (((f['H_mass'] > 70) & (f['H_mass'] < 110)) | ((f['H_mass'] > 140) & (f['H_mass'] < 200))) & (f['pt_vis'] > 30) & ((f['m_vis'] > 75) & (f['m_vis'] < 105))
    masks.append(mask)
    sfs = f['id_wgt_mu_1'][mask]*f['id_wgt_mu_2'][mask]*f['iso_wgt_mu_1'][mask]*f['iso_wgt_mu_2'][mask]*f['trigger_wgt_mu_1'][mask]*f['trigger_wgt_mu_2'][mask]
    if labels[i] == 'Drell-Yan':
        weight = 59e3*np.array(list(map(utility.weights,f['genWeight'][mask])))*f['evtweight'][mask].to_numpy()*idy_zpt_weights[mask]*sfs#f['puweight'][mask]*sfs
        weights.append(weight)
    else:
        weight = 59e3*np.array(list(map(utility.weights,f['genWeight'][mask])))*f['evtweight'][mask].to_numpy()*sfs#f['puweight'][mask]*sfs
        weights.append(weight)

signal_mask = ((signal_branches['H_mass'] > 70) & (signal_branches['H_mass'] < 200)) & (signal_branches["pt_vis"] > 30) & ((signal_branches["m_vis"] > 75) & (signal_branches["m_vis"] < 105))
signal_sfs = signal_branches['id_wgt_mu_1'][signal_mask]*signal_branches['id_wgt_mu_2'][signal_mask]*signal_branches['iso_wgt_mu_1'][signal_mask]*signal_branches['iso_wgt_mu_2'][signal_mask]*signal_branches['trigger_wgt_mu_1'][signal_mask]*signal_branches['trigger_wgt_mu_2'][signal_mask]
signal_weight = np.array(list(map(utility.weights,signal_branches['genWeight'][signal_mask])))*signal_branches['evtweight'][signal_mask]*signal_sfs

datasets = []

idy_pv = idy_branches["PV_npvsGood"][masks[0]].to_numpy()
ttbar_pv = ttbar_branches["PV_npvsGood"][masks[1]].to_numpy()
diboson_pv = diboson_branches["PV_npvsGood"][masks[2]].to_numpy()
wh_pv = wh_branches["PV_npvsGood"][masks[3]].to_numpy()
singletop_pv = singletop_branches["PV_npvsGood"][masks[4]].to_numpy()

datasets.extend([idy_pv, ttbar_pv, diboson_pv, wh_pv, singletop_pv])

signal_pv = signal_branches["PV_npvsGood"][signal_mask].to_numpy()

all_pv = np.concatenate([idy_pv, ttbar_pv, diboson_pv, wh_pv, singletop_pv, signal_pv])
max_pv = np.max(all_pv)

for i, d in enumerate(datasets):
    hist, edges = np.histogram(d, bins=max_pv+1, range=(0,max_pv+1), weights=weights[i]) 
    event_count = np.sum(hist)
    print(labels[i], event_count)
    labels[i] = labels[i] + ", yield: %.2f"%event_count

fig, ax = plt.subplots()

ax.hist([
    idy_pv,
    ttbar_pv,
    diboson_pv,
    wh_pv,
    singletop_pv
    ], 
    bins=np.linspace(0, max_pv+1, max_pv+1), 
    weights=weights, 
    stacked=True, 
    label=labels, 
    color=['seagreen', 'gold', 'darkorange', 'orange', 'mediumspringgreen']
    )

ax.hist(signal_pv, bins=np.linspace(0, max_pv+1, max_pv+1), weights=signal_weight, histtype='step', color='red', label=r'$H\rightarrow aa\rightarrow KKKK$')

hep.cms.label('', year=2018, lumi=59)

handles, labels = ax.get_legend_handles_labels()
handles.append(mpatches.Patch(color='white', label=r'$\mu\mu$ SR'))
handles.append(mpatches.Patch(color='white', label=r'MC blinded and weighted'))
handles.append(mpatches.Patch(color='white', label=r'without puweights'))
ax.legend(handles=handles, prop={'size': 16})

ax.set_xlabel('Number of good Primary Vertices')
ax.set_xlim(0, max_pv+1)
ax.set_ylabel('Events')
ax.set_ylim(bottom=0)
plt.tight_layout()

plt.savefig('/web/jhornung/public_html/discrepancy_study/good_pv_wo_puweights.png')
plt.savefig('/web/jhornung/public_html/discrepancy_study/good_pv_wo_puweights.pdf')

plt.show()
