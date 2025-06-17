import uproot
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import mplhep as hep
import seaborn as sns
import awkward as ak
import scipy as sp
import utility
import sys
import pickle, lz4.frame

plt.style.use(hep.style.CMS)                                                                       
#sns.set_style("ticks")

def trunc(values, decs=0):
    return np.trunc(values*10**decs)/(10**decs)
                                                                                                               
plt.rcParams.update({'font.size': 24})

scope = sys.argv[1]
nbins = int(sys.argv[2])

files = []
labels = ['Drell-Yan', r'$t\,\bar{t}$', r'Diboson', r'$W\,H\rightarrow \ell\,\nu\, b\,\bar{b}$', r'Single Top']

idy = uproot.open("/ceph/jhornung/MC_2018/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8+RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2+MINIAODSIM/%s/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnl\
oFXFX-pythia8+RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2+MINIAODSIM_final.root"%scope)

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

data = uproot.open("/ceph/jhornung/Data_2018/single_muon_data/%s/single_muon_data.root"%scope)

data_ntuples = data['ntuple']
data_branches = data_ntuples.arrays()

res = pickle.load(lz4.frame.open("/work/jhornung/Haa/scetlib_dyturboCorrZ.pkl.lz4"))
h = res["Z"]["scetlib_dyturbo_hist"][{"vars": "pdf0"}]

theory_edges = h.axes[2].edges

theory_zpt_hist = h.project("qT").values()
integral = np.sum(theory_zpt_hist)
theory_dens = theory_zpt_hist/integral

idy_reweighting = uproot.open("/ceph/jhornung/reweighting/reweighting.root")
idy_reweighting_ntuples = idy_reweighting['ntuple']
idy_reweighting_branches = idy_reweighting_ntuples.arrays()

idy_gen_pt_hist, idy_gen_pt_edges = np.histogram(idy_reweighting_branches["z_pt"].to_numpy(),
                                                 bins=theory_edges,
                                                 weights=np.array(list(map(utility.weights,idy_reweighting_branches['genWeight']))),
                                                 density=True
)

theory_weights = theory_dens/idy_gen_pt_hist

idy_zpt = idy_branches["gen_pt_vis"].to_numpy()
idy_idx = np.arange(idy_zpt.shape[0])
idy_zpt_weights = np.ones(idy_zpt.shape[0])

for i in range(len(theory_edges)-1):
    mask = (idy_zpt >= theory_edges[i]) & (idy_zpt < theory_edges[i+1])
    tmp_idx = idy_idx[mask]
    idy_zpt_weights[tmp_idx] = theory_weights[i]

reweightingfile = uproot.open("/work/jhornung/Haa/lindert_qcd_nnlo_sf.root")
eej_histo = reweightingfile["eej"]

eej_weights = eej_histo.values()
eej_weights_axes = eej_histo.axes
eej_edges = eej_weights_axes[0].edges()

for i in range(len(eej_edges)-1):
    mask = (idy_zpt >= eej_edges[i]) & (idy_zpt < eej_edges[i+1])
    tmp_idx = idy_idx[mask]
    idy_zpt_weights[tmp_idx] = eej_weights[i]*idy_zpt_weights[tmp_idx]

weights = []
masks = []
    
if scope == 'em':
    for i, f in enumerate(files):
        if labels[i] == 'Drell-Yan':
            mask = (((f['H_mass'] > 70) & (f['H_mass'] < 200)) & (f['pt_vis'] > 30) & ((f['m_vis'] > 75) & (f['m_vis'] < 105)))
            masks.append(mask)
            sfs = f['id_wgt_ele_1'][mask]*f['id_wgt_mu_2'][mask]*f['iso_wgt_mu_2'][mask]*f['trigger_wgt_ele_1'][mask]*f['trigger_wgt_mu_2'][mask]
            weight = 59e3*np.array(list(map(utility.weights,f['genWeight'][mask])))*f['evtweight'][mask].to_numpy()*idy_zpt_weights[mask]*sfs
            weights.append(weight)
        else:
            mask = (((f['H_mass'] > 70) & (f['H_mass'] < 200)) & (f['pt_vis'] > 30) & ((f['m_vis'] > 75) & (f['m_vis'] < 105)))
            masks.append(mask)
            sfs = f['id_wgt_ele_1'][mask]*f['id_wgt_mu_2'][mask]*f['iso_wgt_mu_2'][mask]*f['trigger_wgt_ele_1'][mask]*f['trigger_wgt_mu_2'][mask]
            weight = 59e3*np.array(list(map(utility.weights,f['genWeight'][mask])))*f['evtweight'][mask].to_numpy()*sfs
            weights.append(weight)
    
    signal_mask = ((signal_branches['H_mass'] > 70) & (signal_branches['H_mass'] < 200)) & (signal_branches['pt_vis'] > 30) & ((signal_branches['m_vis'] > 75) & (signal_branches['m_vis'] < 105))
    signal_sfs = signal_branches['id_wgt_ele_1'][signal_mask]*signal_branches['id_wgt_mu_2'][signal_mask]*signal_branches['iso_wgt_mu_2'][signal_mask]*signal_branches['trigger_wgt_ele_1'][signal_mask]*signal_branches['trigger_wgt_mu_2'][signal_mask]
    signal_weights = np.array(list(map(utility.weights,signal_branches['genWeight'][signal_mask])))*signal_branches['evtweight'][signal_mask].to_numpy()*signal_sfs

    idy_H_phi = idy_branches['H_phi'][masks[0]].to_numpy()
    ttbar_H_phi = ttbar_branches['H_phi'][masks[1]].to_numpy()
    diboson_H_phi = diboson_branches['H_phi'][masks[2]].to_numpy()
    wh_H_phi = wh_branches['H_phi'][masks[3]].to_numpy()
    singletop_H_phi = singletop_branches['H_phi'][masks[4]].to_numpy()
    signal_H_phi = signal_branches['H_phi'][signal_mask].to_numpy()

    data_H_phi = data_branches['H_phi'][((data_branches['H_mass'] > 70) & (data_branches['H_mass'] < 200)) & ((data_branches['m_vis'] > 75) & (data_branches['m_vis'] < 105)) & (data_branches['pt_vis'] > 30)].to_numpy()

    idy_phi_vis = idy_branches['phi_vis'][masks[0]].to_numpy()
    ttbar_phi_vis = ttbar_branches['phi_vis'][masks[1]].to_numpy()
    diboson_phi_vis = diboson_branches['phi_vis'][masks[2]].to_numpy()
    wh_phi_vis = wh_branches['phi_vis'][masks[3]].to_numpy()
    singletop_phi_vis = singletop_branches['phi_vis'][masks[4]].to_numpy()
    signal_phi_vis = signal_branches['phi_vis'][signal_mask].to_numpy()

    data_phi_vis = data_branches['phi_vis'][((data_branches['H_mass'] > 70) & (data_branches['H_mass'] < 200)) & ((data_branches['m_vis'] > 75) & (data_branches['m_vis'] < 105)) & (data_branches['pt_vis'] > 30)].to_numpy()

    idy_delta_phi = np.array(list(map(utility.CheckPhiDiff, idy_H_phi - idy_phi_vis)))
    ttbar_delta_phi = np.array(list(map(utility.CheckPhiDiff, ttbar_H_phi - ttbar_phi_vis)))
    diboson_delta_phi = np.array(list(map(utility.CheckPhiDiff, diboson_H_phi - diboson_phi_vis)))
    wh_delta_phi = np.array(list(map(utility.CheckPhiDiff, wh_H_phi - wh_phi_vis)))
    singletop_delta_phi = np.array(list(map(utility.CheckPhiDiff, singletop_H_phi - singletop_phi_vis)))
    signal_delta_phi = np.array(list(map(utility.CheckPhiDiff, signal_H_phi - signal_phi_vis)))

    data_delta_phi = np.array(list(map(utility.CheckPhiDiff, data_H_phi - data_phi_vis)))

elif scope != 'em':
    for i, f in enumerate(files):
        if labels[i] == 'Drell-Yan':
            mask = (((f['H_mass'] > 70) & (f['H_mass'] < 200))) & (f['pt_vis'] > 30) & ((f['m_vis'] > 75) & (f['m_vis'] < 105))
            masks.append(mask)
            sfs = f['id_wgt_mu_1'][mask]*f['id_wgt_mu_2'][mask]*f['iso_wgt_mu_1'][mask]*f['iso_wgt_mu_2'][mask]*f['trigger_wgt_mu_1'][mask]*f['trigger_wgt_mu_2'][mask]
            weight = 59e3*np.array(list(map(utility.weights,f['genWeight'][mask])))*f['evtweight'][mask].to_numpy()*idy_zpt_weights[mask]*sfs
            weights.append(weight)
        else:
            mask = (((f['H_mass'] > 70) & (f['H_mass'] < 200)) & (f['pt_vis'] > 30) & ((f['m_vis'] > 75) & (f['m_vis'] < 105)))
            masks.append(mask)
            sfs = f['id_wgt_mu_1'][mask]*f['id_wgt_mu_2'][mask]*f['iso_wgt_mu_1'][mask]*f['iso_wgt_mu_2'][mask]*f['trigger_wgt_mu_1'][mask]*f['trigger_wgt_mu_2'][mask]
            weight = 59e3*np.array(list(map(utility.weights,f['genWeight'][mask])))*f['evtweight'][mask].to_numpy()*sfs
            weights.append(weight)

    signal_mask = (((signal_branches['H_mass'] > 70) & (signal_branches['H_mass'] < 110)) | ((signal_branches['H_mass'] > 140) & (signal_branches['H_mass'] < 200))) & (signal_branches['pt_vis'] > 30) & ((signal_branches['m_vis'] > 75) & (signal_branches['m_vis'] < 105))
    signal_sfs = signal_branches['id_wgt_mu_1'][signal_mask]*signal_branches['id_wgt_mu_2'][signal_mask]*signal_branches['iso_wgt_mu_1'][signal_mask]*signal_branches['iso_wgt_mu_2'][signal_mask]*signal_branches['trigger_wgt_mu_1'][signal_mask]*signal_branches['trigger_wgt_mu_2'][signal_mask]
    signal_weights = np.array(list(map(utility.weights,signal_branches['genWeight'][signal_mask])))*signal_branches['evtweight'][signal_mask].to_numpy()*signal_sfs

    idy_H_phi = idy_branches['H_phi'][masks[0]].to_numpy()
    ttbar_H_phi = ttbar_branches['H_phi'][masks[1]].to_numpy()
    diboson_H_phi = diboson_branches['H_phi'][masks[2]].to_numpy()
    wh_H_phi = wh_branches['H_phi'][masks[3]].to_numpy()
    singletop_H_phi = singletop_branches['H_phi'][masks[4]].to_numpy()
    signal_H_phi = signal_branches['H_phi'][signal_mask].to_numpy()

    data_H_phi = data_branches['H_phi'][(((data_branches['H_mass'] > 70) & (data_branches['H_mass'] < 110)) | ((data_branches['H_mass'] > 140)  & (data_branches['H_mass'] < 200))) & ((data_branches['m_vis'] > 75) & (data_branches['m_vis'] < 105)) & (data_branches['pt_vis'] > 30)].to_numpy()

    idy_phi_vis = idy_branches['phi_vis'][masks[0]].to_numpy()
    ttbar_phi_vis = ttbar_branches['phi_vis'][masks[1]].to_numpy()
    diboson_phi_vis = diboson_branches['phi_vis'][masks[2]].to_numpy()
    wh_phi_vis = wh_branches['phi_vis'][masks[3]].to_numpy()
    singletop_phi_vis = singletop_branches['phi_vis'][masks[4]].to_numpy()
    signal_phi_vis = signal_branches['phi_vis'][signal_mask].to_numpy()

    data_phi_vis = data_branches['phi_vis'][(((data_branches['H_mass'] > 70) & (data_branches['H_mass'] < 110)) | ((data_branches['H_mass'] > 140)  & (data_branches['H_mass'] < 200))) & ((data_branches['m_vis'] > 75) & (data_branches['m_vis'] < 105)) & (data_branches['pt_vis'] > 30)].to_numpy()

    idy_delta_phi = np.array(list(map(utility.CheckPhiDiff, idy_H_phi - idy_phi_vis)))
    ttbar_delta_phi = np.array(list(map(utility.CheckPhiDiff, ttbar_H_phi - ttbar_phi_vis)))
    diboson_delta_phi = np.array(list(map(utility.CheckPhiDiff, diboson_H_phi - diboson_phi_vis)))
    wh_delta_phi = np.array(list(map(utility.CheckPhiDiff, wh_H_phi - wh_phi_vis)))
    singletop_delta_phi = np.array(list(map(utility.CheckPhiDiff, singletop_H_phi - singletop_phi_vis)))
    signal_delta_phi = np.array(list(map(utility.CheckPhiDiff, signal_H_phi - signal_phi_vis)))

    data_delta_phi = np.array(list(map(utility.CheckPhiDiff, data_H_phi - data_phi_vis)))

bin_size = (3.2 - (-3.2))/nbins
fig, axs = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [3, 1]})

mc_data = np.concatenate((idy_delta_phi, ttbar_delta_phi, diboson_delta_phi, wh_delta_phi, singletop_delta_phi))
mc_weights = np.concatenate((weights[0], weights[1], weights[2], weights[3], weights[4]))

data_hist, data_edges = np.histogram(
    data_delta_phi,
    bins=nbins,
    range=(-3.2, 3.2)
)

(mc_hist, mc_edges, mc_patches) = axs[0].hist(
    [
        idy_delta_phi, 
        ttbar_delta_phi, 
        diboson_delta_phi, 
        wh_delta_phi, 
        singletop_delta_phi
    ],
    bins=np.linspace(-3.2, 3.2, nbins + 1),
    weights=weights,
    stacked=True,
    label=labels,
    color=['seagreen', 'gold', 'darkorange', 'orange', 'mediumspringgreen']
)

axs[0].hist(
    signal_delta_phi, 
    bins=np.linspace(-3.2, 3.2, nbins + 1), 
    weights=signal_weights, histtype='step', 
    label=r'$H\rightarrow aa\rightarrow KKKK$', 
    color='red'
)

signal_hist, signal_edges = np.histogram(
    signal_delta_phi, 
    bins=nbins, 
    weights=signal_weights, 
    range=(-3.2, 3.2)
)
signal_errors = np.sqrt(
    np.histogram(
        signal_delta_phi, 
        bins=nbins, 
        weights=signal_weights**2, 
        range=(-3.2, 3.2)
    )[0]
)

mc_data_hist, mc_data_edges = np.histogram(
    mc_data, 
    bins=nbins, 
    weights=mc_weights, 
    range=(-3.2, 3.2)
)
mc_errors = np.sqrt(
    np.histogram(
        mc_data, 
        bins=nbins, 
        weights=mc_weights**2, 
        range=(-3.2, 3.2)
    )[0]
)

print(np.sum(data_hist))
print(np.sum(mc_data_hist))

bin_centers = 0.5 * (signal_edges[:-1] + signal_edges[1:])

stat_error = np.sqrt(data_hist[data_hist > 0])

axs[0].errorbar(bin_centers[data_hist > 0], data_hist[data_hist > 0], yerr=stat_error, linestyle='none', marker='.', color='black', label='Single Muon Data')

axs[0].errorbar(bin_centers, signal_hist, yerr=signal_errors, linestyle='none', color='black')

hep.cms.label('Work in Progress', data=True, year=2018, lumi=59, ax=axs[0])
axs[0].legend()
axs[0].set_yscale('linear')
axs[0].set_ylim(bottom=0)
axs[0].set_ylabel('Events/${}$'.format(bin_size))

ratio = trunc((data_hist[data_hist > 0] - mc_data_hist[data_hist > 0])/mc_data_hist[data_hist > 0], decs=11)
stat_ratio_errors = stat_error/mc_data_hist[data_hist > 0]
mc_ratio_errors = data_hist[data_hist > 0]*mc_errors[data_hist > 0]/mc_data_hist[data_hist > 0]**2

axs[1].axhline(0, linestyle=':', color='black')
axs[1].axhline(.25, linestyle=':', color='black')
axs[1].axhline(-.25, linestyle=':', color='black')
axs[1].errorbar(bin_centers[data_hist > 0], ratio, yerr=stat_ratio_errors, linestyle='none', marker='.', color='black')
axs[1].fill_between(bin_centers[data_hist > 0], mc_ratio_errors, -mc_ratio_errors, step='mid', facecolor='gray', edgecolor='white', alpha=0.2, label="MC Error")
axs[1].legend(prop={'size': 12})

axs[1].set_xlabel(r'$\Delta\phi_\mathrm{ZH}$')
axs[1].set_xlim([mc_data_edges[:-1][mc_data_hist > 0][0], mc_data_edges[-1]])

axs[1].set_ylabel(r'$\frac{\mathrm{Data} - \mathrm{MC}}{\mathrm{MC}}$')
axs[1].set_ylim([-0.5, 0.5])
axs[1].set_yticks([-0.25, 0, 0.25])

fig.tight_layout()
fig.subplots_adjust(hspace=0)

plt.savefig('/web/jhornung/public_html/redo_data_to_mc/Delta_Phi_%s_with_sf.png'%(scope))
plt.savefig('/web/jhornung/public_html/redo_data_to_mc/Delta_Phi_%s_with_sf.pdf'%(scope)) 

plt.show()
