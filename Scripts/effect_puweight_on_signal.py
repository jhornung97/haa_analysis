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
from itertools import combinations
import HistLib as hl

plt.style.use(hep.style.CMS)
plt.rcParams.update({'font.size': 14})

mc_loc = "/ceph/jhornung/MC_2018/2018"

scope = "mm"

signal = ["HaaKKKK_M125_13TeV_powheg_pythia8-RunIISummer20UL18MiniAODv2-106X", "DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8+RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2+MINIAODSIM"]
signal_labels = [r'signal', r'DY']#, r"2017"]
mc_signal = hl.data(mc_loc, signal, signal_labels, scope)
signal_branches = mc_signal.branches_dict
'''
branches = signal_branches[r'$H\rightarrow aa \rightarrow 4\,K$, $m_a = 1.5\, GeV$']

comb = list(combinations(["12", "14", "23", "34"], 2))
ps_mass_mask = np.zeros(len(branches['ps_mass_12']), dtype=bool)
for c in comb:
    ps_mass_1 = branches[f'ps_mass_{c[0]}'].to_numpy()
    ps_mass_2 = branches[f'ps_mass_{c[1]}'].to_numpy()
    ps_mass_mask |= (ps_mass_1 < 10) & (ps_mass_2 < 10) & (np.abs(ps_mass_1 - ps_mass_2) < 1)

h_mass_mask = (branches['H_mass'] > 70) & (branches['H_mass'] < 200)
h_eta_mask = (np.abs(branches['H_eta']) < 2.4)
z_pt_mask = (branches['pt_vis'] > 30)
z_mass_mask = (branches['m_vis'] > 75) & (branches['m_vis'] < 105)
mask = h_mass_mask & h_eta_mask & z_pt_mask & z_mass_mask & ps_mass_mask

sfs = branches['id_wgt_mu_1'][mask]*branches['id_wgt_mu_2'][mask]*branches['iso_wgt_mu_1'][mask]*branches['iso_wgt_mu_2'][mask]*branches['trigger_wgt_mu_1'][mask]*branches['trigger_wgt_mu_2'][mask]
'''

for k in signal_branches.keys():
    branches = signal_branches[k]
    comb = list(combinations(["12", "14", "23", "34"], 2))
    ps_mass_mask = np.zeros(len(branches['ps_mass_12']), dtype=bool)
    for c in comb:
        ps_mass_1 = branches[f'ps_mass_{c[0]}'].to_numpy()
        ps_mass_2 = branches[f'ps_mass_{c[1]}'].to_numpy()
        ps_mass_mask |= (ps_mass_1 < 10) & (ps_mass_2 < 10) & (np.abs(ps_mass_1 - ps_mass_2) < 1)

    if k == 'DY':
        h_mass_mask = ((branches['H_mass'] > 70) & (branches['H_mass'] < 110)) | ((branches['H_mass'] > 140) & (branches['H_mass'] < 200))
    else:
        h_mass_mask = (branches['H_mass'] > 70) & (branches['H_mass'] < 200)
    h_eta_mask = (np.abs(branches['H_eta']) < 2.4)
    leading_k_mask = (branches['d1_pt'] > 10) & (branches['d2_pt'] > 10)
    z_pt_mask = (branches['pt_vis'] > 30)
    z_mass_mask = (branches['m_vis'] > 75) & (branches['m_vis'] < 105)
    
    mask = h_mass_mask & h_eta_mask & leading_k_mask & z_pt_mask & z_mass_mask & ps_mass_mask

    zpt = branches['gen_pt_vis'].to_numpy()
    idy_idx = np.arange(zpt.shape[0])
    zpt_weights = np.ones(zpt.shape[0])

    if k == 'DY':
        eejfile = uproot.open("/work/jhornung/Haa/lindert_qcd_nnlo_sf.root")
        eej_histo = eejfile["eej"]
        eej_weights = eej_histo.values()
        eej_weights_axes = eej_histo.axes
        eej_edges = eej_weights_axes[0].edges() 
        for i in range(len(eej_edges)-1):
            idxmask = (zpt >= eej_edges[i]) & (zpt < eej_edges[i+1])
            tmp_idx = idy_idx[idxmask]
            zpt_weights[tmp_idx] = 59e3*eej_weights[i]*zpt_weights[tmp_idx]


    sfs = branches['id_wgt_mu_1'][mask]*branches['id_wgt_mu_2'][mask]*branches['iso_wgt_mu_1'][mask]*branches['iso_wgt_mu_2'][mask]*branches['trigger_wgt_mu_1'][mask]*branches['trigger_wgt_mu_2'][mask]
    weights = branches['evtweight'][mask]*zpt_weights[mask]*sfs*np.array(list(map(lambda x: 1 if x > 0 else -1, branches['genWeight'][mask])))
    hist_wo_puweight, edges = np.histogram(branches["PV_npvsGood"][mask], bins=np.linspace(0, 50, 50), weights=weights)
    print(hist_wo_puweight)
    hist_w_puweight, _ = np.histogram(branches["PV_npvsGood"][mask], bins=np.linspace(0, 50, 50), weights=weights*branches['puweight'][mask])
    print(hist_w_puweight)
    bin_centers = (edges[:-1] + edges[1:]) / 2
    mean_wo_puweight = np.average(bin_centers, weights=hist_wo_puweight)
    mean_w_puweight = np.average(bin_centers, weights=hist_w_puweight)
    delta_mean = mean_w_puweight - mean_wo_puweight
    fig, ax = plt.subplots()
    ax.hist(branches["PV_npvsGood"][mask], bins=np.linspace(0, 50, 50), histtype='step', weights=weights, label=rf"w/o puweight, $\mu={mean_wo_puweight:.2f}$")
    ax.hist(branches["PV_npvsGood"][mask], bins=np.linspace(0, 50, 50), histtype='step', weights=weights*branches['puweight'][mask], label=rf"w puweight, $\mu={mean_w_puweight:.2f}$")
    ax.legend()
    plt.figtext(0.15, 0.8, rf"$\Delta\mu={delta_mean:.2f}$")
    ax.set_xlabel('Number of good primary vertices')
    ax.set_ylabel('Events')
    plt.savefig(f"/web/jhornung/public_html/pustudy/Haa_analysis_puweight_{k}.png")
    plt.savefig(f"/web/jhornung/public_html/pustudy/Haa_analysis_puweight_{k}.pdf")
    #ax.set_title(k)
    fig, ax = plt.subplots()
    branches = signal_branches[k]
    ax.scatter(branches["PV_npvsGood"][mask], branches['puweight'][mask], label=k)
    ax.set_xlabel('Number of good primary vertices')
    ax.set_ylabel('PU weight')
    plt.savefig(f"/web/jhornung/public_html/pustudy/Haa_analysis_puweight_scatter_{k}.png")

mc_loc = "/ceph/jhornung/test_puweights"
signal = ["DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8+RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2+MINIAODSIM"]
signal_labels = [r'DY']#, r"2017"]

mc_signal = hl.data(mc_loc, signal, signal_labels, scope)

signal_branches = mc_signal.branches_dict

for k in signal_branches.keys():
    branches = signal_branches[k]
    weights = np.array(list(map(lambda x: 1 if x > 0 else -1, branches['genWeight'])))
    hist_wo_puweight, edges = np.histogram(branches["PV_npvsGood"], bins=np.linspace(0, 50, 50), weights=weights)
    hist_w_puweight, _ = np.histogram(branches["PV_npvsGood"], bins=np.linspace(0, 50, 50), weights=weights*branches['puweight'])
    bin_centers = (edges[:-1] + edges[1:]) / 2
    mean_wo_puweight = np.average(bin_centers, weights=hist_wo_puweight)
    mean_w_puweight = np.average(bin_centers, weights=hist_w_puweight)
    delta_mean = mean_w_puweight - mean_wo_puweight
    
    fig, ax = plt.subplots()
    ax.hist(branches["PV_npvsGood"], bins=np.linspace(0, 50, 50), histtype='step', weights=weights, label=rf"w/o puweight, $\mu={mean_wo_puweight:.2f}$")
    ax.hist(branches["PV_npvsGood"], bins=np.linspace(0, 50, 50), histtype='step', weights=weights*branches['puweight'], label=rf"w puweight, $\mu={mean_w_puweight:.2f}$")
    ax.legend()
    plt.figtext(0.15, 0.8, rf"$\Delta\mu={delta_mean:.2f}$")
    ax.set_xlabel('Number of good primary vertices')
    ax.set_ylabel('Events')
    plt.savefig(f"/web/jhornung/public_html/pustudy/test_puweight_{k}.png")
    plt.savefig(f"/web/jhornung/public_html/pustudy/test_puweight_{k}.pdf")
    fig, ax = plt.subplots()
    ax.scatter(branches["PV_npvsGood"], branches['puweight'], label=k)
    ax.set_xlabel('Number of good primary vertices')
    ax.set_ylabel('PU weight')
    plt.savefig(f"/web/jhornung/public_html/pustudy/test_puweight_scatter_{k}.png")
    plt.savefig(f"/web/jhornung/public_html/pustudy/test_puweight_scatter_{k}.pdf")


fig, ax = plt.subplots()

samplefile = uproot.open("/work/jhornung/Haa/KingMaker/CROWN/data/pileup/Data_Pileup_2018_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18.root")
hist = samplefile["pileup"]

vals = hist.values()
edges = hist.axes[0].edges()
errors = hist.errors()
bin_centers = (edges[:-1] + edges[1:]) / 2

ax.errorbar(bin_centers, vals, yerr=errors, fmt='o', label='2018')
ax.set_xlabel('Number of PU interactions')
ax.set_xlim(0, 25)
plt.savefig("/web/jhornung/public_html/pustudy/pileup_hist_2018.png")

sfile = uproot.open("/ceph/jhornung/signal/Haa.root")
sbranches = sfile["Events"].arrays()

fig, ax = plt.subplots()

ax.hist(sbranches["Pileup_nTrueInt"], bins=np.linspace(0, 76, 77), histtype='step', label='2018')
ax.set_xlabel('Number of true PU interactions')

plt.savefig("/web/jhornung/public_html/pustudy/nTrueInt_signal.png")

#plt.show()