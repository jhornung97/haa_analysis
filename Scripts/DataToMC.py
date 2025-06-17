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
#sns.set_style("ticks")

def trunc(values, decs=0):
    return np.trunc(values*10**decs)/(10**decs)

plt.rcParams.update({'font.size': 24})

scope = sys.argv[1]
var = sys.argv[2]
lower_bound = int(sys.argv[3])
upper_bound = int(sys.argv[4])
nbins = int(sys.argv[5])

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

data = uproot.open("/ceph/jhornung/Data_2018/single_muon_data/%s/single_muon_data.root"%scope)

data_ntuples = data['ntuple']
data_branches = data_ntuples.arrays()
'''
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

mean = np.mean(theory_weights)

normalized_theory_weights = theory_weights/mean
'''
idy_zpt = idy_branches["gen_pt_vis"].to_numpy()
idy_idx = np.arange(idy_zpt.shape[0])
idy_zpt_weights = np.ones(idy_zpt.shape[0])
'''
for i in range(len(theory_edges)-1):
    mask = (idy_zpt >= theory_edges[i]) & (idy_zpt < theory_edges[i+1])
    tmp_idx = idy_idx[mask]
    idy_zpt_weights[tmp_idx] = normalized_theory_weights[i]

kfactorsfile = uproot.open("/work/jhornung/Haa/merged_kfactors_zjets.root")
kfactors_histo = kfactorsfile["kfactor_monojet_ewk"]
kfactors = kfactors_histo.values()
kfactors_edges = kfactors_histo.axes[0].edges()

for i in range(len(kfactors_edges)-1):
    mask = (idy_zpt >= kfactors_edges[i]) & (idy_zpt < kfactors_edges[i+1])
    tmp_idx = idy_idx[mask]
    idy_zpt_weights[tmp_idx] = kfactors[i]*idy_zpt_weights[tmp_idx]
'''
eejfile = uproot.open("/work/jhornung/Haa/lindert_qcd_nnlo_sf.root")
eej_histo = eejfile["eej"]

eej_weights = eej_histo.values()
eej_weights_axes = eej_histo.axes
eej_edges = eej_weights_axes[0].edges() 

for i in range(len(eej_edges)-1):
    mask = (idy_zpt >= eej_edges[i]) & (idy_zpt < eej_edges[i+1])
    tmp_idx = idy_idx[mask]
    idy_zpt_weights[tmp_idx] = eej_weights[i]*idy_zpt_weights[tmp_idx]

weights = []
masks = []

if var == 'H_mass' and scope != 'em':
    for i, f in enumerate(files):
        mask = (((f['H_mass'] > 70) & (f['H_mass'] < 200))) & (f['pt_vis'] > 30) & ((f['m_vis'] > 75) & (f['m_vis'] < 105)) #(f['H_mass'] < 110)) | ((f['H_mass'] > 140) &
        masks.append(mask)
        if labels[i] == 'Drell-Yan':
            sfs = f['id_wgt_mu_1'][mask]*f['id_wgt_mu_2'][mask]*f['iso_wgt_mu_1'][mask]*f['iso_wgt_mu_2'][mask]*f['trigger_wgt_mu_1'][mask]*f['trigger_wgt_mu_2'][mask]
            weight = 59e3*np.array(list(map(utility.weights,f['genWeight'][mask])))*f['evtweight'][mask].to_numpy()*idy_zpt_weights[mask]*sfs
            weights.append(weight)
        else:
            sfs = f['id_wgt_mu_1'][mask]*f['id_wgt_mu_2'][mask]*f['iso_wgt_mu_1'][mask]*f['iso_wgt_mu_2'][mask]*f['trigger_wgt_mu_1'][mask]*f['trigger_wgt_mu_2'][mask]
            weight = 59e3*np.array(list(map(utility.weights,f['genWeight'][mask])))*f['evtweight'][mask].to_numpy()*sfs
            weights.append(weight)
        
    signal_mask = ((signal_branches['H_mass'] > 70) & (signal_branches['H_mass'] < 200)) & (signal_branches['pt_vis'] > 30) & ((signal_branches['m_vis'] > 75) & (signal_branches['m_vis'] < 105)) 
    signal_sfs = signal_branches['id_wgt_mu_1'][signal_mask]*signal_branches['id_wgt_mu_2'][signal_mask]*signal_branches['iso_wgt_mu_1'][signal_mask]*signal_branches['iso_wgt_mu_2'][signal_mask]*signal_branches['trigger_wgt_mu_1'][signal_mask]*signal_branches['trigger_wgt_mu_2'][signal_mask]
    signal_weights = np.array(list(map(utility.weights,signal_branches['genWeight'][signal_mask])))*signal_branches['evtweight'][signal_mask].to_numpy()*signal_sfs    

    datasets = []

    idy_var = idy_branches[var][masks[0]].to_numpy()
    ttbar_var = ttbar_branches[var][masks[1]].to_numpy()
    diboson_var = diboson_branches[var][masks[2]].to_numpy()
    wh_var = wh_branches[var][masks[3]].to_numpy()
    singletop_var = singletop_branches[var][masks[4]].to_numpy()
    signal_var = signal_branches[var][signal_mask].to_numpy()
    
    datasets.extend([idy_var, ttbar_var, diboson_var, wh_var, singletop_var])

    data_var = data_branches[var][(((data_branches['H_mass'] > 70) & (data_branches['H_mass'] < 110)) | ((data_branches['H_mass'] > 140) & (data_branches['H_mass'] < 200))) & (data_branches['pt_vis'] > 30) & ((data_branches['m_vis'] > 75) & (data_branches['m_vis'] < 105))].to_numpy()
    
elif var == 'H_mass' and scope == 'em':
    for i, f in enumerate(files):
        if labels[i] == 'Drell-Yan':
            mask = (f['pt_vis'] > 30) & ((f['m_vis'] > 75) & (f['m_vis'] < 105))
            masks.append(mask)
            sfs = f['id_wgt_ele_1'][mask]*f['id_wgt_mu_2'][mask]*f['iso_wgt_mu_2'][mask]*f['trigger_wgt_ele_1'][mask]*f['trigger_wgt_mu_2'][mask]
            weight = 59e3*np.array(list(map(utility.weights,f['genWeight'][mask])))*f['evtweight'][mask].to_numpy()*idy_zpt_weights[mask]*sfs
            weights.append(weight)
        else:
            mask = (f['pt_vis'] > 30) & ((f['m_vis'] > 75) & (f['m_vis'] < 105))
            masks.append(mask)
            sfs = f['id_wgt_ele_1'][mask]*f['id_wgt_mu_2'][mask]*f['iso_wgt_mu_2'][mask]*f['trigger_wgt_ele_1'][mask]*f['trigger_wgt_mu_2'][mask]
            weight = 59e3*np.array(list(map(utility.weights,f['genWeight'][mask])))*f['evtweight'][mask].to_numpy()*sfs
            weights.append(weight)
    
    signal_mask = (signal_branches['pt_vis'] > 30) & ((signal_branches['m_vis'] > 75) & (signal_branches['m_vis'] < 105))
    signal_sfs = signal_branches['id_wgt_ele_1'][signal_mask]*signal_branches['id_wgt_mu_2'][signal_mask]*signal_branches['iso_wgt_mu_2'][signal_mask]*signal_branches['trigger_wgt_ele_1'][signal_mask]*signal_branches['trigger_wgt_mu_2'][signal_mask]
    signal_weights = np.array(list(map(utility.weights,signal_branches['genWeight'][signal_mask])))*signal_branches['evtweight'][signal_mask].to_numpy()*signal_sfs

    idy_var = idy_branches[var][masks[0]].to_numpy()
    ttbar_var = ttbar_branches[var][masks[1]].to_numpy()
    diboson_var = diboson_branches[var][masks[2]].to_numpy()
    wh_var = wh_branches[var][masks[3]].to_numpy()
    singletop_var = singletop_branches[var][masks[4]].to_numpy()
    signal_var = signal_branches[var][signal_mask].to_numpy()

    data_var = data_branches[var][(data_branches['pt_vis'] > 30) & ((data_branches['m_vis'] > 75) & (data_branches['m_vis'] < 105))].to_numpy()

elif var != 'H_mass' and scope != 'em':
    for i, f in enumerate(files):
        mask = (((f['H_mass'] > 70) & (f['H_mass'] < 110)) | ((f['H_mass'] > 140) & (f['H_mass'] < 200))) & (f['pt_vis'] > 30) & ((f['m_vis'] > 75) & (f['m_vis'] < 105)) # (f['H_mass'] < 110)) | ((f['H_mass'] > 140) &  
        masks.append(mask)
        if labels[i] == 'Drell-Yan':
            sfs = f['id_wgt_mu_1'][mask]*f['id_wgt_mu_2'][mask]*f['iso_wgt_mu_1'][mask]*f['iso_wgt_mu_2'][mask]*f['trigger_wgt_mu_1'][mask]*f['trigger_wgt_mu_2'][mask]
            weight = 59e3*np.array(list(map(utility.weights,f['genWeight'][mask])))*f['evtweight'][mask].to_numpy()*sfs*idy_zpt_weights[mask]
            weights.append(weight)
        else:
            sfs = f['id_wgt_mu_1'][mask]*f['id_wgt_mu_2'][mask]*f['iso_wgt_mu_1'][mask]*f['iso_wgt_mu_2'][mask]*f['trigger_wgt_mu_1'][mask]*f['trigger_wgt_mu_2'][mask]
            weight = 59e3*np.array(list(map(utility.weights,f['genWeight'][mask])))*f['evtweight'][mask].to_numpy()*sfs
            weights.append(weight)
    
    signal_mask = ((signal_branches['H_mass'] > 70) & (signal_branches['H_mass'] < 200)) & (signal_branches['pt_vis'] > 30) & ((signal_branches['m_vis'] > 75) & (signal_branches['m_vis'] < 105))
    signal_sfs = signal_branches['id_wgt_mu_1'][signal_mask]*signal_branches['id_wgt_mu_2'][signal_mask]*signal_branches['iso_wgt_mu_1'][signal_mask]*signal_branches['iso_wgt_mu_2'][signal_mask]*signal_branches['trigger_wgt_mu_1'][signal_mask]*signal_branches['trigger_wgt_mu_2'][signal_mask]
    signal_weights = 10*np.array(list(map(utility.weights,signal_branches['genWeight'][signal_mask])))*signal_branches['evtweight'][signal_mask].to_numpy()*signal_sfs

    datasets = []

    idy_var = idy_branches[var][masks[0]].to_numpy()
    ttbar_var = ttbar_branches[var][masks[1]].to_numpy()
    #ww_var = ww_branches[var][masks[2]].to_numpy()
    #zz_var = zz_branches[var][masks[3]].to_numpy()
    #wz_three_l_nu_var = wz_three_l_nu_branches[var][masks[4]].to_numpy()
    #wz_l_three_nu_var = wz_l_three_nu_branches[var][masks[5]].to_numpy()
    diboson_var = diboson_branches[var][masks[2]].to_numpy()
    wh_var = wh_branches[var][masks[3]].to_numpy()
    singletop_var = singletop_branches[var][masks[4]].to_numpy()
    signal_var = signal_branches[var][signal_mask].to_numpy()
    
    data_var = data_branches[var][(((data_branches['H_mass'] > 70) & (data_branches['H_mass'] < 110)) | ((data_branches['H_mass'] > 140) & (data_branches['H_mass'] < 200))) & (data_branches['pt_vis'] > 30) & ((data_branches['m_vis'] > 75) & (data_branches['m_vis'] < 105))].to_numpy()

    datasets.extend([idy_var, ttbar_var, diboson_var, wh_var, singletop_var])

elif var != 'H_mass' and scope == 'em':
    for i, f in enumerate(files):
        mask = (((f['H_mass'] > 70) & (f['H_mass'] < 110)) | ((f['H_mass'] > 140) & (f['H_mass'] < 200))) & (f['pt_vis'] > 30) & ((f['m_vis'] > 75) & (f['m_vis'] < 105))
        masks.append(mask)
        if labels[i] == 'Drell-Yan':
            sfs = f['id_wgt_ele_1'][mask]*f['id_wgt_mu_2'][mask]*f['iso_wgt_mu_2'][mask]*f['trigger_wgt_ele_1'][mask]*f['trigger_wgt_mu_2'][mask]
            weight = 59e3*np.array(list(map(utility.weights,f['genWeight'][mask])))*f['evtweight'][mask].to_numpy()*idy_zpt_weights[mask]*sfs
            weights.append(weight)
        else:
            sfs = f['id_wgt_ele_1'][mask]*f['id_wgt_mu_2'][mask]*f['iso_wgt_mu_2'][mask]*f['trigger_wgt_ele_1'][mask]*f['trigger_wgt_mu_2'][mask]
            weight = 59e3*np.array(list(map(utility.weights,f['genWeight'][mask])))*f['evtweight'][mask].to_numpy()*sfs
            weights.append(weight)
            
    signal_mask = ((signal_branches['H_mass'] > 70) & (signal_branches['H_mass'] < 200)) & (signal_branches['pt_vis'] > 30) & ((signal_branches['m_vis'] > 75) & (signal_branches['m_vis'] < 105))
    signal_sfs = signal_branches['id_wgt_ele_1'][signal_mask]*signal_branches['id_wgt_mu_2'][signal_mask]*signal_branches['iso_wgt_mu_2'][signal_mask]*signal_branches['trigger_wgt_ele_1'][signal_mask]*signal_branches['trigger_wgt_mu_2'][signal_mask]
    signal_weights = np.array(list(map(utility.weights,signal_branches['genWeight'][signal_mask])))*signal_branches['evtweight'][signal_mask].to_numpy()

    idy_var = idy_branches[var][masks[0]].to_numpy()
    ttbar_var = ttbar_branches[var][masks[1]].to_numpy()
    diboson_var = diboson_branches[var][masks[2]].to_numpy()
    wh_var = wh_branches[var][masks[3]].to_numpy()
    singletop_var = singletop_branches[var][masks[4]].to_numpy()
    signal_var = signal_branches[var][signal_mask].to_numpy()
    
    data_var = data_branches[var][((data_branches['H_mass'] > 70) & (data_branches['H_mass'] < 200)) & (data_branches['pt_vis'] > 30) & ((data_branches['m_vis'] > 75) & (data_branches['m_vis'] < 105))].to_numpy()
    
bin_size = (upper_bound - lower_bound)/nbins
fig, axs = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [3, 1]})

mc_data = np.concatenate((idy_var, ttbar_var, diboson_var, wh_var, singletop_var)) #ww_var, zz_var, wz_three_l_nu_var, wz_l_three_nu_var,
mc_weights = np.concatenate((weights[0], weights[1], weights[2], weights[3], weights[4])) #weights[5], weights[6], weights[7]))

if var == 'H_mass':

    for i, d in enumerate(datasets):
        hist, edges = np.histogram(d, bins=nbins, range=(lower_bound, upper_bound), weights=weights[i])
        event_sum = np.sum(hist)
        labels[i] = f"{labels[i]}, yield: {event_sum:.2f}"
    data_hist, data_edges = np.histogram(
        data_var,   #mc_data,
        bins=nbins,
        #weights=mc_weights,
        range=(lower_bound, upper_bound)
    )

    data_int = np.sum(data_hist)

    (mc_hist, mc_edges, mc_patches) = axs[0].hist(
        [
            idy_var, 
            ttbar_var, 
            diboson_var,
            wh_var,
            singletop_var
        ],
        bins=np.linspace(lower_bound, upper_bound, nbins + 1),
        weights=weights,
        stacked=True,
        label=labels,
        color=['seagreen', 'gold', 'darkorange', 'orange', 'mediumspringgreen']
    )

    axs[0].hist(
        signal_var, 
        bins=np.linspace(lower_bound, upper_bound, nbins + 1), 
        weights=signal_weights, histtype='step', 
        label=r'$H\rightarrow aa\rightarrow KKKK$', 
        color='red'
    )

    signal_hist, signal_edges = np.histogram(
        signal_var, 
        bins=nbins, 
        weights=signal_weights, 
        range=(lower_bound, upper_bound)
    )
    signal_errors = np.sqrt(
        np.histogram(
            signal_var, 
            bins=nbins, 
            weights=signal_weights**2, 
            range=(lower_bound, upper_bound)
        )[0]
    )

    mc_data_hist, mc_data_edges = np.histogram(
        mc_data, 
        bins=nbins, 
        weights=mc_weights, 
        range=(lower_bound, upper_bound)
    )
    mc_errors = np.sqrt(
        np.histogram(
            mc_data, 
            bins=nbins, 
            weights=mc_weights**2, 
            range=(lower_bound, upper_bound )
        )[0]
    )

elif var != 'H_mass':
    event_sums = []
    for i, d in enumerate(datasets):
        hist, edges = np.histogram(np.clip(d, lower_bound - bin_size/2, upper_bound + bin_size/2), bins=nbins + 2, range=(lower_bound - bin_size, upper_bound + bin_size), weights=weights[i])
        event_sum = np.sum(hist)
        labels[i] = f"{labels[i]}, yield: {event_sum:.2f}"

    data_hist, data_edges = np.histogram(
        np.clip(data_var, lower_bound - bin_size/2, upper_bound + bin_size/2),
        #np.clip(mc_data, lower_bound - bin_size/2, upper_bound + bin_size/2),
        bins=nbins + 2,
        #weights=mc_weights,
        range=(lower_bound - bin_size, upper_bound + bin_size),
    )

    data_int = np.sum(data_hist)

    (mc_hist, mc_edges, mc_patches) = axs[0].hist(
        [
            np.clip(idy_var, lower_bound - bin_size/2, upper_bound + bin_size/2),
            np.clip(ttbar_var, lower_bound - bin_size/2, upper_bound + bin_size/2),
            np.clip(diboson_var, lower_bound - bin_size/2, upper_bound + bin_size/2),
            np.clip(wh_var, lower_bound - bin_size/2, upper_bound + bin_size/2),
            np.clip(singletop_var, lower_bound - bin_size/2, upper_bound + bin_size/2)
        ],
        bins=np.linspace(lower_bound - bin_size, upper_bound + bin_size, nbins + 3),
        weights=weights,
        stacked=True,
        label=labels,
        color=['seagreen', 'gold', 'darkorange', 'orange', 'mediumspringgreen']
    )

    signal_hist, signal_edges = np.histogram(
        np.clip(signal_var, lower_bound - bin_size/2, upper_bound + bin_size/2), 
        bins=nbins+2, 
        weights=signal_weights, 
        range=(lower_bound - bin_size, upper_bound + bin_size)
    )

    axs[0].hist(
        np.clip(signal_var, lower_bound - bin_size/2, upper_bound + bin_size/2), 
        bins=np.linspace(lower_bound - bin_size, upper_bound + bin_size, nbins + 3), 
        weights=signal_weights, histtype='step', 
        label=r'$H\rightarrow aa\rightarrow KKKK$', 
        color='red'
    )
    
    signal_errors = np.sqrt(
        np.histogram(
            np.clip(signal_var, lower_bound - bin_size/2, upper_bound + bin_size/2), 
            bins=nbins + 2, 
            weights=signal_weights**2, 
            range=(lower_bound - bin_size, upper_bound + bin_size)
        )[0]
    )

    mc_data_hist, mc_data_edges = np.histogram(
        np.clip(mc_data, lower_bound - bin_size/2, upper_bound + bin_size/2), 
        bins=nbins + 2, 
        weights=mc_weights, 
        range=(lower_bound - bin_size, upper_bound + bin_size)
    )
    mc_errors = np.sqrt(
        np.histogram(
            np.clip(mc_data, lower_bound - bin_size/2, upper_bound + bin_size/2), 
            bins=nbins + 2, 
            weights=mc_weights**2, 
            range=(lower_bound - bin_size, upper_bound + bin_size)
        )[0]
    )
    

print(bin_size*np.sum(data_hist))
print(bin_size*np.sum(mc_data_hist))

bin_centers = 0.5 * (data_edges[:-1] + data_edges[1:])

stat_error = np.sqrt(data_hist[data_hist > 0])

axs[0].errorbar(bin_centers[data_hist > 0], data_hist[data_hist > 0], yerr=stat_error, linestyle='none', marker='.', color='black', label='Single Muon Data, \nyield: %.2f'%data_int)

hep.cms.label('Work in Progress', data=True, year=2018, lumi=59, ax=axs[0])
handles, labels = axs[0].get_legend_handles_labels()
handles.append(mpatches.Patch(color='white', label=r'$\mu\mu$ SR'))
handles.append(mpatches.Patch(color='white', label=r'Data and MC blinded'))
axs[0].legend(handles=handles, prop={'size': 16})
axs[0].set_yscale('linear')
axs[0].set_ylim(bottom=0)
if var == 'H_mass' or var == 'm_vis' or 'pt' in var: 
    axs[0].set_ylabel('Events/${}\,$GeV'.format(bin_size))
else:
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

if var == 'H_mass':
    axs[1].set_xlabel(r'$m_\mathrm{KKKK}$ [GeV]')
elif var == 'H_eta':
    axs[1].set_xlabel(r'$\eta_\mathrm{KKKK}$')
elif var == 'm_vis':
    axs[1].set_xlabel(r'Dilepton Mass [GeV]')
elif var == 'pt_vis':
    axs[1].set_xlabel(r'Dilepton $p_\mathrm{T}$ [GeV]')
elif 'pt' in var:
    axs[1].set_xlabel(r'$p_\mathrm{T}$ [GeV]')
elif 'eta' in var:
    axs[1].set_xlabel(r'$\eta$')
elif 'njets' == var:
    axs[1].set_xlabel(r'Number of Jets')

axs[1].set_xlim([mc_data_edges[:-1][mc_data_hist > 0][0], mc_data_edges[-1]])

axs[1].set_ylabel(r'$\frac{\mathrm{Data} - \mathrm{MC}}{\mathrm{MC}}$')
axs[1].set_ylim([-0.5, 0.5])
axs[1].set_yticks([-0.25, 0, 0.25])

plt.tight_layout()

fig.subplots_adjust(hspace=0)

#plt.savefig('/web/jhornung/public_html/discrepancy_study/%s_%s_mc_unblinded.png'%(var,scope))
#plt.savefig('/web/jhornung/public_html/discrepancy_study/%s_%s_mc_unblinded.pdf'%(var,scope)) 
#plt.savefig('/web/jhornung/public_html/discrepancy_study/pseudodata_test_%s_%s_mc_unblinded.png'%(var,scope))

plt.show()
