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
labels = ['Inclusive Drell-Yan', r'$t\,\bar{t}\rightarrow 2\ell\,2\nu$', r'$W\,W\rightarrow 2\ell\,2\nu$', r'$Z\,Z\rightarrow 2\ell\,2\nu$', r'$W\,Z\rightarrow 3\ell\,\nu$', r'$W\,Z\rightarrow \ell\,3\nu$', r'$W\,H\rightarrow \ell\,\nu\, b\,\bar{b}$', r'Single Top']

idy = uproot.open("/ceph/jhornung/analysis_old/2018/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8+RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2+MINIAODSIM/%s/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8+RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2+MINIAODSIM.root"%scope)

idy_ntuples = idy['ntuple']
idy_branches = idy_ntuples.arrays()

ttbar = uproot.open("/ceph/jhornung/analysis_old/2018/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8+RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1+MINIAODSIM/%s/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8+RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1+MINIAODSIM.root"%scope)

ttbar_ntuples = ttbar['ntuple']
ttbar_branches = ttbar_ntuples.arrays()

ww = uproot.open("/ceph/jhornung/analysis_old/2018/WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8+RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2+MINIAODSIM/%s/WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8+RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2+MINIAODSIM.root"%scope)

ww_ntuples = ww['ntuple']
ww_branches = ww_ntuples.arrays()

zz = uproot.open("/ceph/jhornung/analysis_old/2018/ZZTo2L2Nu_TuneCP5_13TeV_powheg_pythia8+RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1+MINIAODSIM/%s/ZZTo2L2Nu_TuneCP5_13TeV_powheg_pythia8+RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1+MINIAODSIM.root"%scope)

zz_ntuples = zz['ntuple']
zz_branches = zz_ntuples.arrays()

wz_three_l_nu = uproot.open("/ceph/jhornung/analysis_old/2018/WZTo3LNu_mllmin4p0_TuneCP5_13TeV-powheg-pythia8+RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2+MINIAODSIM/%s/WZTo3LNu_mllmin4p0_TuneCP5_13TeV-powheg-pythia8+RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2+MINIAODSIM.root"%scope)

wz_three_l_nu_ntuples = wz_three_l_nu['ntuple']
wz_three_l_nu_branches = wz_three_l_nu_ntuples.arrays()

wz_l_three_nu = uproot.open("/ceph/jhornung/analysis_old/2018/WZTo1L3Nu_4f_TuneCP5_13TeV-amcatnloFXFX-pythia8+RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v3+MINIAODSIM/%s/WZTo1L3Nu_4f_TuneCP5_13TeV-amcatnloFXFX-pythia8+RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v3+MINIAODSIM.root"%scope)

wz_l_three_nu_ntuples = wz_l_three_nu['ntuple']
wz_l_three_nu_branches = wz_l_three_nu_ntuples.arrays()

wh = uproot.open("/ceph/jhornung/analysis_old/2018/merged_mc/%s/wh.root"%scope)

wh_ntuples = wh['ntuple']
wh_branches = wh_ntuples.arrays()

singletop = uproot.open("/ceph/jhornung/analysis_old/2018/merged_mc/%s/singletop.root"%scope)

singletop_ntuples = singletop['ntuple']
singletop_branches = singletop_ntuples.arrays()

signal = uproot.open("/ceph/jhornung/analysis_old/2018/HaaKKKK_M125_13TeV_powheg_pythia8-RunIISummer20UL18MiniAODv2-106X/%s/HaaKKKK_M125_13TeV_powheg_pythia8-RunIISummer20UL18MiniAODv2-106X.root"%scope)

signal_ntuples = signal['ntuple']
signal_branches = signal_ntuples.arrays()

files.append(idy_branches)
files.append(ttbar_branches)
files.append(ww_branches)
files.append(zz_branches)
files.append(wz_three_l_nu_branches)
files.append(wz_l_three_nu_branches)
files.append(wh_branches)
files.append(singletop_branches)

data = uproot.open("/ceph/jhornung/analysis_old/2018/single_muon_data/%s/single_muon_data.root"%scope)

data_ntuples = data['ntuple']
data_branches = data_ntuples.arrays()

weights = []
masks = []

if var == 'H_mass' and scope != 'em':
    for i, f in enumerate(files):
        mask = ((f['H_mass'] > 100) & (f['H_mass'] < 200))
        masks.append(mask)
        weight = 59e3*np.array(list(map(utility.weights,f['genWeight'][mask])))*f['evtweight'][mask].to_numpy()
        weights.append(weight)
        
    signal_mask = (signal_branches['H_mass'] > 100) & (signal_branches['H_mass'] < 200)
    signal_weights = np.array(list(map(utility.weights,signal_branches['genWeight'][signal_mask])))*signal_branches['evtweight'][signal_mask].to_numpy()
    
    idy_var = idy_branches[var][masks[0]].to_numpy()
    ttbar_var = ttbar_branches[var][masks[1]].to_numpy()
    ww_var = ww_branches[var][masks[2]].to_numpy()
    zz_var = zz_branches[var][masks[3]].to_numpy()
    wz_three_l_nu_var = wz_three_l_nu_branches[var][masks[4]].to_numpy()
    wz_l_three_nu_var = wz_l_three_nu_branches[var][masks[5]].to_numpy()
    wh_var = wh_branches[var][masks[6]].to_numpy()
    singletop_var = singletop_branches[var][masks[7]].to_numpy()
    signal_var = signal_branches[var][signal_mask].to_numpy()
    
    data_var = data_branches[var][((data_branches['H_mass'] > 100) & (data_branches['H_mass'] < 110)) | ((data_branches['H_mass'] > 140) & (data_branches['H_mass'] < 200))].to_numpy()
elif var == 'H_mass' and scope == 'em':
    for i, f in enumerate(files):
        weight = 59e3*np.array(list(map(utility.weights,f['genWeight'])))*f['evtweight'].to_numpy()
        weights.append(weight)

    signal_weights = np.array(list(map(utility.weights,signal_branches['genWeight'])))*signal_branches['evtweight'].to_numpy()

    idy_var = idy_branches[var].to_numpy()
    ttbar_var = ttbar_branches[var].to_numpy()
    ww_var = ww_branches[var].to_numpy()
    zz_var = zz_branches[var].to_numpy()
    wz_three_l_nu_var = wz_three_l_nu_branches[var].to_numpy()
    wz_l_three_nu_var = wz_l_three_nu_branches[var].to_numpy()
    wh_var = wh_branches[var].to_numpy()
    singletop_var = singletop_branches[var].to_numpy()
    signal_var = signal_branches[var].to_numpy()

    data_var = data_branches[var].to_numpy()

elif var != 'H_mass' and scope != 'em':
    for i, f in enumerate(files):
        mask = ((f['H_mass'] > 100) & (f['H_mass'] < 110)) | ((f['H_mass'] > 140) & (f['H_mass'] < 200))
        masks.append(mask)
        weight = 59e3*np.array(list(map(utility.weights,f['genWeight'][mask])))*f['evtweight'][mask].to_numpy()
        weights.append(weight)

    signal_weights = np.array(list(map(utility.weights,signal_branches['genWeight'])))*signal_branches['evtweight'].to_numpy()

    idy_var = idy_branches[var][masks[0]].to_numpy()
    ttbar_var = ttbar_branches[var][masks[1]].to_numpy()
    ww_var = ww_branches[var][masks[2]].to_numpy()
    zz_var = zz_branches[var][masks[3]].to_numpy()
    wz_three_l_nu_var = wz_three_l_nu_branches[var][masks[4]].to_numpy()
    wz_l_three_nu_var = wz_l_three_nu_branches[var][masks[5]].to_numpy()
    wh_var = wh_branches[var][masks[6]].to_numpy()
    singletop_var = singletop_branches[var][masks[7]].to_numpy()
    signal_var = signal_branches[var].to_numpy()
    
    data_var = data_branches[var][((data_branches['H_mass'] > 100) & (data_branches['H_mass'] < 110)) | ((data_branches['H_mass'] > 140) & (data_branches['H_mass'] < 200))].to_numpy()

elif var != 'H_mass' and scope == 'em':
    for i, f in enumerate(files):
        mask = (f['H_mass'] > 100) & (f['H_mass'] < 200)
        masks.append(mask)
        weight = 59e3*np.array(list(map(utility.weights,f['genWeight'][mask])))*f['evtweight'][mask].to_numpy()
        weights.append(weight)

    signal_weights = np.array(list(map(utility.weights,signal_branches['genWeight'])))*signal_branches['evtweight'].to_numpy()

    idy_var = idy_branches[var][masks[0]].to_numpy()
    ttbar_var = ttbar_branches[var][masks[1]].to_numpy()
    ww_var = ww_branches[var][masks[2]].to_numpy()
    zz_var = zz_branches[var][masks[3]].to_numpy()
    wz_three_l_nu_var = wz_three_l_nu_branches[var][masks[4]].to_numpy()
    wz_l_three_nu_var = wz_l_three_nu_branches[var][masks[5]].to_numpy()
    wh_var = wh_branches[var][masks[6]].to_numpy()
    singletop_var = singletop_branches[var][masks[7]].to_numpy()
    signal_var = signal_branches[var].to_numpy()
    
    data_var = data_branches[var][(data_branches['H_mass'] > 100) & (data_branches['H_mass'] < 200)].to_numpy()
    
bin_size = (upper_bound - lower_bound)/nbins
fig, axs = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [3, 1]})

mc_data = np.concatenate((idy_var, ttbar_var, ww_var, zz_var, wz_three_l_nu_var, wz_l_three_nu_var, wh_var, singletop_var))
mc_weights = np.concatenate((weights[0], weights[1], weights[2], weights[3], weights[4], weights[5], weights[6], weights[7]))

if var == 'H_mass':
    data_hist, data_edges = np.histogram(
        data_var,
        bins=nbins,
        range=(lower_bound, upper_bound)
    )

    (mc_hist, mc_edges, mc_patches) = axs[0].hist(
        [
            idy_var, 
            ttbar_var, 
            ww_var, 
            zz_var, 
            wz_three_l_nu_var,
            wz_l_three_nu_var, 
            wh_var,
            singletop_var
        ],
        bins=np.linspace(lower_bound, upper_bound, nbins + 1),
        weights=weights,
        stacked=True,
        label=labels,
        color=['seagreen', 'gold', 'darkorange', 'orange', 'magenta', 'mediumspringgreen', 'dodgerblue', 'mediumvioletred']
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
    data_hist, data_edges = np.histogram(
        np.clip(data_var, lower_bound - bin_size/2, upper_bound + bin_size/2),
        bins=nbins + 2,
        range=(lower_bound - bin_size, upper_bound + bin_size),
    )

    (mc_hist, mc_edges, mc_patches) = axs[0].hist(
        [
            np.clip(idy_var, lower_bound - bin_size/2, upper_bound + bin_size/2),
            np.clip(ttbar_var, lower_bound - bin_size/2, upper_bound + bin_size/2),
            np.clip(ww_var, lower_bound - bin_size/2, upper_bound + bin_size/2),
            np.clip(zz_var, lower_bound - bin_size/2, upper_bound + bin_size/2),
            np.clip(wz_three_l_nu_var, lower_bound - bin_size/2, upper_bound + bin_size/2),
            np.clip(wz_l_three_nu_var, lower_bound - bin_size/2, upper_bound + bin_size/2),
            np.clip(wh_var, lower_bound - bin_size/2, upper_bound + bin_size/2),
            np.clip(singletop_var, lower_bound - bin_size/2, upper_bound + bin_size/2)
        ],
        bins=np.linspace(lower_bound - bin_size, upper_bound + bin_size, nbins + 3),
        weights=weights,
        stacked=True,
        label=labels,
        color=['seagreen', 'gold', 'darkorange', 'orange', 'magenta', 'mediumspringgreen', 'dodgerblue', 'mediumvioletred']
    )

    axs[0].hist(
        np.clip(signal_var, lower_bound - bin_size/2, upper_bound + bin_size/2), 
        bins=np.linspace(lower_bound - bin_size, upper_bound + bin_size, nbins + 3), 
        weights=signal_weights, histtype='step', 
        label=r'$H\rightarrow aa\rightarrow KKKK$', 
        color='red'
    )

    signal_hist, signal_edges = np.histogram(
        np.clip(signal_var, lower_bound - bin_size/2, upper_bound + bin_size/2), 
        bins=nbins+2, 
        weights=signal_weights, 
        range=(lower_bound - bin_size, upper_bound + bin_size)
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

print("Scope %s event count:"%scope)
print("Signal: {}".format(np.sum(signal_hist)))
print("Background: {}".format(np.sum(mc_data_hist)))
print("Data: {}".format(np.sum(data_hist)))

bin_centers = 0.5 * (data_edges[:-1] + data_edges[1:])

stat_error = np.sqrt(data_hist[data_hist > 0])

axs[0].errorbar(bin_centers[data_hist > 0], data_hist[data_hist > 0], yerr=stat_error, linestyle='none', marker='.', color='black', label='Single Muon Data')

axs[0].errorbar(bin_centers, signal_hist, yerr=signal_errors, linestyle='none', color='black')

axs[0].errorbar(bin_centers, mc_data_hist, yerr=mc_errors, linestyle='none', color='black')

hep.cms.label('Work in Progress', data=True, year=2018, lumi=59, ax=axs[0])
axs[0].legend()
axs[0].set_yscale('linear')
axs[0].set_ylim(bottom=0)
if var == 'H_mass' or var == 'm_vis' or 'pt' in var: 
    axs[0].set_ylabel('Events/${}\,$GeV'.format(bin_size))
else:
    axs[0].set_ylabel('Events/${}$'.format(bin_size))
ratio = trunc((data_hist[data_hist > 0] - mc_data_hist[data_hist > 0])/mc_data_hist[data_hist > 0], decs=11)
ratio_errors = np.sqrt((stat_error/mc_data_hist[data_hist > 0])**2 + (data_hist[data_hist > 0]*mc_errors[data_hist > 0]/mc_data_hist[data_hist > 0]**2)**2)

axs[1].axhline(0, linestyle=':', color='black')
axs[1].axhline(.25, linestyle=':', color='black')
axs[1].axhline(-.25, linestyle=':', color='black')
axs[1].errorbar(bin_centers[data_hist > 0], ratio, yerr=ratio_errors, linestyle='none', marker='.', color='black')

if var == 'H_mass':
    axs[1].set_xlabel(r'$m_\mathrm{KKKK}$ [GeV]')
elif 'pt' in var:
    axs[1].set_xlabel(r'$p_\mathrm{T}$ [GeV]')
elif 'eta' in var:
    axs[1].set_xlabel(r'$\eta$')
elif 'njets' == var:
    axs[1].set_xlabel(r'Number of Jets')

axs[1].set_ylabel(r'$\frac{\mathrm{Data} - \mathrm{MC}}{\mathrm{MC}}$')
axs[1].set_ylim([-0.2, 0.2])

fig.tight_layout()

plt.savefig('/web/jhornung/public_html/data_to_mc/%s_%s.png'%(var,scope))
plt.savefig('/web/jhornung/public_html/data_to_mc/%s_%s.pdf'%(var,scope)) 

plt.show()
