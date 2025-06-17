import uproot
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
#import vector
import mplhep as hep
import awkward as ak
import scipy as sp
hep.style.use("CMS")
#plt.style.use(hep.style.CMS)
plt.rcParams.update({'font.size': 12})

def weights(weight):
        
        if weight > 0:
                weight = 1
        elif weight < 0:
                weight = -1

        return weight

def PairDiff(it, kin):
        
	Pairs = np.array(np.meshgrid(kin[it,:],kin[it,:])).T.reshape(-1,2)
        
	Deltas = []
        
	for pair in Pairs:
                
                delta = pair[0] - pair[1]
                Deltas.append(delta)
                
	return Deltas

def CheckPhiDiff(delPhi):
        
        if delPhi < -np.pi:
                delPhi += 2*np.pi
        if delPhi > np.pi:
                delPhi -= 2*np.pi
        
        return delPhi

files = []
labels = [r'$\mathrm{H}\rightarrow \mathrm{aa}$', 'inclusive Drell-Yan', r'$t\,\bar{t}\rightarrow 2\ell 2\nu$']

#signal = uproot.open("/ceph/jhornung/Haa/2018/add/HaaKKKK_M125_13TeV_powheg_pythia8-RunIISummer20UL18MiniAODv2-106X_mm.root")
signal = uproot.open("/ceph/jhornung/analysis/2018/HaaKKKK_M125_13TeV_powheg_pythia8-RunIISummer20UL18MiniAODv2-106X/mm/HaaKKKK_M125_13TeV_powheg_pythia8-RunIISummer20UL18MiniAODv2-106X.root")
#samples/no_cuts/output_with_hadronic_particles_in_HF/output_mm.root")
signal_ntuples = signal['ntuple']
signal_branches = signal_ntuples.arrays()

idy = uproot.open("/ceph/jhornung/analysis/2018/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8+RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2+MINIAODSIM/mm/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8+RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2+MINIAODSIM.root")

idy_ntuples = idy['ntuple']
idy_branches = idy_ntuples.arrays()

ttbar = uproot.open("/ceph/jhornung/analysis/2018/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8+RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1+MINIAODSIM/mm/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8+RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1+MINIAODSIM.root")

ttbar_ntuples = ttbar['ntuple']
ttbar_branches = ttbar_ntuples.arrays()


files.append(signal_branches)
files.append(idy_branches)
files.append(ttbar_branches)

signs = []
masks = []
deltas = []

for i, f in enumerate(files):
           
        Delta_Phi_ZH = f['phi_vis'] - f['H_phi']
        Corrected_Delta_Phi_ZH = ak.Array(list(map(CheckPhiDiff, Delta_Phi_ZH)))
   
#        mask = (f['H_mass'] >= 50) |
        #Delta_Phi_ZH_cut = (np.abs(Corrected_Delta_Phi_ZH) >= 2.8)
        #cut = four_kaon_syst_pt_cut & Delta_Phi_ZH_cut

 #       masks.append(cut)
        
        signs.append(list(map(weights,f['genWeight'])))
        
# kin of higgs candidate

fig = plt.subplots()

plt.hist(signal_branches['H_mass'], bins=np.linspace(100, 200, 50), weights=.0004*signal_branches['evtweight']*signs[0], histtype='step', label=labels[0])
plt.hist(idy_branches['H_mass'], bins=np.linspace(100, 200, 50), weights=idy_branches['evtweight']*signs[1], histtype='step', label=labels[1])
plt.hist(ttbar_branches['H_mass'], bins=np.linspace(100, 200, 50), weights=ttbar_branches['evtweight']*signs[2], histtype='step', label=labels[2])
        
hep.cms.label('Preliminary',  data=False, year=2018)
#plt.text(176.5, 0.105, 'Higgs Candidate')
plt.legend()
plt.xlabel(r'$m_\mathrm{KKKK}$')
plt.ylabel(r'a.u.')

#plt.savefig('/web/jhornung/public_html/algostudy/H_mass_algo_wo_lepton_filters.png')
#plt.savefig('/web/jhornung/public_html/algostudy/H_mass_algo_wo_lepton_filters.pdf')
plt.show()
'''
fig = plt.subplots(figsize=(8,8))

for i, f in enumerate(files):
        plt.hist(f['H_pt'][masks[i]], bins=np.linspace(0, 250, 50), density=True, weights=f['evtweight'][masks[i]], histtype='step', label=labels[i])

plt.legend()
hep.cms.label('Preliminary', data=False, year=2018)
plt.text(190, 0.0615, 'Higgs Candidate')
plt.xlabel(r'$p_\mathrm{T}^\mathrm{KKKK}$')
plt.ylabel(r'a.u.')
plt.savefig('/web/jhornung/public_html/haa_w_cuts/Higgs_Candidate_pt.png')
plt.savefig('/web/jhornung/public_html/haa_w_cuts/Higgs_Candidate_pt.pdf')

fig = plt.subplots(figsize=(8,8))

for i, f in enumerate(files):
        plt.hist(f['H_eta'][masks[i]], bins=np.linspace(-6, 6, 120), density=True, weights=f['evtweight'][masks[i]], histtype='step', label=labels[i])

plt.legend(loc='lower center')
hep.cms.label('Preliminary', data=False, year=2018)
plt.text(-6, 0.22, 'Higgs Candidate')
plt.xlabel(r'$\eta_\mathrm{KKKK}$')
plt.ylabel(r'a.u.')

plt.savefig('/web/jhornung/public_html/haa_w_cuts/Higgs_Candidate_eta.png')
plt.savefig('/web/jhornung/public_html/haa_w_cuts/Higgs_Candidate_eta.pdf')

fig = plt.subplots(figsize=(8,8))

for i, f in enumerate(files):
        plt.hist(f['H_phi'][masks[i]], bins=50, density=True, weights=f['evtweight'][masks[i]], histtype='step', label=labels[i])

plt.legend(loc='lower center')
hep.cms.label('Preliminary', data=False, year=2018)
plt.text(-.81, 0.0275, 'Higgs Candidate')  
plt.xlabel(r'$\phi$')
plt.ylabel(r'a.u.')

plt.savefig('/web/jhornung/public_html/haa_w_cuts/Higgs_Candidate_phi.png')
plt.savefig('/web/jhornung/public_html/haa_w_cuts/Higgs_Candidate_phi.png')

# kin of kaons
fig = plt.subplots(figsize=(8,8))

for i, f in enumerate(files): 
        plt.hist(f['d1_pt'][masks[i]], bins=np.linspace(0, 150, 50), density=True, weights=f['evtweight'][masks[i]], histtype='step', label=labels[i])

plt.xlabel(r'$p_\mathrm{T}$')
plt.ylabel(r'a.u.')
hep.cms.label('Preliminary', data=False, year=2018)
plt.text(115, 0.1425, 'leading Kaon') 
plt.legend()

plt.savefig('/web/jhornung/public_html/haa_w_cuts/leading_Kaon_pt.png')  
plt.savefig('/web/jhornung/public_html/haa_w_cuts/leading_Kaon_pt.pdf')

fig = plt.subplots(figsize=(8,8))                                                                                                                                                                           
for i, f in enumerate(files):
        plt.hist(f['d1_eta'][masks[i]], bins=50, density=True, weights=f['evtweight'][masks[i]], histtype='step', label=labels[i])

plt.xlabel(r'$\eta$')
plt.ylabel(r'a.u.')
hep.cms.label('Preliminary', data=False, year=2018)
plt.text(-3.2, 0.32, 'leading Kaon')
plt.legend(loc='lower center')

plt.savefig('/web/jhornung/public_html/haa_w_cuts/leading_Kaon_eta.png')
plt.savefig('/web/jhornung/public_html/haa_w_cuts/leading_Kaon_eta.pdf')

fig = plt.subplots(figsize=(8,8))                                                                                                                                                                         

for i, f in enumerate(files):
        plt.hist(f['d1_phi'][masks[i]], bins=50, density=True, weights=f['evtweight'][masks[i]], histtype='step', label=labels[i])

plt.xlabel(r'$\phi$')
plt.ylabel(r'a.u.')
hep.cms.label('Preliminary', data=False, year=2018)
plt.text(-.81, 0.0275, 'leading Kaon')
plt.legend(loc='lower center')

plt.savefig('/web/jhornung/public_html/haa_w_cuts/leading_Kaon_phi.png')
plt.savefig('/web/jhornung/public_html/haa_w_cuts/leading_Kaon_phi.pdf')


fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(18,6))

for j, f in enumerate(files):
        axs[i].hist(f['d2_' + kin], bins=50, density=True, weights=signs[j]*f['puweight'], histtype='step', label=labels[j])

axs[0].set_xlabel(r'$p_\mathrm{T}$')
axs[0].set_ylabel(r'a.u.')
hep.cms.label('Preliminary', ax=axs[0], data=False, year=2018)
axs[0].legend()

axs[1].set_xlabel(r'$\eta$')
axs[1].set_ylabel(r'a.u.')
hep.cms.label('Preliminary', ax=axs[1], data=False, year=2018)
axs[1].legend()

axs[2].set_xlabel(r'$\phi$')
axs[2].set_ylabel(r'a.u.')
hep.cms.label('Preliminary', ax=axs[2], data=False, year=2018)
axs[2].legend(loc='lower right')
fig.tight_layout()

plt.savefig('/web/jhornung/public_html/K_2.png')

fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(18,6))

for j, f in enumerate(files):
        axs[i].hist(f['d3_' + kin], bins=50, density=True, weights=signs[j]*f['puweight'], histtype='step', label=labels[j])

axs[0].set_xlabel(r'$p_\mathrm{T}$')
axs[0].set_ylabel(r'a.u.')
hep.cms.label('Preliminary', ax=axs[0], data=False, year=2018)
axs[0].legend()

axs[1].set_xlabel(r'$\eta$')
axs[1].set_ylabel(r'a.u.')
hep.cms.label('Preliminary', ax=axs[1], data=False, year=2018)
axs[1].legend()

axs[2].set_xlabel(r'$\phi$')
axs[2].set_ylabel(r'a.u.')
hep.cms.label('Preliminary', ax=axs[2], data=False, year=2018)
axs[2].legend(loc='lower right')
fig.tight_layout()

plt.savefig('/web/jhornung/public_html/K_3.png')

fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(18,6))

for j, f in enumerate(files):
        axs[i].hist(f['d4_' + kin], bins=50, density=True, weights=signs[j]*f['puweight'], histtype='step', label=labels[j])

axs[0].set_xlabel(r'$p_\mathrm{T}$')
axs[0].set_ylabel(r'a.u.')
hep.cms.label('Preliminary', ax=axs[0], data=False, year=2018)
axs[0].legend()

axs[1].set_xlabel(r'$\eta$')
axs[1].set_ylabel(r'a.u.')
hep.cms.label('Preliminary', ax=axs[1], data=False, year=2018)
axs[1].legend()

axs[2].set_xlabel(r'$\phi$')
axs[2].set_ylabel(r'a.u.')
hep.cms.label('Preliminary', ax=axs[2], data=False, year=2018)
axs[2].legend(loc='lower right')
fig.tight_layout()

plt.savefig('/web/jhornung/public_html/K_4.png')

# max delta R(K,K)
fig = plt.figure(figsize=(8,8))
for i, f in enumerate(files):
        etas = np.array([f['d1_eta'], f['d2_eta'], f['d3_eta'], f['d4_eta']]).T
        phis = np.array([f['d1_phi'], f['d2_phi'], f['d3_phi'], f['d4_phi']]).T
        
        nEvents = etas.shape[0]
        DeltaEtas = []
        DeltaPhis = []

        for j in range(nEvents):
                
                DeltaEtas.append(PairDiff(j,etas))
                DeltaPhis.append(list(map(CheckPhiDiff,PairDiff(j,phis))))
                        

        DeltaR = np.sqrt(np.square(DeltaEtas)+np.square(DeltaPhis))

        maxDeltaR = ak.Array(np.max(np.abs(DeltaR), axis=-1))

        plt.hist(maxDeltaR[masks[i]], bins=50, density=True, weights=f['evtweight'][masks[i]], histtype='step', label=labels[i])
        
plt.xlabel(r'max $\Delta$R(K,K)')
plt.ylabel(r'a.u.')
hep.cms.label('Preliminary', data=False, year=2018)
plt.legend()

plt.savefig('/web/jhornung/public_html/haa_w_cuts/Max_Delta_R.png')
plt.savefig('/web/jhornung/public_html/haa_w_cuts/Max_Delta_R.pdf')

# delta phi zh

fig = plt.figure(figsize=(8,8))

for i, f in enumerate(files):
        Delta_Phi_ZH = f['phi_vis'] - f['H_phi']
        Corrected_Delta_Phi_ZH = ak.Array(list(map(CheckPhiDiff, Delta_Phi_ZH)))
        plt.hist(Corrected_Delta_Phi_ZH[masks[i]], bins=50, density=True, weights=f['evtweight'][masks[i]], histtype='step', label=labels[i])

plt.xlabel(r'$\Delta\phi_\mathrm{ZH}$')
plt.ylabel(r'a.u.')
hep.cms.label('Preliminary', data=False, year=2018)
plt.legend(loc='upper center')

plt.savefig('/web/jhornung/public_html/haa_w_cuts/Delta_phi_ZH.png')
plt.savefig('/web/jhornung/public_html/haa_w_cuts/Delta_phi_ZH.pdf')

# dilepton system
fig = plt.subplots(figsize=(8,8))

for i, f in enumerate(files):
        plt.hist(f['m_vis'][masks[i]], bins=np.linspace(40, 140, 50), density=True, weights=f['evtweight'][masks[i]], histtype='step', label=labels[i])

plt.xlabel(r'$m_{\ell\ell}$')
plt.ylabel(r'a.u.')
hep.cms.label('Preliminary', data=False, year=2018)
plt.text(116, 0.1125, 'Z candidate')
plt.legend()

plt.savefig('/web/jhornung/public_html/haa_w_cuts/dilepton_mass.png')
plt.savefig('/web/jhornung/public_html/haa_w_cuts/dilepton_mass.pdf')

# Jet multiplicity                                                                                                                                                                                         

fig = plt.figure(figsize=(8,8))
for i,f in enumerate(files):
        plt.hist(f['njets'][masks[i]], bins=range(min(f['njets']), max(f['njets']+1), 1), density=True, weights=f['evtweight'][masks[i]], histtype='step', label=labels[i])

plt.xlabel(r'Jet multiplicity')
plt.ylabel(r'a.u.')
hep.cms.label('Preliminary', data=False, year=2018)
plt.legend()

plt.savefig('/web/jhornung/public_html/haa_w_cuts/Jet_Multiplicity.png')
plt.savefig('/web/jhornung/public_html/haa_w_cuts/Jet_Multiplicity.pdf')

# leading jet vars
fig = plt.subplots(figsize=(8,8))

for i,f in enumerate(files):
        plt.hist(f['jpt_1'][masks[i]], bins=np.linspace(0, 400, 50), density=True, weights=f['evtweight'][masks[i]], histtype='step', label=labels[i])

plt.xlabel(r'$p_\mathrm{T}$')
plt.ylabel(r'a.u.')
hep.cms.label('Preliminary', data=False, year=2018)
plt.text(305, 0.0365, 'leading Jet')
plt.legend()

plt.savefig('/web/jhornung/public_html/haa_w_cuts/leading_jet_pt.png')
plt.savefig('/web/jhornung/public_html/haa_w_cuts/leading_jet_pt.pdf')

fig = plt.subplots(figsize=(8,8))

for i,f in enumerate(files):
        plt.hist(f['jeta_1'][masks[i]], bins=np.linspace(-4, 4, 50), density=True, weights=f['evtweight'][masks[i]], histtype='step', label=labels[i])

plt.xlabel(r'$\eta$')
plt.ylabel(r'a.u.')
hep.cms.label('Preliminary', data=False, year=2018)
plt.text(-4, 0.31, 'leading Jet')
plt.legend(loc='lower center')

plt.savefig('/web/jhornung/public_html/haa_w_cuts/leading_jet_eta.png')
plt.savefig('/web/jhornung/public_html/haa_w_cuts/leading_jet_eta.pdf')

fig = plt.subplots(figsize=(8,8))

for i,f in enumerate(files):
        plt.hist(f['jbtagDeepB_1'][masks[i]], bins=np.linspace(0, 1, 50), density=True, weights=f['evtweight'][masks[i]], histtype='step', label=labels[i])

plt.xlabel(r'b tag value')
plt.ylabel(r'a.u.')
hep.cms.label('Preliminary', data=False, year=2018)
plt.text(0.76, 20.25, 'leading Jet')
plt.legend()

plt.savefig('/web/jhornung/public_html/haa_w_cuts/leading_jet_btag.png')
plt.savefig('/web/jhornung/public_html/haa_w_cuts/leading_jet_btag.pdf')
'''
