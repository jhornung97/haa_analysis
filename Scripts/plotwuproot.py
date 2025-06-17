import uproot
import math
import numpy as np
import matplotlib.pyplot as plt
#import vector
#import mplhep as hep
#hep.style.use("CMS")

rfile = uproot.open("/work/jhornung/Haa/CROWN/build/bin/output_mm.root")

ntuples = rfile['ntuple']
branches = ntuples.arrays()

# kin of higgs candidate

fig, axs = plt.subplots(nrows=2, ncols=2)
fig.suptitle('Higgs Candidate')

axs[0,0].hist(branches['H_mass'], bins=20, density=True, histtype='step')
axs[0,0].set_xticks(np.arange(0,800,50))
axs[0,0].set_xlabel(r'$m_\mathrm{KKKK}$')
axs[0,0].set_ylabel(r'a.u.')

axs[0,1].hist(branches['H_pt'], bins=20, density=True, histtype='step')
axs[0,1].set_xlabel(r'$p_\mathrm{T}^\mathrm{KKKK}$')
axs[0,1].set_ylabel(r'a.u.')

axs[1,0].hist(branches['H_eta'], bins=20, density=True, histtype='step')
axs[1,0].set_xlabel(r'$\eta_\mathrm{KKKK}$')
axs[1,0].set_ylabel(r'a.u.')

axs[1,1].hist(branches['H_phi'], bins=20, density=True, histtype='step')
axs[1,1].set_xlabel(r'$\phi_\mathrm{KKKK}$')
axs[1,1].set_ylabel(r'a.u.')

fig.tight_layout(pad=0.5)

plt.show()

# kin of kaons

fig, axs = plt.subplots(nrows=1, ncols=3)
fig.suptitle('Hardest Kaon')

axs[0].hist(branches['d1_pt'], bins=20, density=True, histtype='step')
axs[0].set_xlabel(r'$p_\mathrm{T}$')
axs[0].set_ylabel(r'a.u.')

axs[1].hist(branches['d1_eta'], bins=20, density=True, histtype='step')
axs[1].set_xlabel(r'$\eta$')
axs[1].set_ylabel(r'a.u.')

axs[2].hist(branches['d1_phi'], bins=20, density=True, histtype='step')
axs[2].set_xlabel(r'$\phi$')
axs[2].set_ylabel(r'a.u.')

plt.show()

fig, axs = plt.subplots(nrows=1, ncols=3)
fig.suptitle('2nd hardest Kaon')

axs[0].hist(branches['d2_pt'], bins=20, density=True, histtype='step')
axs[0].set_xlabel(r'$p_\mathrm{T}$')
axs[0].set_ylabel(r'a.u.')

axs[1].hist(branches['d2_eta'], bins=20, density=True, histtype='step')
axs[1].set_xlabel(r'$\eta$')
axs[1].set_ylabel(r'a.u.')

axs[2].hist(branches['d2_phi'], bins=20, density=True, histtype='step')
axs[2].set_xlabel(r'$\phi$')
axs[2].set_ylabel(r'a.u.')

plt.show()

fig, axs = plt.subplots(nrows=1, ncols=3)
fig.suptitle('3rd hardest Kaon')

axs[0].hist(branches['d3_pt'], bins=20, density=True, histtype='step')
axs[0].set_xlabel(r'$p_\mathrm{T}$')
axs[0].set_ylabel(r'a.u.')

axs[1].hist(branches['d3_eta'], bins=20, density=True, histtype='step')
axs[1].set_xlabel(r'$\eta$')
axs[1].set_ylabel(r'a.u.')

axs[2].hist(branches['d3_phi'], bins=20, density=True, histtype='step')
axs[2].set_xlabel(r'$\phi$')
axs[2].set_ylabel(r'a.u.')

plt.show()

fig, axs = plt.subplots(nrows=1, ncols=3)
fig.suptitle('4th hardest Kaon')

axs[0].hist(branches['d4_pt'], bins=20, density=True, histtype='step')
axs[0].set_xlabel(r'$p_\mathrm{T}$')
axs[0].set_ylabel(r'a.u.')

axs[1].hist(branches['d4_eta'], bins=20, density=True, histtype='step')
axs[1].set_xlabel(r'$\eta$')
axs[1].set_ylabel(r'a.u.')

axs[2].hist(branches['d4_phi'], bins=20, density=True, histtype='step')
axs[2].set_xlabel(r'$\phi$')
axs[2].set_ylabel(r'a.u.')

plt.show()

# max delta R(K,K)

etas = np.array([branches['d1_eta'], branches['d2_eta'], branches['d3_eta'], branches['d4_eta']]).T
phis = np.array([branches['d1_phi'], branches['d2_phi'], branches['d3_phi'], branches['d4_phi']]).T

nEvents = etas.shape[0]
DeltaEta = []
DeltaPhi = []

for i in range(nEvents):

	EtaPairs = np.array(np.meshgrid(etas[i,:],etas[i,:])).T.reshape(-1,2)
	PhiPairs = np.array(np.meshgrid(phis[i,:],phis[i,:])).T.reshape(-1,2)
	
	EtaDeltas = []
	PhiDeltas = []
	
	for pair in EtaPairs:
	
		delta = pair[0] - pair[1]
		EtaDeltas.append(delta)
	
	for pair in PhiPairs:
	
		delta = pair[0] - pair[1]
		PhiDeltas.append(delta) 
		
	DeltaEta.append(EtaDeltas)
	DeltaPhi.append(PhiDeltas)

DeltaR = np.sqrt(np.square(DeltaEta)+np.square(DeltaPhi))

maxDeltaR = np.max(DeltaR, axis=-1)

fig = plt.figure()

plt.hist(maxDeltaR, bins=20, density=True, histtype='step')
plt.xlabel(r'max $\Delta$R(K,K)')
plt.ylabel(r'a.u.')

plt.show()

# delta phi zh

Delta_Phi_ZH = branches['phi_vis'] - branches['H_phi']

fig = plt.figure()

plt.hist(Delta_Phi_ZH, bins=20, density=True, histtype='step')
plt.xlabel(r'$\Delta\phi_\mathrm{ZH}$')
plt.ylabel(r'a.u.')

plt.show()

# jet multiplicity

fig = plt.figure()

plt.hist(branches['njets'], bins=range(min(branches['njets']), max(branches['njets']+1), 1), histtype='step')
plt.xlabel(r'Jet multiplicity')
plt.ylabel(r'a.u.')

plt.show()
