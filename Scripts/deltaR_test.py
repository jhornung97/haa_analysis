import sys
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

signal_ifile = uproot.open("/ceph/jhornung/signal/Haa.root")
signal_tree = signal_ifile["Events"]

dy_ifile = uproot.open("/ceph/jhornung/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8+RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2+MINIAODSIM/014406E7-6906-7041-BBB5-E99625AF42E1.root")
dy_tree = dy_ifile["Events"]

# Load the needed data

signal_genpart_pdgId = signal_tree["GenPart_pdgId"].array()
signal_genpart_motherIdx = signal_tree["GenPart_genPartIdxMother"].array()
signal_genpart_eta = signal_tree["GenPart_eta"].array()
signal_genpart_phi = signal_tree["GenPart_phi"].array()

dy_genpart_pdgId = dy_tree["GenPart_pdgId"].array()
dy_genpart_motherIdx = dy_tree["GenPart_genPartIdxMother"].array()
dy_genpart_eta = dy_tree["GenPart_eta"].array()
dy_genpart_phi = dy_tree["GenPart_phi"].array()

# Select the muons and electrons from the Z boson decay

signal_mask = (((np.abs(signal_genpart_pdgId) == 13) | (np.abs(signal_genpart_pdgId) == 11)) & (np.abs(signal_genpart_pdgId[signal_genpart_motherIdx]) == 23))
dy_mask = (((np.abs(dy_genpart_pdgId) == 13) | (np.abs(dy_genpart_pdgId) == 11)) & (np.abs(dy_genpart_pdgId[dy_genpart_motherIdx]) == 23))

signal_lepton_eta = np.array([eta for eta in signal_genpart_eta[signal_mask] if len(eta) == 2])
signal_lepton_phi = np.array([phi for phi in signal_genpart_phi[signal_mask] if len(phi) == 2])

dy_lepton_eta = np.array([eta for eta in dy_genpart_eta[dy_mask] if len(eta) == 2])
dy_lepton_phi = np.array([phi for phi in dy_genpart_phi[dy_mask] if len(phi) == 2]) 

# Calculate the delta R between the two leptons

signal_deltaEta = signal_lepton_eta[:,0] - signal_lepton_eta[:,1]
signal_deltaPhi = signal_lepton_phi[:,0] - signal_lepton_phi[:,1]

dy_deltaEta = dy_lepton_eta[:,0] - dy_lepton_eta[:,1]
dy_deltaPhi = dy_lepton_phi[:,0] - dy_lepton_phi[:,1]

signal_deltaPhi = util.CheckPhiDiff(signal_deltaPhi)
dy_deltaPhi = util.CheckPhiDiff(dy_deltaPhi)

signal_deltaR = np.sqrt(signal_deltaEta**2 + signal_deltaPhi**2)
dy_deltaR = np.sqrt(dy_deltaEta**2 + dy_deltaPhi**2)

# Plot the delta R

plt.axvline(0.5, color='teal', linestyle='--', label=r'Current $\Delta R$ cut')
plt.axvline(3.0, color='k', linestyle='-.', label=r'Proposed $\Delta R$ cut')

#plt.axvspan(0, 0.5, color='teal', alpha=0.1)
#plt.axvspan(3.0, 6, color='k', alpha=0.1)

plt.xlim(0, 6.5)

plt.hist(dy_deltaR, bins=65, range=(0, 6.5), color="seagreen", density=True, label='Gen. Level DY MC')
plt.hist(signal_deltaR, bins=65, range=(0, 6.5), histtype='step', color="r", density=True, label='Gen. Level Signal MC')

plt.xlabel(r'Lepton $\Delta R$')
plt.ylabel('a.u.')
plt.legend()

hep.cms.label("Work in Progress", data=False, loc=0)

plt.savefig("/web/jhornung/public_html/analysis_plots/genlvl/deltaR_test.png")
plt.savefig("/web/jhornung/public_html/analysis_plots/genlvl/deltaR_test.pdf")

plt.show()