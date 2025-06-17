import sys
import numpy as np
np.set_printoptions(threshold=sys.maxsize)
import matplotlib.pyplot as plt
import uproot
import awkward as ak
sys.path.append("/work/jhornung/Haa/Scripts")
import useful_stuff

#ee = uproot.open(sys.argv[1])
mm = uproot.open(sys.argv[1])

#ee_ntuples = ee['ntuple']
mm_ntuples = mm['ntuple']

#ee_branches = ee_ntuples.arrays()
mm_branches = mm_ntuples.arrays()

h_pts = mm_branches['truth_h_pt'].to_numpy()
#h_pts = np.append(h_pts, mm_branches['truth_h_pt'].to_numpy())

etas = np.array([mm_branches['truth_d1_eta'], mm_branches['truth_d2_eta'], mm_branches['truth_d3_eta'], mm_branches['truth_d4_eta']]).T
phis = np.array([mm_branches['truth_d1_phi'], mm_branches['truth_d2_phi'], mm_branches['truth_d3_phi'], mm_branches['truth_d4_phi']]).T

#ee_etas = np.array([ee_branches['truth_d1_eta'], ee_branches['truth_d2_eta'], ee_branches['truth_d3_eta'], ee_branches['truth_d4_eta']]).T
#ee_phis = np.array([ee_branches['truth_d1_phi'], ee_branches['truth_d2_phi'], ee_branches['truth_d3_phi'], ee_branches['truth_d4_phi']]).T

#etas = np.append(mm_etas, ee_etas, axis=0)
#phis = np.append(mm_phis, ee_phis, axis=0)

DeltaRs = []                                                                                                                                                                                          
TrueDeltaRs = []                                                                                                                                                                                           
FalseDeltaRs = []

'''
rootfile = uproot.open(sys.argv[1])
events = rootfile['Events']
branches = events.arrays()

nGenParts = branches['nGenPart']
genpts = branches['GenPart_pt']
genetas = branches['GenPart_eta']
genphis = branches['GenPart_phi']
genpdgids = branches['GenPart_pdgId']
genmotheridx = branches['GenPart_genPartIdxMother']

hpts = []
DeltaRs = []
TrueDeltaRs = []
FalseDeltaRs = []

for i, genpt in enumerate(genpts):
    n = nGenParts[i]
    tmp_pts = np.array(genpt)
    tmp_etas = np.array(genetas[i])
    tmp_phis = np.array(genphis[i])
    tmp_pdgids = np.array(genpdgids[i])
    tmp_mothers = np.array(genmotheridx[i])
    H_pt_already_in = False
    for i in range(n):
        if tmp_pdgids[tmp_mothers[i]] == 25 & H_pt_already_in == False:
            hpts.append(tmp_pts[tmp_mothers[i]])
            H_pt_already_in = True

hpts = np.array(hpts)
print(hpts.shape)

'''
for i in range(h_pts.shape[0]):
    DeltaRpos = np.round(np.sqrt((etas[i,0] - etas[i,1])**2 + useful_stuff.CheckPhiDiff(phis[i,0] - phis[i,1])**2), 4)
    DeltaRneg = np.round(np.sqrt((etas[i,2] - etas[i,3])**2 + useful_stuff.CheckPhiDiff(phis[i,2] - phis[i,3])**2), 4)
    TrueDeltaRs.append([DeltaRpos, DeltaRneg])
    DeltaEtas = useful_stuff.PairDiff(i, etas, etas)
    DeltaPhis = np.array(list(map(useful_stuff.CheckPhiDiff, useful_stuff.PairDiff(i, phis, phis))))
    EventDeltaRs = np.round(np.sqrt(np.square(DeltaEtas) + np.square(DeltaPhis)), 4)
    DeltaRs.append(EventDeltaRs)
    FalseDeltaRs.append(EventDeltaRs[(EventDeltaRs != DeltaRpos) & (EventDeltaRs != DeltaRneg)])
    if EventDeltaRs[(EventDeltaRs != DeltaRpos) & (EventDeltaRs != DeltaRneg)].shape[0] != 4:
        print(EventDeltaRs[(EventDeltaRs != DeltaRpos) & (EventDeltaRs != DeltaRneg)].shape[0])
        print([DeltaRpos, DeltaRneg])
        print(EventDeltaRs[(EventDeltaRs != DeltaRpos) & (EventDeltaRs != DeltaRneg)])
        
DeltaRs = np.array(DeltaRs)
meanDeltaRs = np.mean(DeltaRs, axis=-1)

fig, ax = plt.subplots()

h = ax.hist2d(h_pts, meanDeltaRs, range=[[0,450], [0,4]])
fig.colorbar(h[3], ax=ax)
plt.xlabel(r'$p_\mathrm{T}^\mathrm{H}$')
plt.ylabel(r'mean $\Delta R_\mathrm{KK}$')
plt.savefig('/web/jhornung/public_html/algostudy/Delta_R_means.png')
plt.savefig('/web/jhornung/public_html/algostudy/Delta_R_means.pdf')

TrueDeltaRs = np.array(TrueDeltaRs)
TrueDeltaRMeans = np.mean(TrueDeltaRs, axis=-1)

fig, ax = plt.subplots()

h = ax.hist2d(h_pts, TrueDeltaRMeans, range=[[0,450], [0,0.15]])
fig.colorbar(h[3], ax=ax)
plt.xlabel(r'$p_\mathrm{T}^\mathrm{H}$')
plt.ylabel(r'mean $\Delta R_\mathrm{KK}^\mathrm{true}$')
plt.savefig('/web/jhornung/public_html/algostudy/Delta_R_means_True_Combinations.png')
plt.savefig('/web/jhornung/public_html/algostudy/Delta_R_means_True_Combinations.pdf')

FalseDeltaRs = np.array(FalseDeltaRs)
FalseDeltaRMeans = np.mean(FalseDeltaRs, axis=-1)

fig, ax = plt.subplots()

h = ax.hist2d(h_pts, FalseDeltaRMeans, range=[[0,450], [0,4]])
fig.colorbar(h[3], ax=ax)
plt.xlabel(r'$p_\mathrm{T}^\mathrm{H}$')
plt.ylabel(r'mean $\Delta R_\mathrm{KK}^\mathrm{false}$')
plt.savefig('/web/jhornung/public_html/algostudy/Delta_R_means_False_Combinations.png')
plt.savefig('/web/jhornung/public_html/algostudy/Delta_R_means_False_Combinations.pdf')

plt.show()

