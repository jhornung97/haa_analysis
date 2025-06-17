import sys
import numpy as np
np.set_printoptions(threshold=sys.maxsize)
import matplotlib.pyplot as plt
import uproot
import awkward as ak
sys.path.append("/work/jhornung/Haa/Scripts")
import useful_stuff

ee = uproot.open(sys.argv[1])
mm = uproot.open(sys.argv[2])

ee_ntuples = ee['ntuple']
mm_ntuples = mm['ntuple']

ee_branches = ee_ntuples.arrays()
mm_branches = mm_ntuples.arrays()

files = []
files.append(ee_branches)
files.append(mm_branches)

H_mass_cut_one = (mm_branches['H_mass'] < 110)
H_mass_cut_two = ((mm_branches['H_mass'] >= 110) & (mm_branches['H_mass'] < 140))
H_mass_cut_three  = (mm_branches['H_mass'] >= 140)

masks = []
masks.append(H_mass_cut_one)
masks.append(H_mass_cut_two)
masks.append(H_mass_cut_three)

labels = [r'$m_\mathrm{KKKK} < 110 \mathrm{GeV}$', r'$110 \mathrm{GeV} \leq m_\mathrm{KKKK} < 140 \mathrm{GeV}$', r'$m_\mathrm{KKKK} \geq 140 \mathrm{GeV}$']

fig = plt.subplots()
hadron_delPt = []
pdgid_min_delR_delPt = []

for i, mask in enumerate(masks):

    print(i)

    true_pts = np.array([mm_branches['truth_d1_pt'][mask], mm_branches['truth_d2_pt'][mask], mm_branches['truth_d3_pt'][mask], mm_branches['truth_d4_pt'][mask]]).T    
    true_etas = np.array([mm_branches['truth_d1_eta'][mask], mm_branches['truth_d2_eta'][mask], mm_branches['truth_d3_eta'][mask], mm_branches['truth_d4_eta'][mask]]).T
    true_phis =np.array([mm_branches['truth_d1_phi'][mask], mm_branches['truth_d2_phi'][mask], mm_branches['truth_d3_phi'][mask], mm_branches['truth_d4_phi'][mask]]).T

    pts = mm_branches['nano_pfcands_pt'][mask]
    etas = mm_branches['nano_pfcands_eta'][mask]    
    phis = mm_branches['nano_pfcands_phi'][mask]
    pdgids = mm_branches['nano_pfcands_pdgId'][mask]

    pts_with_match = []

    n = true_etas.shape[0]

    for j in range(n):
        tmp_pts = np.array(pts[j])
        tmp_etas = np.array(etas[j])
        tmp_phis = np.array(phis[j])
        tmp_pdgids = np.array(pdgids[j])
        
        hadron_pts = tmp_pts[np.abs(tmp_pdgids) == 211]
        hadron_etas = tmp_etas[np.abs(tmp_pdgids) == 211]
        hadron_phis = tmp_phis[np.abs(tmp_pdgids) == 211]

        nHadrons = hadron_etas.shape[0]
        nPFCands = tmp_etas.shape[0]

        DeltaEtas = np.array(useful_stuff.PairDiffAk(j, true_etas, hadron_etas))
        AllDeltaEtas = np.array(useful_stuff.PairDiffAk(j, true_etas, tmp_etas))
        corrected_DeltaPhis = np.array(list(map(useful_stuff.CheckPhiDiff, useful_stuff.PairDiffAk(j, true_phis, hadron_phis))))
        AllDeltaPhis = np.array(list(map(useful_stuff.CheckPhiDiff, useful_stuff.PairDiffAk(j, true_phis, tmp_phis))))
        DeltaRs = np.sqrt(np.square(DeltaEtas)+np.square(corrected_DeltaPhis))
        AllDeltaRs = np.sqrt(np.square(AllDeltaEtas)+np.square(AllDeltaPhis))
        
        for k in range(true_pts.shape[1]):
            count = 0            
            for l in range(k*nHadrons, (k+1)*nHadrons, 1):
                if DeltaRs[l] < 0.2:
                    count += 1
            if count >= 1:
                pts_with_match.append(true_pts[j,k])
                minIdx = np.argmin(DeltaRs[k*nHadrons:(k+1)*nHadrons])
                hadron_delPt.append(np.abs(true_pts[j,k] - hadron_pts[minIdx]))
            if count == 0:
                minIdx = np.argmin(AllDeltaRs[k*nPFCands:(k+1)*nPFCands])
                pdgid_min_delR_delPt.append([np.abs(tmp_pdgids[minIdx]), AllDeltaRs[k*nPFCands:(k+1)*nPFCands][minIdx], np.abs(true_pts[j,k] - tmp_pts[minIdx])])

    print(true_pts.flatten().shape)
    pts_with_match = np.array(pts_with_match)
    print(pts_with_match.shape)

    pt_hist, pt_edges = np.histogram(true_pts.flatten(), 40, range=(0,200))
    matched_pt_hist, matched_pt_edges = np.histogram(pts_with_match, 40, range=(0,200))
    pt_match_percentage = matched_pt_hist/pt_hist

    plt.plot(pt_edges[:-1], pt_match_percentage, '-', drawstyle='steps', label=labels[i] + ', ' + str('%.1f'%(100*pts_with_match.shape[0]/true_pts.flatten().shape[0])) + '%')

plt.xlabel(r'$p_\mathrm{T}^{true}$')
plt.ylabel(r'percentage') 
plt.legend()

#plt.savefig('/web/jhornung/public_html/matching/match_higgs_daughters_with_all_hadrons.png') 
#plt.savefig('/web/jhornung/public_html/matching/match_higgs_daughters_with_all_hadrons.pdf') 

fig, ax = plt.subplots()

pdgid_min_delR_delPt = np.array(pdgid_min_delR_delPt)
print(pdgid_min_delR_delPt.shape)

h = ax.hist2d(pdgid_min_delR_delPt[:,0], pdgid_min_delR_delPt[:,1], bins=[10,10], range=[[0,10],[0,0.2]])
fig.colorbar(h[3], ax=ax)
plt.xlabel(r'pdgid')
plt.xticks(np.linspace(0.5, 9.5, 10), np.linspace(0, 9, 10).astype(int))
plt.ylabel(r'$\Delta R$')

#plt.savefig('/web/jhornung/public_html/matching/no_match_pdgid_0_to_10_minDelR.png')
#plt.savefig('/web/jhornung/public_html/matching/no_match_pdgid_0_to_10_minDelR.pdf')

fig = plt.subplots()

hadron_delPt = np.array(hadron_delPt)

plt.hist(hadron_delPt, bins=np.linspace(0, 100, 50), histtype='step')
plt.xlabel(r'$|\Delta p_\mathrm{T}|$')
plt.ylabel(r'count')
plt.savefig('/web/jhornung/public_html/matching/match_delPt.png')
plt.savefig('/web/jhornung/public_html/matching/match_delPt.pdf')

fig, ax = plt.subplots()

h = ax.hist2d(pdgid_min_delR_delPt[:,0], pdgid_min_delR_delPt[:,2], bins=[10,10], range=[[0,10],[0,40]])
fig.colorbar(h[3], ax=ax)
plt.xlabel(r'pdgid')
plt.xticks(np.linspace(0.5, 9.5, 10), np.linspace(0, 9, 10).astype(int))
plt.ylabel(r'$|\Delta p_\mathrm{T}|$')
plt.savefig('/web/jhornung/public_html/matching/no_match_delPt.png')                                                                                                                       
plt.savefig('/web/jhornung/public_html/matching/no_match_delPt.pdf')

plt.show()


