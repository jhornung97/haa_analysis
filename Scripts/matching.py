import sys
import numpy as np
np.set_printoptions(threshold=sys.maxsize)
import matplotlib.pyplot as plt
import uproot
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

H_mass_cut_one = ((mm_branches['H_mass'] < 110) & (mm_branches['H_mass'] > 1))
H_mass_cut_two = ((mm_branches['H_mass'] >= 110) & (mm_branches['H_mass'] < 140))
H_mass_cut_three  = (mm_branches['H_mass'] >= 140)

masks = []
masks.append(H_mass_cut_one)
masks.append(H_mass_cut_two)
masks.append(H_mass_cut_three)

labels = [r'$m_\mathrm{KKKK} < 110 \mathrm{GeV}$', r'$110 \mathrm{GeV} \leq m_\mathrm{KKKK} < 140 \mathrm{GeV}$', r'$m_\mathrm{KKKK} \geq 140 \mathrm{GeV}$']

print(labels[0])

fig = plt.subplots()

for i, mask in enumerate(masks):

    print(i)

    true_pts = np.array([mm_branches['truth_d1_pt'][mask], mm_branches['truth_d2_pt'][mask], mm_branches['truth_d3_pt'][mask], mm_branches['truth_d4_pt'][mask]]).T
    
    true_etas = np.array([mm_branches['truth_d1_eta'][mask], mm_branches['truth_d2_eta'][mask], mm_branches['truth_d3_eta'][mask], mm_branches['truth_d4_eta'][mask]]).T
    etas = np.array([mm_branches['d1_eta'][mask], mm_branches['d2_eta'][mask], mm_branches['d3_eta'][mask], mm_branches['d4_eta'][mask]]).T

    true_phis =np.array([mm_branches['truth_d1_phi'][mask], mm_branches['truth_d2_phi'][mask], mm_branches['truth_d3_phi'][mask], mm_branches['truth_d4_phi'][mask]]).T
    phis = np.array([mm_branches['d1_phi'][mask], mm_branches['d2_phi'][mask], mm_branches['d3_phi'][mask], mm_branches['d4_phi'][mask]]).T

    DeltaEtas = []
    corrected_DeltaPhis = []
    n = true_etas.shape[0]
    
    for j in range(n):
        DeltaEtas.append(useful_stuff.PairDiff(j, true_etas, etas))
        corrected_DeltaPhi = map(useful_stuff.CheckPhiDiff, useful_stuff.PairDiff(j, true_phis, phis))
        corrected_DeltaPhis.append(list(corrected_DeltaPhi))

    DeltaEtas = np.array(DeltaEtas).flatten()
    corrected_DeltaPhis = np.array(corrected_DeltaPhis).flatten()

    DeltaRs = np.sqrt(np.square(DeltaEtas)+np.square(corrected_DeltaPhis))

    pts_with_match = []

    for j in range(true_pts.flatten().shape[0]):
        count = 0
        for k in range(4*j, 4*(j+1), 1):
            if DeltaRs[k] < 0.2:
                count += 1
        if count >= 1:
            pts_with_match.append(true_pts.flatten()[j])

    print(true_pts.flatten().shape)
    pts_with_match = np.array(pts_with_match)
    print(pts_with_match.shape)

    pt_hist, pt_edges = np.histogram(true_pts.flatten(), 40, range=(0,200))
    matched_pt_hist, matched_pt_edges = np.histogram(pts_with_match, 40, range=(0,200))

    pt_match_percentage = matched_pt_hist/pt_hist

    plt.plot(pt_edges[:-1], pt_match_percentage, '-', drawstyle='steps', label=labels[i]  + ', ' + str('%.1f'%(100*pts_with_match.shape[0]/true_pts.flatten().shape[0])) + '%')

#plt.plot(pt_edges[:-1], np.ones(pt_edges[:-1].shape[0]), color='red')

plt.xlabel(r'$p_\mathrm{T}^{true}$')
plt.ylabel(r'percentage') 
plt.legend()

plt.savefig('/web/jhornung/public_html/matching/match_within_Delta_R_0.2_percentages_charge_pairs.png')
plt.savefig('/web/jhornung/public_html/matching/match_within_Delta_R_0.2_percentages_charge_pairs.pdf')
plt.show()
    
