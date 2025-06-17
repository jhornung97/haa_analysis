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

H_mass_cut = (mm_branches['H_mass'] > 0)

true_pts = np.array([mm_branches['truth_d1_pt'][H_mass_cut], mm_branches['truth_d2_pt'][H_mass_cut], mm_branches['truth_d3_pt'][H_mass_cut], mm_branches['truth_d4_pt'][H_mass_cut]]).T
pts = np.array([mm_branches['d1_pt'][H_mass_cut], mm_branches['d2_pt'][H_mass_cut], mm_branches['d3_pt'][H_mass_cut], mm_branches['d4_pt'][H_mass_cut]]).T
    
true_etas = np.array([mm_branches['truth_d1_eta'][H_mass_cut], mm_branches['truth_d2_eta'][H_mass_cut], mm_branches['truth_d3_eta'][H_mass_cut], mm_branches['truth_d4_eta'][H_mass_cut]]).T
etas = np.array([mm_branches['d1_eta'][H_mass_cut], mm_branches['d2_eta'][H_mass_cut], mm_branches['d3_eta'][H_mass_cut], mm_branches['d4_eta'][H_mass_cut]]).T

true_phis =np.array([mm_branches['truth_d1_phi'][H_mass_cut], mm_branches['truth_d2_phi'][H_mass_cut], mm_branches['truth_d3_phi'][H_mass_cut], mm_branches['truth_d4_phi'][H_mass_cut]]).T
phis = np.array([mm_branches['d1_phi'][H_mass_cut], mm_branches['d2_phi'][H_mass_cut], mm_branches['d3_phi'][H_mass_cut], mm_branches['d4_phi'][H_mass_cut]]).T

print(true_pts.flatten().shape)

DeltaEtas = []
corrected_DeltaPhis = []
pTs = []

n = true_etas.shape[0]

for i in range(n):
    DeltaEtas.append(useful_stuff.PairDiff(i, true_etas, etas))
    corrected_DeltaPhi = map(useful_stuff.CheckPhiDiff, useful_stuff.PairDiff(i, true_phis, phis))
    corrected_DeltaPhis.append(list(corrected_DeltaPhi))
#    pTs.append()


DeltaEtas = np.array(DeltaEtas)

corrected_DeltaPhis = np.array(corrected_DeltaPhis)

DeltaRs = np.sqrt(np.square(DeltaEtas)+np.square(corrected_DeltaPhis))

DelR_cut = (DeltaRs < 0.2)
#print(DelR_cut)
check = []

for i in range(4): 
    psum = np.sum(DelR_cut[:,(4*i):(4*(i+1))], axis=-1)
    pcheck = (psum >= 1)
#    print(psum)
    check.append(pcheck)

check = np.array(check).T.flatten()

print(true_pts.flatten()[check].shape)

pt_hist, pt_edges = np.histogram(true_pts.flatten(), 40, range=(0,200))
matched_pt_hist, matched_pt_edges = np.histogram(true_pts.flatten()[check], 40, range=(0,200))

pt_match_percentage = matched_pt_hist/pt_hist

print(pt_edges)

fig = plt.subplots()
plt.plot(pt_edges[:-1], pt_match_percentage, '-', drawstyle='steps')
#plt.xticks(pt_edges)
plt.ylim(0,1.1)
#plt.show()

'''
n, bins, patches = plt.hist(DeltaRs.flatten()[DelR_cut], bins=np.linspace(0, 6, 60), histtype='step')
percent = np.sum(n)/DeltaRs.flatten().shape[0]

plt.xlabel(r'$\Delta R$')
plt.ylabel(r'count')
plt.text(4.25, 250, r'$m_\mathrm{KKKK} < 125$GeV')
plt.text(4.25, 240, r'$\Delta R > 0.2$')
plt.text(4.25, 230, r'Percentage of wrongly') 
plt.text(4.25, 220, r'identified daughters: ' + str('%.3f'%(percent)))

plt.savefig('/web/jhornung/public_html/matching/DeltaR_mKKKK_less_than_125.png')
plt.savefig('/web/jhornung/public_html/matching/DeltaR_mKKKK_less_than_125.pdf')

for i in range(4):
    fig = plt.subplots(figsize=(8,8))
    
    DelR_cut = (DeltaRs[:,i] > 0.2)

    plt.hist(DeltaRs[:,i][DelR_cut], bins=np.linspace(0, 6, 60), histtype='step')
    plt.xlabel(r'$\Delta R$')
    plt.ylabel(r'count')
    if i == 0:
        plt.text(4.75, 50, r'Pair ' + str(i+1))
        plt.text(4.75, 48.5, r'$m_\mathrm{KKKK} < 125$GeV')
        plt.text(4.75, 47, r'$\Delta R > 0.2$')
    if i == 1:
        plt.text(4.75, 70, r'Pair ' + str(i+1))
        plt.text(4.75, 68, r'$m_\mathrm{KKKK} < 125$GeV')
        plt.text(4.75, 66, r'$\Delta R > 0.2$')
    if i == 2:
        plt.text(4.75, 82, r'Pair ' + str(i+1))
        plt.text(4.75, 80, r'$m_\mathrm{KKKK} < 125$GeV')
        plt.text(4.75, 77.5, r'$\Delta R > 0.2$')
    if i == 3:
        plt.text(4.75, 82, r'Pair ' + str(i+1))
        plt.text(4.75, 80, r'$m_\mathrm{KKKK} < 125$GeV')
        plt.text(4.75, 77.5, r'$\Delta R > 0.2$')
#    plt.show()
    plt.savefig('/web/jhornung/public_html/matching/DeltaR_of_Pair'+ str(i+1) + '_mKKKK_less_than_125.png')
    plt.savefig('/web/jhornung/public_html/matching/DeltaR_of_Pair'+ str(i+1) + '_mKKKK_less_than_125.pdf')
    
    fig = plt.subplots(figsize=(8,8))

    x = np.linspace(0,np.max(np.maximum(true_pts[:,i][DelR_cut], pts[:,i][DelR_cut])), 50)
    plt.plot(x, x, color='red')

    plt.scatter(true_pts[:,i][DelR_cut], pts[:,i][DelR_cut], marker='x')
    plt.xlabel(r'$p_\mathrm{T}^\mathrm{truth}$')
    plt.ylabel(r'$p_\mathrm{T}^\mathrm{reco}$')
    if i == 0:
        plt.ylim(0, np.max(pts[:,i][DelR_cut]) + 2.5)
        plt.text(0, 160, r'Pair ' + str(i+1))
        plt.text(0, 155, r'$m_\mathrm{KKKK} < 125$GeV')
        plt.text(0, 150, r'$\Delta R > 0.2$')
    if i == 1:
        plt.ylim(0, np.max(pts[:,i][DelR_cut]) + 2.5)
        plt.text(0, 150, r'Pair ' + str(i+1))
        plt.text(0, 145, r'$m_\mathrm{KKKK} < 125$GeV')
        plt.text(0, 140, r'$\Delta R > 0.2$')
    if i == 2:
        plt.ylim(0, np.max(pts[:,i][DelR_cut]) + 2.5)
        plt.text(0, 70, r'Pair ' + str(i+1))
        plt.text(0, 68, r'$m_\mathrm{KKKK} < 125$GeV')
        plt.text(0, 66, r'$\Delta R > 0.2$')
    if i == 3:
        plt.ylim(0, np.max(pts[:,i][DelR_cut]) + 2.5)
        plt.text(0, 42.5, r'Pair ' + str(i+1))
        plt.text(0, 41.5, r'$m_\mathrm{KKKK} < 125$GeV')
        plt.text(0, 40, r'$\Delta R > 0.2$')

#    plt.show()
    
    plt.savefig('/web/jhornung/public_html/matching/pt_scatter_Pair'+ str(i+1) + '.png')
    plt.savefig('/web/jhornung/public_html/matching/pt_scatter_Pair'+ str(i+1) + '.pdf')

    fig = plt.subplots(figsize=(8,8))

    x = np.linspace(np.min(np.minimum(true_etas[:,i][DelR_cut], etas[:,i][DelR_cut])),np.max(np.maximum(true_etas[:,i][DelR_cut], etas[:,i][DelR_cut])), 50)
    plt.plot(x, x, color='red')

    plt.scatter(true_etas[:,i][DelR_cut], etas[:,i][DelR_cut], marker='x')
    plt.xlabel(r'$\eta^\mathrm{truth}$')
    plt.ylabel(r'$\eta^\mathrm{reco}$')
    if i == 0:
        plt.ylim(np.min(etas[:,i][DelR_cut] - .25), np.max(etas[:,i][DelR_cut]) + .25)
        plt.text(0, 160, r'Pair ' + str(i+1))
        plt.text(0, 155, r'$m_\mathrm{KKKK} < 125$GeV')
        plt.text(0, 150, r'$\Delta R > 0.2$')
    if i == 1:
        plt.ylim(np.min(etas[:,i][DelR_cut] - .25), np.max(etas[:,i][DelR_cut]) + .25)
        plt.text(0, 150, r'Pair ' + str(i+1))
        plt.text(0, 145, r'$m_\mathrm{KKKK} < 125$GeV')
        plt.text(0, 140, r'$\Delta R > 0.2$')
    if i == 2:
        plt.ylim(np.min(etas[:,i][DelR_cut] - .25), np.max(etas[:,i][DelR_cut]) + .25)
        plt.text(0, 70, r'Pair ' + str(i+1))
        plt.text(0, 68, r'$m_\mathrm{KKKK} < 125$GeV')
        plt.text(0, 66, r'$\Delta R > 0.2$')
    if i == 3:
        plt.ylim(np.min(etas[:,i][DelR_cut] - .25), np.max(etas[:,i][DelR_cut]) + .25)
        plt.text(0, 42.5, r'Pair ' + str(i+1))
        plt.text(0, 41.5, r'$m_\mathrm{KKKK} < 125$GeV')
        plt.text(0, 40, r'$\Delta R > 0.2$')

#    plt.show()                                                                                                                                                                                            

    plt.savefig('/web/jhornung/public_html/matching/eta_scatter_Pair'+ str(i+1) + '.png')
    plt.savefig('/web/jhornung/public_html/matching/eta_scatter_Pair'+ str(i+1) + '.pdf')

    fig = plt.subplots(figsize=(8,8))

    x = np.linspace(np.min(np.minimum(true_phis[:,i][DelR_cut], phis[:,i][DelR_cut])),np.max(np.maximum(true_phis[:,i][DelR_cut], phis[:,i][DelR_cut])), 50)
    plt.plot(x, x, color='red')

    plt.scatter(true_phis[:,i][DelR_cut], phis[:,i][DelR_cut], marker='x')
    plt.xlabel(r'$\phi^\mathrm{truth}$')
    plt.ylabel(r'$\phi^\mathrm{reco}$')
    if i == 0:
        plt.ylim(np.min(phis[:,i][DelR_cut] - .25), np.max(phis[:,i][DelR_cut]) + .25)
        plt.text(0, 160, r'Pair ' + str(i+1))
        plt.text(0, 155, r'$m_\mathrm{KKKK} < 125$GeV')
        plt.text(0, 150, r'$\Delta R > 0.2$')
    if i == 1:
        plt.ylim(np.min(phis[:,i][DelR_cut] - .25), np.max(phis[:,i][DelR_cut]) + .25)
        plt.text(0, 150, r'Pair ' + str(i+1))
        plt.text(0, 145, r'$m_\mathrm{KKKK} < 125$GeV')
        plt.text(0, 140, r'$\Delta R > 0.2$')
    if i == 2:
        plt.ylim(np.min(phis[:,i][DelR_cut] - .25), np.max(phis[:,i][DelR_cut]) + .25)
        plt.text(0, 70, r'Pair ' + str(i+1))
        plt.text(0, 68, r'$m_\mathrm{KKKK} < 125$GeV')
        plt.text(0, 66, r'$\Delta R > 0.2$')
    if i == 3:
        plt.ylim(np.min(phis[:,i][DelR_cut] - .25), np.max(phis[:,i][DelR_cut]) + .25)
        plt.text(0, 42.5, r'Pair ' + str(i+1))
        plt.text(0, 41.5, r'$m_\mathrm{KKKK} < 125$GeV')
        plt.text(0, 40, r'$\Delta R > 0.2$')

#    plt.show()                                                                                                                                                                                             

    plt.savefig('/web/jhornung/public_html/matching/phi_scatter_Pair'+ str(i+1) + '.png')
    plt.savefig('/web/jhornung/public_html/matching/phi_scatter_Pair'+ str(i+1) + '.pdf')
'''    
    
