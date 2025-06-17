import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
import uproot
from matplotlib.patches import Circle

def avg_iso(pfcands_iso):
    pfcands_iso = np.array(pfcands_iso) 
    iso_weight = np.where(pfcands_iso > 1, 1/pfcands_iso.shape[0], 1)
    return np.average(pfcands_iso, weights=iso_weight)

hep.style.use(hep.style.CMS)

ifile = uproot.open("/work/jhornung/Haa/KingMaker/CROWN/build/bin/output_mm.root")
tree = ifile["ntuple"]
branches = tree.arrays()

pfcands_pt = branches["pfcands_pt"][1000].to_numpy()
pfcands_eta = branches["pfcands_eta"][1000].to_numpy()
pfcands_phi = branches["pfcands_phi"][1000].to_numpy()
pfcands_mass = branches["pfcands_mass"][1000].to_numpy()
pfcands_iso = branches["pfcands_iso"][1000].to_numpy()
pfcands_pdgid = branches["pfcands_pdgid"][1000].to_numpy()

genpart_pdgid = branches["genpart_pdgid"][1000].to_numpy()
genpart_status = branches["genpart_status"][1000].to_numpy()
genpart_eta = branches["genpart_eta"][1000].to_numpy()
genpart_phi = branches["genpart_phi"][1000].to_numpy()

mask_iso = (pfcands_iso < 1) 
mask_pt = (pfcands_pt > 1)
mask_eta = (np.abs(pfcands_eta) < 2.5) 
mask_mass = (pfcands_mass > 0.139) & (pfcands_mass < 0.14)
mask_pdgid = (np.abs(pfcands_pdgid) == 211)

genpart_pdgid_mask = (np.abs(genpart_pdgid) == 9000006)
#genpart_status_mask = (genpart_status == 62)

#print("gen part pdgid: ", genpart_pdgid[genpart_pdgid_mask])

fig, ax = plt.subplots()

norm = plt.Normalize(np.min(pfcands_iso[mask_eta]), np.max(pfcands_iso[mask_eta]))
cmap = plt.get_cmap("viridis")

sc = ax.scatter(pfcands_eta[mask_eta],
                pfcands_phi[mask_eta], 
                c=pfcands_iso[mask_eta],
                cmap=cmap,
                norm=norm,
                s=10)

ax.scatter(genpart_eta[genpart_pdgid_mask],
           genpart_phi[genpart_pdgid_mask],
           c="red",
           s=10,
           label="Gen lvl $a$ bosons")

cbar = fig.colorbar(sc, ax=ax)
cbar.set_label("PFCand Isolation")

for eta, phi in zip(genpart_eta[genpart_pdgid_mask], genpart_phi[genpart_pdgid_mask]):
    ax.add_patch(Circle((eta, phi), 0.25, color="red", fill=False))
    deta = eta - pfcands_eta[mask_eta]
    dphi = np.abs(phi - pfcands_phi[mask_eta])
    dphi = np.where(dphi > np.pi, dphi - 2*np.pi, dphi)
    dr = np.sqrt(deta**2 + dphi**2)
    indices_sorted_by_pt = np.argsort(pfcands_pt[mask_eta][(dr < 0.25)])[::-1]
    sorted_pts = pfcands_pt[mask_eta][(dr < 0.25)][indices_sorted_by_pt]
    print("Two leading PFCand pt inside cone: ", sorted_pts[:2])
    print("Their isos: ", pfcands_iso[mask_eta][(dr < 0.25)][indices_sorted_by_pt][:2])
    print("Their pdgids: ", pfcands_pdgid[mask_eta][(dr < 0.25)][indices_sorted_by_pt][:2])
    print("PT values larger than the leading a daughter PFCand: ", pfcands_pt[mask_eta][pfcands_pt[mask_eta] > sorted_pts[0]])
    print("Their pdgids: ", pfcands_pdgid[mask_eta][pfcands_pt[mask_eta] > sorted_pts[0]])

#for pt, eta, phi in zip(pfcands_pt[mask_iso & mask_eta], pfcands_eta[mask_iso & mask_eta], pfcands_phi[mask_iso & mask_eta]):
#    ax.add_patch(Circle((eta, phi), 0.25, color="blue", fill=False))
#    ax.add_patch(Circle((eta, phi), 0.05, color="blue", fill=False))
#    print(f"PFCand pt: {pt}, PFCand eta: {eta}, PFCand phi: {phi}")
#    deta = eta - pfcands_eta[mask_eta]
#    dphi = np.abs(phi - pfcands_phi[mask_eta])
#    dphi = np.where(dphi > np.pi, dphi - 2*np.pi, dphi)
#    dr = np.sqrt(deta**2 + dphi**2)
#    sorted_pts = np.sort(pfcands_pt[mask_eta][(dr > 0.05) & (dr < 0.25)])[::-1]
#    print("Two leading PFCand pt inside cone: ", sorted_pts[:2])

hep.cms.label("Private Work", data=False, year="2018", lumi=59.83, loc=1)

#plt.legend()

ax.set_xlabel(r"$\eta$")
ax.set_ylabel(r"$\phi$")

plt.tight_layout()
#plt.savefig("/web/jhornung/public_html/analysis_plots/2018/signal_sim_pf_cands_a_highlighted_eta_phi_plot.png")
#plt.savefig("/web/jhornung/public_html/analysis_plots/2018/signal_sim_pf_cands_a_highlighted_eta_phi_plot.pdf")
plt.show()

pfcands_iso = branches["pfcands_iso"]
pfcands_eta = branches["pfcands_eta"]

mask_eta = (np.abs(pfcands_eta) < 2.5)

average_iso = np.array([avg_iso(pfcands_iso[i][mask_eta[i]]) for i in range(len(pfcands_iso))])

fig, ax = plt.subplots()
counts, bins, _ = plt.hist(average_iso, bins=50, histtype="step", color="black")
hep.cms.label("Private Work", data=False, year="2018", lumi=59.83, loc=3)
ax.set_xlabel("Average PFCand Isolation")
ax.set_ylabel("Events/Iso")
ax.set_ylim(0, np.max(counts)*1.15)
plt.tight_layout()
#plt.savefig("/web/jhornung/public_html/analysis_plots/2018/signal_sim_pf_cands_avg_iso_plot.png")
#plt.savefig("/web/jhornung/public_html/analysis_plots/2018/signal_sim_pf_cands_avg_iso_plot.pdf")
#plt.show()

print("Mean: ", np.mean(average_iso))

genpart_pdgid = branches["genpart_pdgid"]
genpart_status = branches["genpart_status"]
genpart_eta = branches["genpart_eta"]
genpart_phi = branches["genpart_phi"]

genpart_pdgid_mask = (np.abs(genpart_pdgid) == 9000006)
genpart_eta_mask = (np.abs(genpart_eta) < 2.5)

sample_leading_isos = []

for event in range(len(genpart_pdgid)):
    if event > 1000:
        break
    a_boson_eta = genpart_eta[event][genpart_pdgid_mask[event] & genpart_eta_mask[event]].to_numpy()
    a_boson_phi = genpart_phi[event][genpart_pdgid_mask[event] & genpart_eta_mask[event]].to_numpy()

    if a_boson_eta.shape[0] != 2: continue

    a_boson_eta = a_boson_eta[:, np.newaxis]
    a_boson_phi = a_boson_phi[:, np.newaxis]

    pfcands_pt = branches["pfcands_pt"][event][mask_eta[event]].to_numpy()
    pfcands_eta = branches["pfcands_eta"][event][mask_eta[event]].to_numpy()
    pfcands_phi = branches["pfcands_phi"][event][mask_eta[event]].to_numpy()
    pfcands_iso = branches["pfcands_iso"][event][mask_eta[event]].to_numpy()

    deta = a_boson_eta - pfcands_eta
    dphi = np.abs(a_boson_phi - pfcands_phi)
    dphi = np.where(dphi > np.pi, dphi - 2*np.pi, dphi)
    dr = np.sqrt(deta**2 + dphi**2)

    dr_mask = dr < 0.25

    pf_cands_pt = np.array([pfcands_pt, pfcands_pt])
    pf_cands_iso = np.array([pfcands_iso, pfcands_iso])

    pt_sorted_indices = np.argsort(pf_cands_pt[dr_mask])[::-1]
    leading_isos = pf_cands_iso[dr_mask][pt_sorted_indices][:4]

    if np.sum(leading_isos > 1) > 0:
        print("Event: ", event)
        print("Leading pt: ", pf_cands_pt[dr_mask][pt_sorted_indices][:4])
        print("Leading isos: ", leading_isos)

    sample_leading_isos.append(leading_isos)

sample_leading_isos = np.array(sample_leading_isos)

fig, ax = plt.subplots()
counts, bins, _ = plt.hist(sample_leading_isos.flatten(), bins=50, histtype="step", color="black")

hep.cms.label("Private Work", data=False, year="2018", lumi=59.83, loc=3)

ax.set_xlabel("Leading PFCand Isolation")
ax.set_ylabel("Events/Iso")

ax.set_ylim(0, np.max(counts)*1.15)

plt.tight_layout()
#plt.savefig("/web/jhornung/public_html/analysis_plots/2018/signal_sim_pf_cands_leading_iso_plot.png")
#plt.savefig("/web/jhornung/public_html/analysis_plots/2018/signal_sim_pf_cands_leading_iso_plot.pdf")

plt.show()
