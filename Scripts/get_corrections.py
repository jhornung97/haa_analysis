import ROOT
import pickle, lz4.frame
import uproot
import numpy as np
import sys
from tqdm import tqdm
from array import array
import awkward as ak
import matplotlib.pyplot as plt
import utility
'''
ROOT.gInterpreter.Declare("""
    int get_index(
        ROOT::VecOps::RVec<Int_t> *GenPart_pdgId,                                                       
        ROOT::VecOps::RVec<Int_t> *GenPart_status                                                       
        ){
        for (int i = 0; i < GenPart_pdgId->size(); i++){
                if (GenPart_pdgId->at(i) == 23 && GenPart_status->at(i)==62) return i;
        }
        int j = -99;
        return j;
    }
""")
'''
res = pickle.load(lz4.frame.open("/work/jhornung/Haa/scetlib_dyturboCorrZ.pkl.lz4"))
h = res["Z"]["scetlib_dyturbo_hist"][{"vars": "pdf0"}]

bin_edges = h.axes[2].edges

theory_zpt_hist = h.project("qT").values()

integral = np.sum(theory_zpt_hist)
theory_zpt_dens = theory_zpt_hist/integral
'''
rdf = ROOT.RDataFrame("Events", "/work/jhornung/Haa/FE16FF04-741B-7E4B-8C4B-A81C14017B4A.root")
totalEvents = rdf.Count().GetValue()
print(totalEvents)
rdf = rdf.Define("indexZ", "(int) get_index(&GenPart_pdgId, &GenPart_status)")
rdf = rdf.Filter("indexZ > 0")
rdf = rdf.Define("genbosonpt", "GenPart_pt[indexZ]")
rdf = rdf.Define("genbosonrapidity", "GenPart_eta[indexZ]")
rdf.Snapshot("ntuple", "/work/jhornung/Haa/output.root", ["genbosonpt", "genbosonrapidity", "genWeight"])
'''
idy_reweighting = uproot.open("/ceph/jhornung/reweighting/reweighting.root")
#idy_reweighting = uproot.open("/work/jhornung/Haa/output.root")
idy_reweighting_ntuples = idy_reweighting['ntuple']
idy_reweighting_branches = idy_reweighting_ntuples.arrays()

weights = np.array(list(map(utility.weights,idy_reweighting_branches['genWeight'])))#*6473./idy_reweighting_ntuples.num_entries

idy_gen_pt_dens, idy_gen_pt_edges = np.histogram(idy_reweighting_branches["z_pt"].to_numpy(), bins=bin_edges, weights=weights, density=True)

zpt_weights = theory_zpt_dens/idy_gen_pt_dens

z_pt_mean = np.mean(zpt_weights)

print(z_pt_mean)

normalized_weights = zpt_weights/z_pt_mean

print(np.mean(normalized_weights))

ifile = uproot.open('/work/jhornung/Haa/merged_kfactors_zjets.root')
other_weights_histo = ifile['kfactor_monojet_ewk']
other_weights_edges = other_weights_histo.axes[0].edges()
other_weights_centers = (other_weights_edges[:-1] + other_weights_edges[1:])/2
other_weights = other_weights_histo.values()

fig, ax = plt.subplots()
plt.plot(h.axes[2].centers, theory_zpt_dens, "-", drawstyle="steps-mid", label="Theory")
plt.hist(idy_reweighting_branches["z_pt"].to_numpy(), bins=bin_edges, weights=weights, density=True, label="Gen lvl MC")
plt.legend()
ax.set_xlabel(r"$p_\mathrm{T}$")

plt.savefig('/web/jhornung/public_html/redo_data_to_mc/reweighting_histos.png')
plt.savefig('/web/jhornung/public_html/redo_data_to_mc/reweighting_histos.pdf')

fig, ax = plt.subplots()
plt.plot(h.axes[2].centers, zpt_weights, "x", label=r'Calculated weights')
plt.plot(other_weights_centers, other_weights, "x", label="Weights from k-factors file")
ax.axvspan(100, 150, facecolor='grey', edgecolor='white', alpha=0.2, label="Area with no weights")
ax.set_xlabel(r"$p_\mathrm{T}$")
ax.set_ylabel("weights")
plt.legend()

plt.savefig('/web/jhornung/public_html/redo_data_to_mc/reweighting_weights.png')
plt.savefig('/web/jhornung/public_html/redo_data_to_mc/reweighting_weights.pdf')

fig, ax = plt.subplots()
plt.plot(h.axes[2].centers, normalized_weights, "x", label=r'Calculated weights normalized to mean')
plt.plot(other_weights_centers, other_weights, "x", label="Weights from k-factors file")
limits = ax.get_ylim()
ax.axvspan(100, 150, facecolor='grey', edgecolor='white', alpha=0.2, label="Area with no weights")
ax.set_ylabel("weights")
plt.legend()
plt.savefig('/web/jhornung/public_html/redo_data_to_mc/reweighting_weights_normalized.png')
plt.savefig('/web/jhornung/public_html/redo_data_to_mc/reweighting_weights_normalized.pdf')

plt.show()



