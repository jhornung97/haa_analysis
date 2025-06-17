import uproot
import numpy as np
import matplotlib.pyplot as plt
import utility as util

# Open the files
nlo = uproot.open("/ceph/jhornung/reweighting/2018/nlo.root")  # hier sind die gen_z_pts (pgdId == 23 & status == 62) + genWeight von der private production, die ich auch so verwende, drin
nnlo = uproot.open("/ceph/jhornung/reweighting/2017/nnlo.root")# hier die von der offiziellen production, M-10-50 und M-50 hadded

print("files opened")

# Get the trees
nlo_tree = nlo["ntuple"]
nnlo_tree = nnlo["ntuple"]

nlo_branches = nlo_tree.arrays()
nnlo_branches = nnlo_tree.arrays()

print("branches loaded")

nlo_z_pt = nlo_branches["z_pt"].to_numpy()
nnlo_z_pt = nnlo_branches["z_pt"].to_numpy()

print("variables loaded")

nlo_signs = np.array(list(map(lambda x: 1 if x > 0 else -1, nlo_branches['genWeight'].to_numpy())))    # signs of genWeights needed for shape
nnlo_signs = np.array(list(map(lambda x: 1 if x > 0 else -1, nnlo_branches['genWeight'].to_numpy()))) 

print("weights calculated")

## Plot the Z pT distributions
#
#edges = np.linspace(0, 100,101)
#fig, ax = plt.subplots()
#ax.hist(nlo_z_pt, bins=edges, histtype='step', label='NLO', weights=nlo_signs, density=True)
#ax.hist(nnlo_z_pt, bins=edges, histtype='step', label='NNLO', weights=nnlo_signs, density=True)
#ax.legend()
#ax.set_xlim(0, 100)
#ax.set_xlabel("Z pT [GeV]")
#ax.set_ylabel("Normalized counts")
##plt.savefig("/web/jhornung/consistent_k_factors/only_M-50_z_pt_distrs.png")
#
## Calculate the weights
#
#counts_nlo, _ = np.histogram(nlo_z_pt, bins=edges, weights=nlo_signs, density=True)
#counts_nnlo, _ = np.histogram(nnlo_z_pt, bins=edges, weights=nnlo_signs, density=True)
#
##print(counts_nlo)
##print(counts_nnlo)
#
#z_pt_weights = counts_nnlo / counts_nlo
#
#mean = np.mean(z_pt_weights)
##print(mean)
#normalized_weights = z_pt_weights / mean
##print(normalized_weights)
#
## Plot the weights
#
#fig, ax = plt.subplots()
#bin_centers = (edges[1:] + edges[:-1]) / 2
#
#ax.plot(bin_centers, z_pt_weights, linestyle=None, marker='x')
#ax.plot(bin_centers, np.ones_like(bin_centers), linestyle='--', color='red')
#ax.set_xlabel("Z pT [GeV]")
#ax.set_ylabel("Normalized weights")
#ax.set_xlim(0, 100)
##plt.savefig("/web/jhornung/consistent_k_factors/full_z_pt_weights.png")
#plt.show()

edges_0_100 = np.linspace(0, 100, 21)
edges_100_200 = np.linspace(100, 200, 11)
edges_200_1000 = np.linspace(200, 1000, 6)
edges_1000_1300 = np.linspace(1000, 1300, 3)

edges_list = [edges_0_100, edges_100_200, edges_200_1000, edges_1000_1300]
n_edges = [0, 20, 30, 35, 37]

edges = np.concatenate((edges_0_100, edges_100_200, edges_200_1000, edges_1000_1300))
edges, indices = np.unique(edges, return_index=True)
edges = edges[np.argsort(indices)]

counts_nlo, _ = np.histogram(nlo_z_pt, bins=edges, weights=nlo_signs, density=True)
counts_nnlo, _ = np.histogram(nnlo_z_pt, bins=edges, weights=nnlo_signs, density=True)

z_pt_weights = counts_nnlo / counts_nlo

#mean = np.mean(z_pt_weights)

#normalized_weights = z_pt_weights / mean

fig, ax = plt.subplots()

ax.hist(nlo_z_pt, bins=edges, histtype='step', label='NLO', weights=nlo_signs, density=True)
ax.hist(nnlo_z_pt, bins=edges, histtype='step', label='NNLO', weights=nnlo_signs, density=True)
ax.legend()
ax.set_xlim(edges[0], edges[-1])
ax.set_yscale('log')
ax.set_xlabel("Z pT [GeV]")
ax.set_ylabel("Normalized counts")
plt.savefig("/web/jhornung/public_html/consistent_k_factors/z_pt_distrs_full_range.png")

fig, ax = plt.subplots()
bin_centers = (edges[1:] + edges[:-1]) / 2

ax.plot(bin_centers, z_pt_weights, linestyle=None, marker='x')
ax.plot(bin_centers, np.ones_like(bin_centers), linestyle='--', color='red')
ax.set_xlabel("Z pT [GeV]")
ax.set_ylabel("Weights")
ax.set_xlim(edges[0], edges[-1])
plt.savefig("/web/jhornung/public_html/consistent_k_factors/z_pt_weights_full_range.png")
plt.show()

#for i in range(len(n_edges) - 1):
#    fig, ax = plt.subplots()
#    ax.hist(nlo_z_pt, bins=edges, histtype='step', label='NLO', weights=nlo_signs, density=True)
#    ax.hist(nnlo_z_pt, bins=edges, histtype='step', label='NNLO', weights=nnlo_signs, density=True)
#    ax.legend()
#
#    hist_max = max(counts_nlo[n_edges[i]:n_edges[i+1]+1].max(), counts_nnlo[n_edges[i]:n_edges[i+1]+1].max())
#    ax.set_ylim(0, hist_max * 1.1)
#
#    ax.set_xlim(edges[n_edges[i]:n_edges[i+1]+1][0], edges[n_edges[i]:n_edges[i+1]+1][-1])
#
#    plt.savefig(f"/web/jhornung/public_html/consistent_k_factors/z_pt_distrs_{i}.png")
#
#    fig, ax = plt.subplots()
#    bin_centers = (edges[n_edges[i]:n_edges[i+1]+1][1:] + edges[n_edges[i]:n_edges[i+1]+1][:-1]) / 2
#
#    ax.plot(bin_centers, z_pt_weights[n_edges[i]:n_edges[i+1]], marker='x')
#    ax.set_xlim(edges[n_edges[i]:n_edges[i+1]+1][0], edges[n_edges[i]:n_edges[i+1]+1][-1])
#
#    plt.savefig(f"/web/jhornung/public_html/consistent_k_factors/z_pt_weights_{i}.png")
#
#
#
#