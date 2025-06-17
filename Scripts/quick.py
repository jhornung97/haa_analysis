import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
import uproot

hep.style.use(hep.style.CMS)

mcfile = uproot.open("/work/jhornung/Haa/KingMaker/CROWN/build/bin/output_mm.root")
mc = mcfile["ntuple"].arrays()

#datafile = uproot.open("/ceph/jhornung/Data_2018/2018/single_muon_data/mm/single_muon_data.root")
#data = datafile["ntuple"].arrays()

fig, ax = plt.subplots()

data_min = np.min(mc["PV_npvsGood"])
data_max = np.max(mc["PV_npvsGood"])
bins = np.arange(data_min, data_max, 1)

hist, edges, _  = ax.hist(mc["PV_npvsGood"], bins=bins, histtype="step", label="sig MC")
hist, edges, _  = ax.hist(mc["PV_npvsGood"], bins=bins, histtype="step", weights=mc['puweight'], label="sig MC w PU weights")
#ax.hist(data["PV_npvsGood"], bins=bins, histtype="step", density=True, label="data")

print(edges)
ax.set_xlabel("PV_npvsGood")
ax.set_ylabel("Events")

ax.legend()
#plt.savefig("/web/jhornung/public_html/analysis_plots/2018/pileup_nPU.png")
plt.show()