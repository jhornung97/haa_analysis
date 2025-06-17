import uproot
import matplotlib.pyplot as plt
import numpy as np
import mplhep as hep

plt.style.use(hep.style.CMS)

ifile_KIT = uproot.open("/work/jhornung/Haa/KingMaker/CROWN/data/pileup/Data_Pileup_2018_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18.root")
pileup_hist_KIT = ifile_KIT["pileup"]

ifile_mit = uproot.open("/work/jhornung/Haa/MitAnalysisRunII/data/90x/pu/puWeights_90x_2018.root")
pileup_hist_mit = ifile_mit["puWeights"]

pileup_edges_KIT = pileup_hist_KIT.axis().edges()
pileup_values_KIT = pileup_hist_KIT.values()

print(pileup_edges_KIT)

pileup_centers_KIT = (pileup_edges_KIT[1:] + pileup_edges_KIT[:-1]) / 2

pileup_edges_mit = pileup_hist_mit.axis().edges()
pileup_values_mit = pileup_hist_mit.values()

pileup_centers_mit = (pileup_edges_mit[1:] + pileup_edges_mit[:-1]) / 2

fig, ax = plt.subplots()
#plt.scatter(pileup_centers_KIT[249], pileup_values_KIT[249], label='Our pick')
#plt.scatter(pileup_centers_KIT[250], pileup_values_KIT[250], label='Guillelmo\'s pick')

plt.step(pileup_centers_KIT, pileup_values_KIT, where='mid', label='KIT')
#plt.step(pileup_centers_mit, pileup_values_mit, where='mid', label='MIT')
plt.xlabel('Number of interactions')
plt.ylabel('weight')
plt.legend()
ax.set_xlim(24,26)
ax.set_ylim(bottom=0)
#ticks = np.arange(0, 21, 1)
#ax.set_xticks(ticks)
#ax.tick_params(axis='x', which='minor', bottom=False)
plt.savefig("/web/jhornung/public_html/analysis_plots/2018/pileup_weights.png")

plt.show()