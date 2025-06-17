import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
import uproot
import awkward as ak
import sys

hep.style.use(hep.style.CMS)

ifile = uproot.open("/work/jhornung/Haa/limits/all_years_combined/higgsCombine.alpha_2018_ee.MultiDimFit.mH125.root")
tree = ifile['limit']
branches = tree.arrays()

r = ak.to_numpy(branches["alpha_2018_ee"])
deltaNLL = ak.to_numpy(branches["deltaNLL"])
quantileExpected = ak.to_numpy(branches["quantileExpected"])
print(np.min(deltaNLL[quantileExpected != -1]))
plt.figure()

plt.plot(r[2*deltaNLL < 5][1:], 2*deltaNLL[2*deltaNLL < 5][1:], "-", color="black", label="Profile Likelihood")
#plt.plot(r[quantileExpected == -1], 2*deltaNLL[quantileExpected == -1], "x", color="red", label="Best Fit Value")
#plt.legend()
plt.axhline(1, linestyle="--", color="red")
plt.axhline(4, linestyle="--", color="red")
plt.xlabel("r")
plt.ylabel(r"$-2 \Delta \mathrm{ln}(L)$")
#plt.ylim(bottom=2*np.min(deltaNLL[2*deltaNLL < 5][1:]))
hep.cms.label("Private Work", data=True, lumi=41.48+59.83, loc=0)

#plt.savefig("/web/jhornung/public_html/analysis_plots/all_years/profile_likelihood_scan.pdf")
#plt.savefig("/web/jhornung/public_html/analysis_plots/all_years/profile_likelihood_scan.png")

plt.show()