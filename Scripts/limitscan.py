import numpy as np
import matplotlib.pyplot as plt
import os
import uproot
import mplhep
mplhep.style.use('CMS')

data = np.loadtxt('/ceph/jhornung/iso_limitscan_2018', delimiter=' ')
sorted_idx = np.argsort(data[:, 0])
data = data[sorted_idx]
iso_cuts, unique = np.unique(data[:, 0], return_index=True)
limits_18 = data[unique, 1]

data = np.loadtxt('/ceph/jhornung/iso_limitscan_2017', delimiter=' ')
sorted_idx = np.argsort(data[:, 0])
data = data[sorted_idx]
iso_cuts_17, unique = np.unique(data[:, 0], return_index=True)
limits_17 = data[unique, 1]

data = np.loadtxt('/ceph/jhornung/iso_limitscan_2016postVFP', delimiter=' ')
sorted_idx = np.argsort(data[:, 0])
data = data[sorted_idx]
iso_cuts_16post, unique = np.unique(data[:, 0], return_index=True)
limits_16post = data[unique, 1]

data = np.loadtxt('/ceph/jhornung/iso_limitscan_2016preVFP', delimiter=' ')
sorted_idx = np.argsort(data[:, 0])
data = data[sorted_idx]
iso_cuts_16pre, unique = np.unique(data[:, 0], return_index=True)
limits_16pre = data[unique, 1]

fig, axs = plt.subplots()

plt.plot(iso_cuts_16pre[(limits_16pre > 0)], limits_16pre[(limits_16pre > 0)]/.88/.101, '-', label='2016preVFP')
plt.plot(iso_cuts_16post[(limits_16post > 0)], limits_16post[(limits_16post > 0)]/.88/.101, '-', label='2016postVFP')
plt.plot(iso_cuts_17[(limits_17 > 0)], limits_17[(limits_17 > 0)]/.88/.101, '-', label='2017')
plt.plot(iso_cuts[(limits_18 > 0)], limits_18[(limits_18 > 0)]/.88/.101, '-', label='2018')

plt.legend(fontsize=20)

axs.set_xlabel('Iso cut', fontsize=20)
axs.set_ylabel(r"$BR(H\rightarrow aa\rightarrow K^+K^-K^+K^-)$", fontsize=20)

axs.set_title(r"$\it{" + "Private \: work" + "}$" +" "+ r"$\bf{" + "(CMS \: data/simulation)" + "}$",loc="left", fontsize=20)
axs.set_title(f"{101.3+36.3} fb$^{{-1}}$ (13 $\mathrm{{{{TeV}}}}$)",loc="right", fontsize=20)

print(iso_cuts_17[(iso_cuts_17 < 5) & (limits_17 > 0.02)])

print(iso_cuts_16pre[limits_16pre < 0])
print(iso_cuts_16post[limits_16post < 0])
print(iso_cuts_17[limits_17 < 0])
print(iso_cuts[limits_18 < 0])

plt.savefig('/web/jhornung/public_html/analysis_plots/limitscan_iso_cut.pdf', bbox_inches='tight')
plt.savefig('/web/jhornung/public_html/analysis_plots/limitscan_iso_cut.png', bbox_inches='tight')

plt.show()

