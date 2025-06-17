import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.patches as mpatches
import ROOT

def exponential(x, N, a):
    return N * np.exp(a * x)

def normalize_exponential(a, x_min, x_max, total_events):
    integral = (np.exp(a * x_max) - np.exp(a * x_min)) / a
    N = total_events / integral
    return 2*N

def normalize_error(a, a_err, x_min, x_max, total_events, total_events_err):
    return np.sqrt((2*a/(np.exp(a*x_max) - np.exp(a*x_min)))**2 * total_events_err**2 + ((2*total_events)/(np.exp(a*x_max) - np.exp(a*x_min))-2*total_events*a*(x_max*np.exp(a*x_max)-x_min*np.exp(a*x_min))/(np.exp(a*x_max)-np.exp(a*x_min))**2)**2 * a_err**2)

fig, ax = plt.subplots(2,4, figsize=(24, 12))

for i, year in enumerate(["2016preVFP", "2016postVFP", "2017", "2018"]):
    for j, scope in enumerate(["mm", "ee"]):

        datafile = ROOT.TFile(f"/work/jhornung/Haa/limits/{year}/{scope}/workspace_bkg_{scope}.root")
        w = datafile.Get(f"decorrelated_workspace_bkg_{scope}")
        data = w.data(f"data_H_mass_{year}_{scope}")
        nEntries = data.numEntries()

        print("Entries: ", nEntries)

        fitResultFile = ROOT.TFile("/work/jhornung/Haa/limits/all_years_combined/fitDiagnosticsTest.root")
        fitResult = fitResultFile.Get("fit_s")
        param_alpha = fitResult.floatParsFinal().find(f"alpha_{year}_{scope}")
        param_norm = fitResult.floatParsFinal().find(f"shapeBkg_bkg_mass_{year}_H_mass_{year}_{scope}__norm")
        alpha, alpha_err = param_alpha.getVal(), param_alpha.getError()
        norm, norm_err = param_norm.getVal(), param_norm.getError()

        print("Norm: ", norm, " +/- ", norm_err)

        data_array = np.array([data.get(n).getRealValue(f"H_mass_{year}_{scope}") for n in range(data.numEntries())])

        mask = (data_array < 110) | (data_array > 140)
        hist, bins = np.histogram(data_array[mask],  bins=65, range=(70, 200))
        integral = np.sum(hist)
        stat_err = np.sqrt(hist[hist > 0])

        bin_centers = (bins[:-1] + bins[1:]) / 2

        plot_range = np.linspace(70, 200, 1000)
        N = normalize_exponential(alpha, plot_range[0], plot_range[-1], norm)
        N_up = N + norm_err
        N_down = N - norm_err
        N_two_up = N + 2*norm_err
        N_two_down = N - 2*norm_err
        
        alpha_up = alpha + alpha_err
        alpha_down = alpha - alpha_err
        alpha_two_up = alpha + 2*alpha_err
        alpha_two_down = alpha - 2*alpha_err
        
        ax[j,i].errorbar(bin_centers[hist > 0], hist[hist > 0], yerr=stat_err, marker='o', linestyle='none', color='black', markersize=2,  label=f'{year} {scope} Data')
        ax[j,i].plot(plot_range, exponential(plot_range, N, alpha), color='red', label=r'$\alpha_{%s}^{%s} = %.4f \pm %.4f$' % (year, scope, alpha, alpha_err))

        #ax[j,i].fill_between(plot_range, exponential(plot_range, N_two_up, alpha_two_up), exponential(plot_range, N_two_down, alpha_two_down), facecolor='gold', edgecolor='gold', alpha=1, label=r'$2\sigma$ Fit uncertainty')
        #ax[j,i].fill_between(plot_range, exponential(plot_range, N_up, alpha_up), exponential(plot_range, N_down, alpha_down), facecolor='limegreen', edgecolor='limegreen', alpha=1, label=r'$1\sigma$ Fit uncertainty')

        ax[j,i].legend()
        ax[j,i].set_xlabel(r'$m_{4K}$ [GeV]')

plt.tight_layout()
plt.savefig("/web/jhornung/public_html/analysis_plots/check_bkg_fits.pdf")
plt.savefig("/web/jhornung/public_html/analysis_plots/check_bkg_fits.png")
plt.show()