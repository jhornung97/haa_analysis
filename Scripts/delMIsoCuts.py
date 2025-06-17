import uproot
import numpy as np
import awkward as ak
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
import mplhep
import sys
import utility as util
import HistLib as hl

mc = sys.argv[1]
scope = sys.argv[2]

if mc == 'signal':
    signal_loc = "/ceph/jhornung/MC_2018/2018"
    signal = ["HaaKKKK_M125_13TeV_powheg_pythia8-RunIISummer20UL18MiniAODv2-106X"]
    signal_labels = [r'$Z(\ell\ell)H(aa\rightarrow 4K)$']
    mc_signal = hl.data(signal_loc, signal, signal_labels, scope)
    signal_branches = mc_signal.branches_dict

    signal_masks = {}
    signal_weights = {}

    for key, value in signal_branches.items():
        H_m_mask = ((value['H_mass'] > 70) & (value['H_mass'] < 200))
        H_eta_mask = (np.abs(value['H_eta']) < 2.4)
        leading_k_mask = (value['d1_pt'] > 10) & (value['d2_pt'] > 10)
        z_pt_mask = (value['pt_vis'] > 30)
        z_m_mask = ((value['m_vis'] > 75) & (value['m_vis'] < 105))
        ps_mask = (value['ps_1_mass'] < 3) & (value['ps_2_mass'] < 3) & (np.abs(value['ps_1_mass'] - value['ps_2_mass']) < 1)
        iso_mask = ((value['d1_iso'] < 3) & (value['d2_iso'] < 3) & (value['d3_iso'] < 3) & (value['d4_iso'] < 3))
        mask = H_m_mask & H_eta_mask & z_m_mask & ps_mask & z_pt_mask & leading_k_mask #& iso_mask
        signal_masks[key] = mask
    
    if scope == "mm":
        sfs = value['id_wgt_mu_1'][mask]*value['id_wgt_mu_2'][mask]*value['iso_wgt_mu_1'][mask]*value['iso_wgt_mu_2'][mask]*value['trigger_wgt_mu_1'][mask]*value['trigger_wgt_mu_2'][mask]
        signs = np.array(list(map(util.weights, value['genWeight'][mask])))
        pu_mean = np.mean(value['puweight'][mask])
        pu_normalized = value['puweight'][mask]/pu_mean
        weight = 59.8e3*value['evtweight'][mask]*signs*sfs*value['puweight'][mask]
        signal_weights[key] = weight
    elif scope == "ee":
        sfs = value['id_wgt_ele_1'][mask]*value['id_wgt_ele_2'][mask]*value['trigger_wgt_ele_1'][mask]*value['trigger_wgt_ele_2'][mask]
        signs = np.array(list(map(util.weights, value['genWeight'][mask])))
        weight = 59.8e3*value['evtweight'][mask]*signs*sfs*value['puweight'][mask]
        signal_weights[key] = weight
    
    signal_d1_iso = [tree['d1_iso'][signal_masks[key]].to_numpy() for key, tree in signal_branches.items()]
    signal_d2_iso = [tree['d2_iso'][signal_masks[key]].to_numpy() for key, tree in signal_branches.items()]
    signal_d3_iso = [tree['d3_iso'][signal_masks[key]].to_numpy() for key, tree in signal_branches.items()]
    signal_d4_iso = [tree['d4_iso'][signal_masks[key]].to_numpy() for key, tree in signal_branches.items()]
    signal_ps_1_mass = [tree['ps_1_mass'][signal_masks[key]].to_numpy() for key, tree in signal_branches.items()]
    signal_ps_2_mass = [tree['ps_2_mass'][signal_masks[key]].to_numpy() for key, tree in signal_branches.items()]
    signal_weights = list(signal_weights.values())
    signal_labels = list(signal_branches.keys())

    delta_mass = signal_ps_1_mass[0] - signal_ps_2_mass[0]
    delta_mass = np.abs(delta_mass)
    delta_mass = np.concatenate((delta_mass, delta_mass))

    leading_isos = np.concatenate((signal_d1_iso[0], signal_d2_iso[0]))
    subleading_isos = np.concatenate((signal_d3_iso[0], signal_d4_iso[0]))
    iso_weights = np.concatenate((signal_weights[0].to_numpy(), signal_weights[0].to_numpy()))

    norm = np.sum(iso_weights)

    fig, ax = plt.subplots()
    mplhep.style.use("CMS")
    
    h, xedges, y_edges, image = ax.hist2d(delta_mass, leading_isos, range=((0,0.1),(0,0.8)), bins=(10, 20), weights=iso_weights)
    
    ax.set_xlabel(r"$|\Delta m_{K^+K^-}|$", fontsize=20)
    ax.set_ylabel(r"Leading $K^\pm$ Isolation", fontsize=20)
    ax.tick_params(axis='both', which='major', labelsize=20)

    rate = h / norm
    pcm = ax.pcolormesh(xedges, y_edges, rate.T, cmap='viridis')

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    cbar = fig.colorbar(pcm, cax=cax)
    cbar.set_label("Per-Kaon Rate", fontsize=20)
    cbar.ax.tick_params(labelsize=20)

    ax.set_title(r"$\it{" + "Private \: work" + "}$" +" "+ r"$\bf{" + "(CMS \: simulation)" + "}$",loc="left", fontsize=20)
    ax.set_title(f"59.8 fb$^{{-1}}$, 2018 (13 $\mathrm{{{{TeV}}}}$)",loc="right", fontsize=20)

    plt.savefig(f"/web/jhornung/public_html/analysis_plots/2018/Delta_mKK_iso_2D_{scope}_AN_{mc}_leading_iso.png", bbox_inches='tight')
    plt.savefig(f"/web/jhornung/public_html/analysis_plots/2018/Delta_mKK_iso_2D_{scope}_AN_{mc}_leading_iso.pdf", bbox_inches='tight')

    fig, ax = plt.subplots()
    mplhep.style.use("CMS")

    h, xedges, y_edges, image = ax.hist2d(delta_mass, subleading_isos, range=((0,0.1),(0,0.8)), bins=(10, 20), weights=iso_weights)

    ax.set_xlabel(r"$|\Delta m_{K^+K^-}|$", fontsize=20)
    ax.set_ylabel(r"Subleading $K^\pm$ Isolation", fontsize=20)
    ax.tick_params(axis='both', which='major', labelsize=20)

    rate = h / norm
    pcm = ax.pcolormesh(xedges, y_edges, rate.T, cmap='viridis')
    
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    cbar = fig.colorbar(pcm, cax=cax)
    cbar.set_label("Per-Kaon Rate", fontsize=20)
    cbar.ax.tick_params(labelsize=20)

    ax.set_title(r"$\it{" + "Private \: work" + "}$" +" "+ r"$\bf{" + "(CMS \: simulation)" + "}$",loc="left", fontsize=20)
    ax.set_title(f"59.8 fb$^{{-1}}$, 2018 (13 $\mathrm{{{{TeV}}}}$)",loc="right", fontsize=20)

    plt.savefig(f"/web/jhornung/public_html/analysis_plots/2018/Delta_mKK_iso_2D_{scope}_AN_{mc}_subleading_iso.png", bbox_inches='tight')
    plt.savefig(f"/web/jhornung/public_html/analysis_plots/2018/Delta_mKK_iso_2D_{scope}_AN_{mc}_subleading_iso.pdf", bbox_inches='tight')

    plt.show()


elif mc == 'bkg':
    mc_loc = f"/ceph/jhornung/MC_2018/2018"
    bkgs = [
                "DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8+RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2+MINIAODSIM",
                "ttbar",
                "diboson",
                "wh",
                "single_top"
                #"qcd"
                ]
    bkg_labels = [
        "Drell-Yan",
        r"$t\bar{t}$",
        "Diboson",
        r'$W\,H\rightarrow \ell\,\nu\, b\,\bar{b}$',
        "Single Top"
    ]
    mc_bkg = hl.data(mc_loc, bkgs, bkg_labels, scope)
    bkg_branches = mc_bkg.branches_dict
    z_pt_weights = mc_bkg.get_Z_pt_weights() 

    bkg_masks = {}
    bkg_weights = {}

    for key, value in bkg_branches.items():
        H_m_mask = ((value['H_mass'] > 70) & (value['H_mass'] < 200))
        H_eta_mask = (np.abs(value['H_eta']) < 2.4)
        leading_k_mask = (value['d1_pt'] > 10) & (value['d2_pt'] > 10)
        z_pt_mask = (value['pt_vis'] > 30)
        z_m_mask = ((value['m_vis'] > 75) & (value['m_vis'] < 105))
        ps_mask = (value['ps_1_mass'] < 3) & (value['ps_2_mass'] < 3) & (np.abs(value['ps_1_mass'] - value['ps_2_mass']) < 1)
        iso_mask = ((value['d1_iso'] < 3) & (value['d2_iso'] < 3) & (value['d3_iso'] < 3) & (value['d4_iso'] < 3))
        combined_mask = (((np.abs(value['ps_1_mass'] - value['ps_2_mass'])/0.06)**2 + (value['d1_iso']/3.5)**2) < 1)
        combined_mask = combined_mask & (((np.abs(value['ps_1_mass'] - value['ps_2_mass'])/0.06)**2 + (value['d2_iso']/3.5)**2) < 1)
        mask = H_m_mask & H_eta_mask & z_m_mask & ps_mask & z_pt_mask & leading_k_mask #& combined_mask
        bkg_masks[key] = mask
        
        if key == "DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8+RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2+MINIAODSIM":
            weight = 59.8e3*value['evtweight'][mask]*z_pt_weights[mask]*value['puweight'][mask]
        else:
            weight = 59.8e3*value['evtweight'][mask]*value['puweight'][mask]
        bkg_weights[key] = weight

    bkg_d1_iso = [tree['d1_iso'][bkg_masks[key]].to_numpy() for key, tree in bkg_branches.items()]
    bkg_d2_iso = [tree['d2_iso'][bkg_masks[key]].to_numpy() for key, tree in bkg_branches.items()]
    bkg_d3_iso = [tree['d3_iso'][bkg_masks[key]].to_numpy() for key, tree in bkg_branches.items()]
    bkg_d4_iso = [tree['d4_iso'][bkg_masks[key]].to_numpy() for key, tree in bkg_branches.items()]
    bkg_ps_1_mass = [tree['ps_1_mass'][bkg_masks[key]].to_numpy() for key, tree in bkg_branches.items()]
    bkg_ps_2_mass = [tree['ps_2_mass'][bkg_masks[key]].to_numpy() for key, tree in bkg_branches.items()]
    bkg_weights = list(bkg_weights.values())
    bkg_labels = list(bkg_branches.keys())

    bkg_d1_iso = np.concatenate([bkg_d1_iso[i] for i in range(len(bkg_labels))])
    bkg_d2_iso = np.concatenate([bkg_d2_iso[i] for i in range(len(bkg_labels))])
    bkg_d3_iso = np.concatenate([bkg_d3_iso[i] for i in range(len(bkg_labels))])
    bkg_d4_iso = np.concatenate([bkg_d4_iso[i] for i in range(len(bkg_labels))])
    bkg_ps_1_mass = np.concatenate([bkg_ps_1_mass[i] for i in range(len(bkg_labels))])
    bkg_ps_2_mass = np.concatenate([bkg_ps_2_mass[i] for i in range(len(bkg_labels))])
    bkg_weights = np.concatenate([bkg_weights[i].to_numpy() for i in range(len(bkg_labels))])

    delta_mass = bkg_ps_1_mass - bkg_ps_2_mass
    delta_mass = np.abs(delta_mass)
    delta_mass = np.concatenate((delta_mass, delta_mass))

    leading_isos = np.concatenate((bkg_d1_iso, bkg_d2_iso))
    subleading_isos = np.concatenate((bkg_d3_iso, bkg_d4_iso))
    iso_weights = np.concatenate((bkg_weights, bkg_weights))

    norm = np.sum(iso_weights)

    fig, ax = plt.subplots()
    mplhep.style.use("CMS")

    h, xedges, y_edges, image = ax.hist2d(delta_mass, leading_isos, range=((0,.1),(0,6)), bins=(10, 24), weights=iso_weights)
    ax.set_xlabel(r"$|\Delta m_{K^+K^-}|$", fontsize=20)
    ax.set_ylabel(r"Leading $K^\pm$ Isolation", fontsize=20)
    ax.tick_params(axis='both', which='major', labelsize=20)

    #theta = np.linspace(0, np.pi/2, 200)
    #x = 0.06*np.cos(theta)
    x = np.linspace(0, 0.1, 10)
    #y = 3.5*np.sin(theta)
    y = -(66+2./3.)*x+4

    ax.plot(x, y, color='black', linewidth=2, linestyle='--')

    rate = h / norm

    pcm = ax.pcolormesh(xedges, y_edges, rate.T)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    cbar = fig.colorbar(pcm, cax=cax)
    cbar.set_label("Per-Kaon Rate", fontsize=20)
    cbar.ax.tick_params(labelsize=20)

    ax.set_title(r"$\it{" + "Private \: work" + "}$" +" "+ r"$\bf{" + "(CMS \: simulation)" + "}$",loc="left", fontsize=20)
    ax.set_title(f"59.8 fb$^{{-1}}$, 2018 (13 $\mathrm{{{{TeV}}}}$)",loc="right", fontsize=20)

    plt.savefig(f"/web/jhornung/public_html/analysis_plots/2018/Delta_mKK_iso_2D_{scope}_AN_{mc}_leading_iso_line_cut.png", bbox_inches='tight')
    plt.savefig(f"/web/jhornung/public_html/analysis_plots/2018/Delta_mKK_iso_2D_{scope}_AN_{mc}_leading_iso_line_cut.pdf", bbox_inches='tight')

    fig, ax = plt.subplots()
    mplhep.style.use("CMS")
    h, xedges, y_edges, image = ax.hist2d(delta_mass, subleading_isos, range=((0,0.1),(0,6)), bins=(10, 24), weights=iso_weights)
    ax.set_xlabel(r"$|\Delta m_{K^+K^-}|$", fontsize=20)
    ax.set_ylabel(r"Subleading $K^\pm$ Isolation", fontsize=20)
    ax.tick_params(axis='both', which='major', labelsize=20)
    
    plt.plot(x, y, color='black', linewidth=2, linestyle='--')

    rate = h / norm
    pcm = ax.pcolormesh(xedges, y_edges, rate.T)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    cbar = fig.colorbar(pcm, cax=cax)
    cbar.set_label("Per-Kaon Rate", fontsize=20)
    cbar.ax.tick_params(labelsize=20)
    
    ax.set_title(r"$\it{" + "Private \: work" + "}$" +" "+ r"$\bf{" + "(CMS \: simulation)" + "}$",loc="left", fontsize=20)
    ax.set_title(f"59.8 fb$^{{-1}}$, 2018 (13 $\mathrm{{{{TeV}}}}$)",loc="right", fontsize=20)

    plt.savefig(f"/web/jhornung/public_html/analysis_plots/2018/Delta_mKK_iso_2D_{scope}_AN_{mc}_subleading_iso_line_cut.png", bbox_inches='tight')
    plt.savefig(f"/web/jhornung/public_html/analysis_plots/2018/Delta_mKK_iso_2D_{scope}_AN_{mc}_subleading_iso_line_cut.pdf", bbox_inches='tight')

    plt.show()