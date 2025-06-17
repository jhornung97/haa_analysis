import uproot
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.patches as mpatches
from matplotlib.legend_handler import HandlerPatch
import matplotlib.colors as mcolors
import mplhep as hep
import utility as util
import sys
import pickle, lz4.frame
from itertools import combinations

plt.style.use(hep.style.CMS)
plt.rcParams.update({'font.size': 20})

class HandlerPatchWithMargin(HandlerPatch):
    def create_artists(self, legend, orig_handle, xdescent, ydescent, width, height, fontsize, trans):
        patch = super().create_artists(legend, orig_handle, xdescent, ydescent, width, height, fontsize, trans)
        for p in patch:
            p.set_linewidth(1)
            p.set_edgecolor('black')
        return patch

class data:
    def __init__(self, loc, datasets, labels, scope):
        self.loc = loc
        self.datasets = datasets
        self.labels = labels
        self.scope = scope
        self.branches_dict = self.read_in()
        
    def read_in(self):
        branches_dict = {}

        for dataset, label in zip(self.datasets, self.labels):
            file_path = f'{self.loc}/{dataset}/{self.scope}/{dataset}.root'
            print(f"Reading in {file_path}")
            file = uproot.open(f'{self.loc}/{dataset}/{self.scope}/{dataset}.root')
            branches = file['ntuple'].arrays()
            branches_dict[label] = branches
            print('Done')
        
        return branches_dict
    
    def get_key(self, dict, possible_keys):
        for key in possible_keys:
            if key in dict:
                return key
        raise KeyError(f"None of the possible keys {possible_keys} found in the dictionary")
    
    def get_Z_pt_weights(self):
        '''
        res = pickle.load(lz4.frame.open("/work/jhornung/Haa/scetlib_dyturboCorrZ.pkl.lz4"))
        h = res["Z"]["scetlib_dyturbo_hist"][{"vars": "pdf0"}]

        theory_edges = h.axes[2].edges

        theory_zpt_hist = h.project("qT").values()
        integral = np.sum(theory_zpt_hist)
        theory_dens = theory_zpt_hist/integral

        idy_reweighting = uproot.open("/ceph/jhornung/reweighting/reweighting.root")
        idy_reweighting_ntuples = idy_reweighting['ntuple']
        idy_reweighting_branches = idy_reweighting_ntuples.arrays()

        idy_gen_pt_hist, idy_gen_pt_edges = np.histogram(idy_reweighting_branches["z_pt"].to_numpy(), 
                                                         bins=theory_edges, 
                                                         weights=np.array(list(map(util.weights,idy_reweighting_branches['genWeight']))),
                                                         density=True
        )

        theory_weights = theory_dens/idy_gen_pt_hist

        mean = np.mean(theory_weights)

        normalized_theory_weights = theory_weights/mean
        '''

        possible_keys = ["Drell-Yan", "MC"]
        key = self.get_key(self.branches_dict, possible_keys)
        
        idy_zpt = self.branches_dict[key]["gen_pt_vis"].to_numpy()
        idy_idx = np.arange(idy_zpt.shape[0])
        idy_zpt_weights = np.ones(idy_zpt.shape[0])
        '''
        for i in range(len(theory_edges)-1):
            mask = (idy_zpt >= theory_edges[i]) & (idy_zpt < theory_edges[i+1])
            tmp_idx = idy_idx[mask]
            idy_zpt_weights[tmp_idx] = normalized_theory_weights[i]
        
        kfactorsfile = uproot.open("/work/jhornung/Haa/merged_kfactors_zjets.root")
        kfactors_histo = kfactorsfile["kfactor_monojet_ewk"]
        kfactors = kfactors_histo.values()
        kfactors_edges = kfactors_histo.axes[0].edges()

        for i in range(len(kfactors_edges)-1):
            mask = (idy_zpt >= kfactors_edges[i]) & (idy_zpt < kfactors_edges[i+1])
            tmp_idx = idy_idx[mask]
            idy_zpt_weights[tmp_idx] = kfactors[i]*idy_zpt_weights[tmp_idx]
        '''
        eejfile = uproot.open("/work/jhornung/Haa/lindert_qcd_nnlo_sf.root")
        eej_histo = eejfile["eej"]
        eej_weights = eej_histo.values()
        eej_weights_axes = eej_histo.axes
        eej_edges = eej_weights_axes[0].edges() 

        print("Assigning Z pt weights")

        for i in range(len(eej_edges)-1):
            mask = (idy_zpt >= eej_edges[i]) & (idy_zpt < eej_edges[i+1])
            tmp_idx = idy_idx[mask]
            idy_zpt_weights[tmp_idx] = eej_weights[i]*idy_zpt_weights[tmp_idx]
        
        print("Done")

        return idy_zpt_weights
    
class hist:
    def __init__(self, config):
        self.scope = config['scope']
        self.era = config['era']
        if self.era == "2018":
            self.lumi = 59.83e3
        elif self.era == "2017":
            self.lumi = 41.48e3
        elif self.era == "2016preVFP":
            self.lumi = 19.5e3
        elif self.era == "2016postVFP":
            self.lumi = 16.8e3
        self.mode = config['mode']
        self.selection = config['selection']
        self.beeeater = self.beeeater_cmap()
        if self.mode == "data to mc":
            self.mc_branches = config['mc']
            self.signal_branches = config['signal']
            self.data_branches = config['data']
            self.idy_zpt_weights = config['z_pt_weights']
            self.var = config['var']
            #self.test_mask = config['test_mask']
            if self.signal_branches == None and len(self.mc_branches) > 0:
                self.mc_masks, self.mc_weights, self.mc_weights_up, self.mc_weights_down = self.prep_mc_data()
                self.data_masks = self.prep_data()
            else:
                self.mc_masks, self.mc_weights, self.mc_weights_up, self.mc_weights_down  = self.prep_mc_data()
                self.signal_masks, self.signal_weights = self.prep_signal_data()
                self.data_masks = self.prep_data()
            self.diff = False
            self.append = False
        elif self.mode == "diff":
            self.mc_branches = config['mc']
            self.signal_branches = config['signal']
            self.data_branches = config['data']
            self.idy_zpt_weights = config['z_pt_weights']
            self.var = config['var_1']
            self.var_1 = config['var_1']
            self.var_2 = config['var_2']
            self.mc_masks, self.mc_weights, self.mc_weights_up, self.mc_weights_down = self.prep_mc_data()
            self.signal_masks, self.signal_weights = self.prep_signal_data()
            self.data_masks = self.prep_data()
            self.diff = True
            self.append = False
        elif self.mode == "append":
            self.mc_branches = config['mc']
            self.signal_branches = config['signal']
            self.data_branches = config['data']
            self.idy_zpt_weights = config['z_pt_weights']
            self.var = config['var_1']
            self.var_1 = config['var_1']
            self.var_2 = config['var_2']
            self.mc_masks, self.mc_weights, self.mc_weights_up, self.mc_weights_down = self.prep_mc_data()
            self.signal_masks, self.signal_weights = self.prep_signal_data()
            self.data_masks = self.prep_data()
            self.diff = False
            self.append = True
        elif self.mode == "reco to truth":
            self.signal_branches = config['signal']
            self.var = config['var']
            self.signal_masks, self.signal_weights = self.prep_signal_data()
        elif self.mode == "reco to truth 2D":
            self.signal_branches = config['signal']
            self.var = config['var']
            self.var2 = config['var2']
            self.signal_masks, self.signal_weights = self.prep_signal_data()
        elif self.mode == "mc":
            self.mc_branches = config['mc']
            self.signal_branches = config['signal']
            self.var = config['var']
            self.idy_zpt_weights = config['z_pt_weights']
            self.mc_ps_mass_masks = self.prep_ps_mass_mask(self.mc_branches)
            self.signal_ps_mass_masks = self.prep_ps_mass_mask(self.signal_branches)
            self.mc_masks, self.mc_weights = self.prep_mc_data()
            self.signal_masks, self.signal_weights = self.prep_signal_data()
        elif self.mode == "no signal":
            self.mc_branches = config['mc']
            self.data_branches = config['data']
            self.idy_zpt_weights = config['z_pt_weights']
            self.mc_ps_mass_masks = self.prep_ps_mass_mask(self.mc_branches)
            self.data_ps_mass_masks = self.prep_ps_mass_mask(self.data_branches)
            self.mc_masks, self.mc_weights = self.prep_mc_data()
            self.data_masks = self.prep_data()
        elif self.mode == "check sfs":
            self.signal_branches = config['signal']
            self.var = config['var']
            self.signal_masks, self.signal_weights = self.prep_signal_data()
            self.diff = False
            self.append = False
        self.heatmap = config['heatmap']
        if self.heatmap == True or self.mode == "reco to truth 2D":
            self.xlower = config['xlower']
            self.xupper = config['xupper']
            self.ylower = config['ylower']
            self.yupper = config['yupper']
        else:
            self.lower = config['lower']
            self.upper = config['upper']
        self.bins = config['bins']
        self.region = {"mm": r"$\mu\mu$ SR", "ee": r"ee SR", "em": r"e$\mu$ CR"}
        if self.diff == True:
            self.xlabels = {
                "ps_1_mass": r"$\Delta m_{K^{+}K^{-}}$ [GeV]",
                "d1_pt": r"Leading $K^{+}$ $\Delta p_{T}^{reco - truth}$ [GeV]",
                "d2_pt": r"Leading $K^{-}$ $\Delta p_{T}^{reco - truth}$ [GeV]",
                "d3_pt": r"Subleading $K^{+}$ $\Delta p_{T}^{reco - truth}$ [GeV]",
                "d4_pt": r"Subleading $K^{-}$ $\Delta p_{T}^{reco - truth}$ [GeV]",
                "d1_eta": r"Leading $K^{+}$ $\Delta\eta_{reco - truth}$",
                "d2_eta": r"Leading $K^{-}$ $\Delta\eta_{reco - truth}$",
                "d3_eta": r"Subleading $K^{+}$ $\Delta\eta_{reco - truth}$",
                "d4_eta": r"Subleading $K^{-}$ $\Delta\eta_{reco - truth}$",
                "H_eta": r"$\Delta\eta_{4K}^{reco - truth}$",
                "H_pt": r"$\Delta p_{T}^{4K}$ [GeV]",
            }
        elif self.append == True:
            self.xlabels = {
                "ps_1_mass": r"$m_{K^{+}K^{-}}$ [GeV]",
            }
        else: 
            self.xlabels = {
                "H_mass": r"$m_\mathrm{4K}$ [GeV]",
                "H_pt": r"$p_{T}^\mathrm{4K}$ [GeV]",
                "H_eta": r"$\eta_\mathrm{4K}$",
                "pt_vis": r"$p_{T}^{\ell\ell}$ [GeV]",
                "m_vis": r"$m_{\ell\ell}$ [GeV]",
                "eta_vis": r"Dilepton $\eta$",
                "d1_pt": r"Leading $K^{+}$ $p_{T}$ [GeV]",
                "d2_pt": r"Leading $K^{-}$ $p_{T}$ [GeV]",
                "d3_pt": r"Subleading $K^{+}$ $p_{T}$ [GeV]",
                "d4_pt": r"Subleading $K^{-}$ $p_{T}$ [GeV]",
                "d1_eta": r"Leading $K^{+}$ $\eta$",
                "d2_eta": r"Leading $K^{-}$ $\eta$",
                "d3_eta": r"Subleading $K^{+}$ $\eta$",
                "d4_eta": r"Subleading $K^{-}$ $\eta$",
                "d1_phi": r"Leading $K^{+}$ $\phi$",
                "d2_phi": r"Leading $K^{-}$ $\phi$",
                "d3_phi": r"Subleading $K^{+}$ $\phi$",
                "d4_phi": r"Subleading $K^{-}$ $\phi$",
                "d1_iso": r"Leading $K^{+}$ Isolation",
                "d2_iso": r"Leading $K^{-}$ Isolation",
                "d3_iso": r"Subleading $K^{+}$ Isolation",
                "d4_iso": r"Subleading $K^{-}$ Isolation",
                "njets": r"Number of jets",
                "pt_1": r"Leading lepton $p_{T}$ [GeV]",
                "pt_2": r"Subleading lepton $p_{T}$ [GeV]",
                "PV_npvs": r"Number of primary vertices",
                "PV_npvsGood": r"Number of primary vertices",
                "puweight": r"Pileup weight",
                "ps_mass": r"Mean $m_{KK}$ [GeV]",
                "ps_mass_12": r"$m_{K^{+}_{lead}K^{-}_{lead}}$ [GeV]",
                "ps_mass_14": r"$m_{K^{+}_{lead}K^{-}_{sublead}}$ [GeV]",
                "ps_mass_23": r"$m_{K^{-}_{lead}K^{+}_{sublead}}$ [GeV]",
                "ps_mass_34": r"$m_{K^{-}_{sublead}K^{+}_{sublead}}$ [GeV]",
                "ps_1_mass": r"$m_{K^{+}K^{-}}^{1st\ pair}$ [GeV]",
                "ps_2_mass": r"$m_{K^{+}K^{-}}^{2nd\ pair}$ [GeV]",
                "d1_prompt": r"Leading $K^{+}$ prompt flag",
                "d2_prompt": r"Leading $K^{-}$ prompt flag",
                "d3_prompt": r"Subleading $K^{+}$ prompt flag",
                "d4_prompt": r"Subleading $K^{-}$ prompt flag",
                }
    
    def beeeater_cmap(self):
        colors = ['seagreen', 'gold', 'darkorange', 'orange', 'mediumspringgreen']
        cmap = mcolors.ListedColormap(colors)
        return cmap

    def prep_ps_mass_mask(self, dict):
        comb = list(combinations(["12", "14", "23", "34"], 2))
        ps_masks = {}
        for key in dict.keys():
            mask = np.zeros(dict[key]["ps_mass_12"].to_numpy().shape[0], dtype=bool)
            for c in comb:
                ps_mass_1 = dict[key][f"ps_mass_{c[0]}"].to_numpy()
                ps_mass_2 = dict[key][f"ps_mass_{c[1]}"].to_numpy()
                mask |= (ps_mass_1 < 10) & (ps_mass_2 < 10) & (np.abs(ps_mass_1 - ps_mass_2) < 1)
            ps_masks[key] = mask
        return ps_masks
    
    def prep_prompt_mask(self, dict):
        prompt_masks = {}
        d_list = ["d1", "d2", "d3", "d4"]
        for key in dict.keys():
            masks = np.zeros([4,dict[key]["d1_prompt"].to_numpy().shape[0]], dtype=bool)
            for d in d_list:
                mask = dict[key][f"{d}_prompt"].to_numpy() == 3
                masks[int(d[-1])-1,:] = mask
                #print(masks)
            masks = np.sum(masks, axis=0) == 4
            prompt_masks[key] = masks
        return prompt_masks
    
    def prep_mc_data(self):
        mc_masks = {}
        mc_weights = {}
        mc_weights_up = {}
        mc_weights_down = {}
        for key, value in self.mc_branches.items():
            print(f"Preparing {key} samples")
            mask = np.ones(value['H_mass'].to_numpy().shape[0], dtype=bool)
            if self.selection == "loose":
                if ((self.scope == "mm") | (self.scope == "ee")) & (self.var == "H_mass"):
                    H_m_mask = ((value['H_mass'] > 60) & (value['H_mass'] < 220))

                if ((self.scope  == "mm") | (self.scope == "ee" )) & (self.var != "H_mass"):
                    H_m_mask = (((value['H_mass'] > 60) & (value['H_mass'] < 110)) | ((value['H_mass'] < 220) & (value['H_mass'] > 140)))

                if self.scope == "em":
                    H_m_mask = ((value['H_mass'] > 60) & (value['H_mass'] < 220))

                H_eta_mask = (np.abs(value['H_eta']) < 3)
                leading_k_mask = (value['d1_pt'] > 5) & (value['d2_pt'] > 5)
                z_pt_mask = (value['pt_vis'] > 20)
                z_eta_mask = (np.abs(value['eta_vis']) < 3)
                z_m_mask = ((value['m_vis'] > 70) & (value['m_vis'] < 110))
                if self.scope == "mm" or self.scope == "ee":
                    ps_mask = (value['ps_1_mass'] < 5) & (value['ps_2_mass'] < 5) & (np.abs(value['ps_1_mass'] - value['ps_2_mass']) < 1.5)
                else:
                    ps_mask = True
                mask = H_m_mask & H_eta_mask & z_m_mask  & ps_mask & z_pt_mask & leading_k_mask 
                mc_masks[key] = mask
            if self.selection == "tight":
                if ((self.scope == "mm") | (self.scope == "ee")) & (self.var == "H_mass"):
                    H_m_mask = ((value['H_mass'] > 70) & (value['H_mass'] < 200))

                if ((self.scope  == "mm") | (self.scope == "ee" )) & (self.var != "H_mass"):
                    H_m_mask = (((value['H_mass'] > 70) & (value['H_mass'] < 110)) | ((value['H_mass'] < 200) & (value['H_mass'] > 140)))

                H_eta_mask = (np.abs(value['H_eta']) < 2.4)
                leading_k_mask = (value['d1_pt'] > 10) & (value['d2_pt'] > 10)
                z_pt_mask = (value['pt_vis'] > 30)
                z_m_mask = ((value['m_vis'] > 75) & (value['m_vis'] < 105))
                #iso_mask = ((value['d1_iso'] < self.test_mask) & (value['d2_iso'] < self.test_mask) & (value['d3_iso'] < self.test_mask) & (value['d4_iso'] < self.test_mask))
                if self.scope == "mm" or self.scope == "ee":
                    ps_mask = (value['ps_1_mass'] < 3) & (value['ps_2_mass'] < 3) & (np.abs(value['ps_1_mass'] - value['ps_2_mass']) < 0.06)
                else:
                    raise ValueError("No tight mask for em")
                mask = H_m_mask & H_eta_mask & z_m_mask & ps_mask & z_pt_mask & leading_k_mask #& iso_mask
                mc_masks[key] = mask
            #if len(self.mc_branches) == 1:
            #    mask = np.ones(value['H_mass'].to_numpy().shape[0], dtype=bool)
            #    mc_masks[key] = mask

            if self.scope == "mm":
                trigger_wgt = 1-(1-value['trigger_wgt_mu_1'])*(1-value['trigger_wgt_mu_2'])
                #sfs = trigger_wgt[mask]*(value['id_wgt_mu_1'][mask]*value['iso_wgt_mu_1'][mask]+value['id_wgt_mu_2'][mask]*value['iso_wgt_mu_2'][mask])
                sfs = value['id_wgt_mu_1'][mask]*value['id_wgt_mu_2'][mask]*value['iso_wgt_mu_1'][mask]*value['iso_wgt_mu_2'][mask]*value['trigger_wgt_mu_1'][mask]*value['trigger_wgt_mu_2'][mask]
                signs = np.array(list(map(util.weights, value['genWeight'][mask])))
                mean = np.mean(value['puweight'][mask])
                if key == "Drell-Yan":
                    weight = self.lumi*signs*value['evtweight'][mask]*sfs*self.idy_zpt_weights[mask]*value['puweight'][mask]
                    weight_up = self.lumi*signs*value['evtweight'][mask]*sfs*self.idy_zpt_weights[mask]*value['puweight_up'][mask]
                    weight_down = self.lumi*signs*value['evtweight'][mask]*sfs*self.idy_zpt_weights[mask]*value['puweight_down'][mask]
                    mc_weights[key] = weight
                    mc_weights_up[key] = weight_up
                    mc_weights_down[key] = weight_down
                else:
                    weight = self.lumi*signs*sfs*value['evtweight'][mask]*value['puweight'][mask]
                    weight_up = self.lumi*signs*sfs*value['evtweight'][mask]*value['puweight_up'][mask]
                    weight_down = self.lumi*signs*sfs*value['evtweight'][mask]*value['puweight_down'][mask]
                    mc_weights[key] = weight
                    mc_weights_up[key] = weight_up
                    mc_weights_down[key] = weight_down

            if self.scope == "ee":
                sfs = value['id_wgt_ele_1'][mask]*value['id_wgt_ele_2'][mask]*value['trigger_wgt_ele_1'][mask]*value['trigger_wgt_ele_2'][mask]
                signs = np.array(list(map(util.weights, value['genWeight'][mask])))
                if key == "Drell-Yan":
                    weight = self.lumi*signs*value['evtweight'][mask]*self.idy_zpt_weights[mask]*sfs*value['puweight'][mask]
                    mc_weights[key] = weight
                else:
                    weight = self.lumi*signs*value['evtweight'][mask]*sfs*value['puweight'][mask]
                    mc_weights[key] = weight
            
            if self.scope == "em":
                sfs = value['id_wgt_ele_1'][mask]*value['id_wgt_mu_2'][mask]*value['iso_wgt_mu_2'][mask]*value['trigger_wgt_ele_1'][mask]*value['trigger_wgt_mu_2'][mask]
                signs = np.array(list(map(util.weights, value['genWeight'][mask])))
                if key == "Drell-Yan":
                    weight = self.lumi*signs*value['evtweight'][mask]*self.idy_zpt_weights[mask]*sfs*value['puweight'][mask]
                    mc_weights[key] = weight
                else:
                    weight = self.lumi*signs*value['evtweight'][mask]*sfs*value['puweight'][mask]
                    mc_weights[key] = weight
                
        return mc_masks, mc_weights, mc_weights_up, mc_weights_down
    
    def prep_signal_data(self):
        signal_masks = {}
        signal_weights = {}
        for key, value in self.signal_branches.items():
            mask = np.ones(value['H_mass'].to_numpy().shape[0], dtype=bool)
            if self.selection == "loose":
                H_m_mask = ((value['H_mass'] > 60) & (value['H_mass'] < 220))
                H_eta_mask = (np.abs(value['H_eta']) < 3)
                leading_k_mask = (value['d1_pt'] > 5) & (value['d2_pt'] > 5)
                z_pt_mask = (value['pt_vis'] > 20)
                z_eta_mask = (np.abs(value['eta_vis']) < 3)
                z_m_mask = ((value['m_vis'] > 70) & (value['m_vis'] < 110))
                if self.scope == "mm" or self.scope == "ee":
                    ps_mask = (value['ps_1_mass'] < 5) & (value['ps_2_mass'] < 5) & (np.abs(value['ps_1_mass'] - value['ps_2_mass']) < 1.5)
                else:
                    ps_mask = True
                mask = H_m_mask & z_m_mask & ps_mask & z_pt_mask & H_eta_mask & leading_k_mask
                signal_masks[key] = mask
            if self.selection == "tight":
                H_m_mask = ((value['H_mass'] > 70) & (value['H_mass'] < 200))
                H_eta_mask = (np.abs(value['H_eta']) < 2.4)
                leading_k_mask = (value['d1_pt'] > 10) & (value['d2_pt'] > 10)
                z_pt_mask = (value['pt_vis'] > 30)
                z_m_mask = ((value['m_vis'] > 75) & (value['m_vis'] < 105))
                #iso_mask = ((value['d1_iso'] < self.test_mask) & (value['d2_iso'] < self.test_mask) & (value['d3_iso'] < self.test_mask) & (value['d4_iso'] < self.test_mask))
                if self.scope == "mm" or self.scope == "ee":
                    ps_mask = (value['ps_1_mass'] < 3) & (value['ps_2_mass'] < 3) & (np.abs(value['ps_1_mass'] - value['ps_2_mass']) < 0.06)
                else:
                    raise ValueError("No tight mask for em")
                mask = H_m_mask & H_eta_mask & z_m_mask & ps_mask & z_pt_mask & leading_k_mask #& iso_mask
                signal_masks[key] = mask

            if self.scope == "mm":
                sfs = value['id_wgt_mu_1'][mask]*value['id_wgt_mu_2'][mask]*value['iso_wgt_mu_1'][mask]*value['iso_wgt_mu_2'][mask]*value['trigger_wgt_mu_1'][mask]*value['trigger_wgt_mu_2'][mask]
                signs = np.array(list(map(util.weights, value['genWeight'][mask])))
                pu_mean = np.mean(value['puweight'][mask])
                pu_normalized = value['puweight'][mask]/pu_mean
                weight = 0.1*self.lumi*value['evtweight'][mask]*signs*sfs*value['puweight'][mask]
                signal_weights[key] = weight
            if self.scope == "ee":
                sfs = value['id_wgt_ele_1'][mask]*value['id_wgt_ele_2'][mask]*value['trigger_wgt_ele_1'][mask]*value['trigger_wgt_ele_2'][mask]
                signs = np.array(list(map(util.weights, value['genWeight'][mask])))
                weight = self.lumi*value['evtweight'][mask]*signs*sfs*value['puweight'][mask]
                signal_weights[key] = weight
            if self.scope == "em":
                sfs = value['id_wgt_ele_1'][mask]*value['id_wgt_mu_2'][mask]*value['iso_wgt_mu_2'][mask]*value['trigger_wgt_ele_1'][mask]*value['trigger_wgt_mu_2'][mask]
                signs = np.array(list(map(util.weights, value['genWeight'][mask])))
                weight = signs*sfs*value['puweight'][mask]
                signal_weights[key] = weight
        return signal_masks, signal_weights
    
    def prep_data(self):
        data_masks = {}
        for key, value in self.data_branches.items():
            mask = np.ones(value['H_mass'].to_numpy().shape[0], dtype=bool)
            if (self.scope == "mm") | (self.scope == "ee"):
                if self.selection == "loose":
                    H_m_mask = (((value['H_mass'] > 60) & (value['H_mass'] < 110)) | ((value['H_mass'] > 140) & (value['H_mass'] < 220)))
                    H_eta_mask = (np.abs(value['H_eta']) < 3)
                    leading_k_mask = (value['d1_pt'] > 5) & (value['d2_pt'] > 5)
                    z_pt_mask = (value['pt_vis'] > 20)
                    z_eta_mask = (np.abs(value['eta_vis']) < 3)
                    z_m_mask = ((value['m_vis'] > 70) & (value['m_vis'] < 110))
                    ps_mask = (value['ps_1_mass'] < 5) & (value['ps_2_mass'] < 5) & (np.abs(value['ps_1_mass'] - value['ps_2_mass']) < 1.5)
                    mask = H_m_mask & z_m_mask & ps_mask  & z_pt_mask & H_eta_mask & leading_k_mask
                    data_masks[key] = mask
                if self.selection == "tight":
                    H_m_mask = (((value['H_mass'] > 70) & (value['H_mass'] < 110)) | ((value['H_mass'] > 140) & (value['H_mass'] < 200)))
                    H_eta_mask = (np.abs(value['H_eta']) < 2.4)
                    leading_k_mask = (value['d1_pt'] > 10) & (value['d2_pt'] > 10)
                    z_pt_mask = (value['pt_vis'] > 30)
                    z_m_mask = ((value['m_vis'] > 75) & (value['m_vis'] < 105))
                    #iso_mask = ((value['d1_iso'] < self.test_mask) & (value['d2_iso'] < self.test_mask) & (value['d3_iso'] < self.test_mask) & (value['d4_iso'] < self.test_mask))
                    ps_mask = (value['ps_1_mass'] < 3) & (value['ps_2_mass'] < 3) & (np.abs(value['ps_1_mass'] - value['ps_2_mass']) < 0.06)
                    mask = H_m_mask & H_eta_mask & z_m_mask & ps_mask & z_pt_mask & leading_k_mask #& iso_mask
                    data_masks[key] = mask

            if self.scope == "em":
                H_m_mask = ((value['H_mass'] > 60) & (value['H_mass'] < 220))
                H_eta_mask = (np.abs(value['H_eta']) < 3)
                z_pt_mask = (value['pt_vis'] > 20)
                z_m_mask = ((value['m_vis'] > 70) & (value['m_vis'] < 110))
                leading_k_mask = (value['d1_pt'] > 5) & (value['d2_pt'] > 5)
                mask = H_m_mask & z_pt_mask & z_m_mask & H_eta_mask & leading_k_mask #& self.data_prompt_masks[key]
                data_masks[key] = mask

            #if len(self.mc_branches) == 1:
            #    data_masks[key] = np.ones(value['H_mass'].to_numpy().shape[0], dtype=bool)

        return data_masks
    
    def plot_data_to_mc_full(self, save_as=None):
        fig, axs = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [3, 1]})

        bkg_m = [tree[self.var][self.mc_masks[key]].to_numpy() for key, tree in self.mc_branches.items()]
        bkg_weights = list(self.mc_weights.values())
        bkg_labels = list(self.mc_branches.keys())

        signal_m = [tree[self.var][self.signal_masks[key]].to_numpy() for key, tree in self.signal_branches.items()]
        signal_weights = list(self.signal_weights.values())
        signal_labels = list(self.signal_branches.keys())

        data_m = [tree[self.var][self.data_masks[key]].to_numpy() for key, tree in self.data_branches.items()]
        data_labels = list(self.data_branches.keys())

        bin_size = (self.upper - self.lower) / self.bins

        if self.var == "H_mass":

#            for i, (data, weights) in enumerate(zip(bkg_m, bkg_weights)):
#                hist, _ = np.histogram(data, bins=self.bins, range=(self.lower, self.upper), weights=weights)
#                event_count = np.sum(hist)
#                bkg_labels[i] = f"{bkg_labels[i]}, yield: {event_count:.2f}"

            axs[0].hist(
                bkg_m, 
                bins=np.linspace(self.lower, self.upper, self.bins + 1),
                weights=bkg_weights,
                stacked=True, 
                label=bkg_labels, 
                color=[self.beeeater(i/len(bkg_labels)) for i in range(len(bkg_labels))] #, 'lightcoral']
                )

            axs[0].hist(
                signal_m, 
                bins=np.linspace(self.lower, self.upper, self.bins + 1),
                weights=signal_weights,
                histtype='step', 
                label=signal_labels, 
                color=['red']
                )

            data_hist, data_edges = np.histogram(
                data_m, 
                bins=self.bins,
                range=(self.lower, self.upper)
                )

            bkg_hist, _ = np.histogram(np.concatenate(bkg_m), bins=self.bins, range=(self.lower, self.upper), weights=np.concatenate(bkg_weights))
            bkg_stat_error = np.sqrt(
                np.histogram(
                    np.concatenate(bkg_m), 
                    bins=self.bins, 
                    range=(self.lower, self.upper), 
                    weights=np.concatenate(bkg_weights)**2
                    )[0]
            )
        else:
 #           for i, (data, weights) in enumerate(zip(bkg_m, bkg_weights)):
 #               hist, _ = np.histogram(np.clip(data, self.lower - bin_size/2, self.upper + bin_size/2), 
 #                                      bins=self.bins +2, 
 #                                      range=(self.lower - bin_size, self.upper + bin_size), 
 #                                      weights=weights
 #                                      )
 #               event_count = np.sum(hist)
 #               bkg_labels[i] = f"{bkg_labels[i]}, yield: {event_count:.2f}"

            clipped_bkg_m = [np.clip(data, self.lower - bin_size/2, self.upper + bin_size/2) for data in bkg_m]

            n, bins, patches = axs[0].hist(
                clipped_bkg_m, 
                bins=np.linspace(self.lower - bin_size, self.upper + bin_size, self.bins + 3),
                weights=bkg_weights,
                stacked=True, 
                label=bkg_labels, 
                color=[self.beeeater(i/len(bkg_labels)) for i in range(len(bkg_labels))] #, 'lightcoral']
                )
            
            clipped_signal_m = [np.clip(data, self.lower - bin_size/2, self.upper + bin_size/2) for data in signal_m]

            axs[0].hist(
                clipped_signal_m,
                bins=np.linspace(self.lower - bin_size, self.upper + bin_size, self.bins + 3),
                weights=signal_weights,
                histtype='step', 
                label=signal_labels, 
                color=['red']
                )

            clipped_data_m = [np.clip(data, self.lower - bin_size/2, self.upper + bin_size/2) for data in data_m]

            data_hist, data_edges = np.histogram(
                clipped_data_m,
                bins=self.bins + 2,
                range=(self.lower - bin_size, self.upper + bin_size)
                )

            bkg_hist, _ = np.histogram(
                np.clip(np.concatenate(bkg_m), self.lower - bin_size/2, self.upper + bin_size/2),
                bins=self.bins + 2, 
                range=(self.lower - bin_size, self.upper + bin_size), 
                weights=np.concatenate(bkg_weights))
            bkg_stat_error = np.sqrt(
                np.histogram(
                    np.clip(np.concatenate(bkg_m), self.lower - bin_size/2, self.upper + bin_size/2), 
                    bins=self.bins + 2, 
                    range=(self.lower - bin_size, self.upper + bin_size), 
                    weights=np.concatenate(bkg_weights)**2
                    )[0]
            )

        bin_centers = 0.5 * (data_edges[:-1] + data_edges[1:])

        stat_error = np.sqrt(data_hist[data_hist > 0])
        
        data_int = np.sum(data_hist)
        data_labels = f"{data_labels[0]}"#, \nyield: {data_int:.2f}"

        axs[0].errorbar(bin_centers[data_hist > 0], data_hist[data_hist > 0], yerr=stat_error, linestyle='none', marker='.', color='black', label=data_labels)

        axs[0].set_title(r"$\it{" + "Private \: work" + "}$" +" "+ r"$\bf{" + "(CMS \: data/simulation)" + "}$",loc="left", fontsize=20)

        axs[0].set_title(f"{self.lumi/1000:.1f} fb$^{{-1}}$, {self.era} (13 $\mathrm{{{{TeV}}}}$)",loc="right", fontsize=20)

        #hep.cms.label('Private Work (CMS data/simulation)', data=True, year=self.era, lumi=self.lumi/1000, ax=axs[0])
        
        if "m" in self.var or "pt" in self.var:
            axs[0].set_ylabel(f'Events / {bin_size:.1f} GeV', fontsize=20)
        elif "iso" in self.var:
            axs[0].set_ylabel(f'Events / {bin_size:.2f}', fontsize=20)
        else:
            axs[0].set_ylabel(f'Events / {bin_size:.1f}', fontsize=20)

        obscured_patch = None

        if (self.scope == "mm" or self.scope == "ee") and (self.var == "H_eta"):
            name = 'Drell-Yan'
            for patch, label in zip(patches, bkg_labels):
                if label == name:
                    obscured_patch = patch
                    break
        elif (self.scope == "em") and (self.var == "H_eta" or self.var == "m_vis"):
            name = r'$t\,\bar{t}$'
            for patch, label in zip(patches, bkg_labels):
                if label == name:
                    obscured_patch = patch
                    break

        handles, labels = axs[0].get_legend_handles_labels()
        handles.append(mpatches.Patch(color='white', label=self.region[self.scope], alpha=0))
        #handles.append(mpatches.Patch(color='white', label=r'$60 < m_\mathrm{KKKK} < 220\, GeV$', alpha=0))
        #handles.append(mpatches.Patch(color='white', label=r'$|\eta_\mathrm{KKKK}| < 3$', alpha=0))
        #handles.append(mpatches.Patch(color='white', label=r'$70 < m_{\ell\ell} < 110\, GeV$', alpha=0))
        #handles.append(mpatches.Patch(color='white', label=r'$p_{T}^{\ell\ell} > 20\, GeV$', alpha=0))
        #handles.append(mpatches.Patch(color='white', label=r'Leading $K^{\pm}$ $p_{T} > 5\, GeV$', alpha=0))
        #if self.scope == "mm" or self.scope == "ee":
        #    handles.append(mpatches.Patch(color='white', label=r'$m_{K^+K^-} < 5\, GeV$', alpha=0))
        #    handles.append(mpatches.Patch(color='white', label=r'$|\Delta m_{K^+K^-, K^+K^-}| < 1.5\, GeV$', alpha=0))
        #handles.append(mpatches.Patch(color='white', label=r'Blinded Region:', alpha=0))
        #handles.append(mpatches.Patch(color='white', label=r'$110 < m_\mathrm{KKKK} < 140\, GeV$', alpha=0))
        
        handler_map = {}
        if obscured_patch is not None:
            for bar in obscured_patch:
                handler_map[bar] = HandlerPatchWithMargin()

        axs[0].legend(handles=handles, handler_map=handler_map, prop={'size': 20})
        
        if "iso" in self.var:
            axs[0].set_yscale('log')
            axs[0].set_ylim(top=2e8)
            axs[0].axvline(3, linestyle='--', color='black')
        else:
            axs[0].set_ylim(bottom=0)

        #axs[0].set_ylim(bottom=0)

        ratio = (data_hist[data_hist > 0] - bkg_hist[data_hist > 0]) / bkg_hist[data_hist > 0]
        stat_ratio_errors = np.abs(stat_error/bkg_hist[data_hist > 0])
        mc_ratio_errors = np.abs(data_hist[data_hist > 0]*bkg_stat_error[data_hist > 0]/bkg_hist[data_hist > 0]**2)
        axs[1].errorbar(bin_centers[data_hist > 0], ratio, yerr=stat_ratio_errors, linestyle='none', marker='.', color='black')
        axs[1].fill_between(bin_centers[data_hist > 0], - mc_ratio_errors, mc_ratio_errors, step='mid', facecolor='gray', edgecolor='white', alpha=0.2, label='MC stat. Error')
        axs[1].legend(prop={'size': 10})

        axs[1].axhline(0, linestyle=':', color='black')
        axs[1].axhline(.25, linestyle=':', color='black')
        axs[1].axhline(-.25, linestyle=':', color='black')
        axs[1].set_ylim(-0.5, 0.5)
        axs[1].set_yticks([-0.25, 0, 0.25])

        axs[1].set_ylabel(r'$\frac{\mathrm{Data} - \mathrm{MC}}{\mathrm{MC}}$')

        nonzero = np.where(data_hist > 0)[0]

        if nonzero.size > 0:
            axs[1].set_xlim(data_edges[:-1][data_hist > 0][0], data_edges[1:][data_hist > 0][-1])
        axs[1].set_xlabel(self.xlabels[self.var], fontsize=20)
        
        plt.tight_layout()
        fig.subplots_adjust(hspace=0)
        
        if save_as is not None:
            plt.savefig(f"{save_as}.png")
            plt.savefig(f"{save_as}.pdf")

        return fig, axs
    
    def plot_data_to_mc_no_signal(self, save_as=None):
        fig, axs = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [3, 1]})

        bkg_m = [tree[self.var][self.mc_masks[key]].to_numpy() for key, tree in self.mc_branches.items()]
        bkg_weights = list(self.mc_weights.values())
        bkg_weights_up = list(self.mc_weights_up.values())
        bkg_weights_down = list(self.mc_weights_down.values())
        bkg_labels = list(self.mc_branches.keys())

        data_m = [tree[self.var][self.data_masks[key]].to_numpy() for key, tree in self.data_branches.items()]
        data_labels = list(self.data_branches.keys())

        bin_size = (self.upper - self.lower) / self.bins

        if self.var == "H_mass":

#            for i, (data, weights) in enumerate(zip(bkg_m, bkg_weights)):
#                hist, _ = np.histogram(data, bins=self.bins, range=(self.lower, self.upper), weights=weights)
#                event_count = np.sum(hist)
#                bkg_labels[i] = f"{bkg_labels[i]}, yield: {event_count:.2f}"

            bkg_hist, _ = np.histogram(np.concatenate(bkg_m), bins=self.bins, range=(self.lower, self.upper), weights=np.concatenate(bkg_weights))

            axs[0].hist(
                bkg_m[0], 
                bins=np.linspace(self.lower, self.upper, self.bins + 1),
                weights=bkg_weights[0],
                stacked=True, 
                density=True,
                histtype='step',
                label=bkg_labels, 
                color=[self.beeeater(i/len(bkg_labels)) for i in range(len(bkg_labels))] #, 'lightcoral']
                )

            data_hist, data_edges = np.histogram(
                data_m, 
                bins=self.bins,
                density=True,
                range=(self.lower, self.upper)
                )
            
            bkg_stat_error = np.sqrt(
                np.histogram(
                    np.concatenate(bkg_m), 
                    bins=self.bins, 
                    range=(self.lower, self.upper), 
                    weights=np.concatenate(bkg_weights)**2
                    )[0]
            )

            bin_centers = 0.5 * (data_edges[:-1] + data_edges[1:])

        else:
 #           for i, (data, weights) in enumerate(zip(bkg_m, bkg_weights)):
 #               hist, _ = np.histogram(np.clip(data, self.lower - bin_size/2, self.upper + bin_size/2), 
 #                                      bins=self.bins +2, 
 #                                      range=(self.lower - bin_size, self.upper + bin_size), 
 #                                      weights=weights
 #                                      )
 #               event_count = np.sum(hist)
 #               bkg_labels[i] = f"{bkg_labels[i]}, yield: {event_count:.2f}"

            clipped_bkg_m = [np.clip(data, self.lower - bin_size/2, self.upper + bin_size/2) for data in bkg_m]

            bkg_hist, _ = np.histogram(
                np.clip(np.concatenate(bkg_m), self.lower - bin_size/2, self.upper + bin_size/2),
                bins=self.bins + 2, 
                range=(self.lower - bin_size, self.upper + bin_size), 
                weights=np.concatenate(bkg_weights))

            bkg_hist = bkg_hist.to_numpy()/np.sum(bkg_hist.to_numpy())

            bkg_hist_up, _ = np.histogram(
                np.clip(np.concatenate(bkg_m), self.lower - bin_size/2, self.upper + bin_size/2),
                bins=self.bins + 2, 
                range=(self.lower - bin_size, self.upper + bin_size), 
                weights=np.concatenate(bkg_weights_up))
            
            bkg_hist_up = bkg_hist_up.to_numpy()/np.sum(bkg_hist_up.to_numpy())

            bkg_hist_down, _ = np.histogram(
                np.clip(np.concatenate(bkg_m), self.lower - bin_size/2, self.upper + bin_size/2),
                bins=self.bins + 2, 
                range=(self.lower - bin_size, self.upper + bin_size), 
                weights=np.concatenate(bkg_weights_down))
            
            bkg_hist_down = bkg_hist_down.to_numpy()/np.sum(bkg_hist_down.to_numpy())

            clipped_data_m = [np.clip(data, self.lower - bin_size/2, self.upper + bin_size/2) for data in data_m]
        
            data_hist, data_edges = np.histogram(
                clipped_data_m,
                bins=self.bins + 2,
                range=(self.lower - bin_size, self.upper + bin_size),
                )
            
            data_hist = data_hist/np.sum(data_hist)

            print(type(data_hist))

            bkg_stat_error = np.sqrt(
                np.histogram(
                    np.clip(np.concatenate(bkg_m), self.lower - bin_size/2, self.upper + bin_size/2), 
                    bins=self.bins + 2, 
                    range=(self.lower - bin_size, self.upper + bin_size), 
                    weights=np.concatenate(bkg_weights)**2
                    )[0]
            )

            bin_centers = 0.5 * (data_edges[:-1] + data_edges[1:])
            axs[0].step(bin_centers[data_hist > 0], bkg_hist[data_hist > 0], where='mid', label='MC')#, color=[self.beeeater(i/len(bkg_labels)) for i in range(len(bkg_labels))])
            axs[0].fill_between(bin_centers[data_hist > 0], bkg_hist_down[data_hist > 0], bkg_hist_up[data_hist > 0], step='mid', facecolor='gray', edgecolor='white', alpha=0.2, label='MC stat. Error')
        
        print(type(bin_centers))

        stat_error = np.sqrt(data_hist[data_hist > 0])
        
        data_int = np.sum(data_hist)
        data_labels = f"{data_labels[0]}"#, \nyield: {data_int:.2f}"

        axs[0].errorbar(bin_centers[data_hist > 0], data_hist[data_hist > 0], yerr=0, linestyle='none', marker='.', color='black', label=data_labels)

        hep.cms.label('Private Work', data=True, year=self.era, lumi=self.lumi/1000, ax=axs[0])
        if "m" in self.var or "pt" in self.var:
            axs[0].set_ylabel(f'Events / {bin_size:.1f} GeV')
        else:
            axs[0].set_ylabel(f'Events / {bin_size:.1f}')

        #obscured_patch = None
#
        #if (self.scope == "mm" or self.scope == "ee") and self.var == "H_eta":
        #    name = 'Drell-Yan'
        #    for patch, label in zip(patches, bkg_labels):
        #        if label == name:
        #            obscured_patch = patch
        #            break
        #elif (self.scope == "em") and (self.var == "H_eta" or self.var == "m_vis"):
        #    name = r'$t\,\bar{t}$'
        #    for patch, label in zip(patches, bkg_labels):
        #        if label == name:
        #            obscured_patch = patch
        #            break
#
        #if isinstance(obscured_patch, plt.Rectangle):
        #    obscured_patch = [obscured_patch]

        handles, labels = axs[0].get_legend_handles_labels()
        #handles.append(mpatches.Patch(color='white', label=self.region[self.scope], alpha=0))
        #handles.append(mpatches.Patch(color='white', label=r'$60 < m_\mathrm{KKKK} < 220\, GeV$', alpha=0))
        #handles.append(mpatches.Patch(color='white', label=r'$|\eta_\mathrm{KKKK}| < 3$', alpha=0))
        #handles.append(mpatches.Patch(color='white', label=r'$70 < m_{\ell\ell} < 110\, GeV$', alpha=0))
        #handles.append(mpatches.Patch(color='white', label=r'$p_{T}^{\ell\ell} > 20\, GeV$', alpha=0))
        #handles.append(mpatches.Patch(color='white', label=r'Leading $K^{\pm}$ $p_{T} > 5\, GeV$', alpha=0))
        #if self.scope == "mm" or self.scope == "ee":
        #    handles.append(mpatches.Patch(color='white', label=r'$m_{K^+K^-} < 5\, GeV$', alpha=0))
        #    handles.append(mpatches.Patch(color='white', label=r'$|\Delta m_{K^+K^-, K^+K^-}| < 1.5\, GeV$', alpha=0))
        #handles.append(mpatches.Patch(color='white', label=r'Blinded Region:', alpha=0))
        #handles.append(mpatches.Patch(color='white', label=r'$110 < m_\mathrm{KKKK} < 140\, GeV$', alpha=0))
        
        #handler_map = {}
        #if obscured_patch is not None:
        #    for bar in obscured_patch:
        #        handler_map[bar] = HandlerPatchWithMargin()

        axs[0].legend(handles=handles) #, handler_map=handler_map)
        axs[0].set_ylim(bottom=0)

        ratio = (data_hist[data_hist > 0] - bkg_hist[data_hist > 0]) / bkg_hist[data_hist > 0]
        ratio_up = (data_hist[data_hist > 0] - bkg_hist_up[data_hist > 0]) / bkg_hist_up[data_hist > 0]
        ratio_down = (data_hist[data_hist > 0] - bkg_hist_down[data_hist > 0]) / bkg_hist_down[data_hist > 0]
        stat_ratio_errors = np.abs(stat_error/bkg_hist[data_hist > 0])
        mc_ratio_errors = np.abs(data_hist[data_hist > 0]*bkg_stat_error[data_hist > 0]/bkg_hist[data_hist > 0]**2)
        axs[1].errorbar(bin_centers[data_hist > 0], ratio, yerr=0, linestyle='none', marker='.', color='black')
        axs[1].fill_between(bin_centers[data_hist > 0], ratio_down, ratio_up, step='mid', facecolor='gray', edgecolor='white', alpha=0.2, label='MC stat. Error')
        #axs[1].fill_between(bin_centers[data_hist > 0], - mc_ratio_errors, mc_ratio_errors, step='mid', facecolor='gray', edgecolor='white', alpha=0.2, label='MC stat. Error')
        #axs[1].legend(prop={'size': 10})

        axs[1].axhline(0, linestyle=':', color='black')
        axs[1].axhline(.25, linestyle=':', color='black')
        axs[1].axhline(-.25, linestyle=':', color='black')
        axs[1].set_ylim(-0.5, 0.5)
        axs[1].set_yticks([-0.25, 0, 0.25])

        axs[1].set_ylabel(r'$\frac{\mathrm{Data} - \mathrm{MC}}{\mathrm{MC}}$')

        axs[1].set_xlim(data_edges[:-1][data_hist > 0][0], data_edges[1:][data_hist > 0][-1])
        axs[1].set_xlabel(self.xlabels[self.var])
        
        plt.tight_layout()
        fig.subplots_adjust(hspace=0)
        
        if save_as is not None:
            plt.savefig(f"{save_as}.png")
            plt.savefig(f"{save_as}.pdf")

        return fig, axs

    def plot_data_to_mc(self, save_as=None):

        if self.signal_branches == None and len(self.mc_branches) > 0:
            return self.plot_data_to_mc_no_signal(save_as)
        else:
            return self.plot_data_to_mc_full(save_as)

    def plot_mc(self, save_as=None):
        fig, axs = plt.subplots()
        #bkg_m = [tree[self.var][self.mc_masks[key]].to_numpy() for key, tree in self.mc_branches.items()]
        #bkg_weights = list(self.mc_weights.values())
        #bkg_labels = list(self.mc_branches.keys())

        signal_m = [tree[self.var][self.signal_masks[key]].to_numpy() for key, tree in self.signal_branches.items()]
        signal_weights = list(self.signal_weights.values())
        signal_labels = list(self.signal_branches.keys())

        bin_size = (self.upper - self.lower) / self.bins
        '''
        clipped_bkg_m = [np.clip(data, self.lower - bin_size/2, self.upper + bin_size/2) for data in bkg_m]

        for i, (data, weights) in enumerate(zip(bkg_m, bkg_weights)):
            hist, _ = np.histogram(np.clip(data, self.lower - bin_size/2, self.upper + bin_size/2), 
                                   bins=self.bins +2, 
                                   range=(self.lower - bin_size, self.upper + bin_size), 
                                   weights=weights
                                   )
            event_count = np.sum(hist)
            bkg_labels[i] = f"{bkg_labels[i]}, yield: {event_count:.2f}"

        hist, edges, _ = axs.hist(
                clipped_bkg_m, 
                bins=np.linspace(self.lower - bin_size, self.upper + bin_size, self.bins + 3),
                weights=bkg_weights,
                stacked=True, 
                label=bkg_labels, 
                color=['seagreen', 'gold', 'darkorange', 'orange', 'mediumspringgreen']
                )

        flattened_bkg_m = np.concatenate(bkg_m)
        flattened_bkg_weights = np.concatenate(bkg_weights)

        mean = np.average(flattened_bkg_m, weights=flattened_bkg_weights)
        '''
        clipped_signal_m = [np.clip(data, self.lower - bin_size/2, self.upper + bin_size/2) for data in signal_m]

        axs.hist(
            clipped_signal_m,
            bins=np.linspace(self.lower - bin_size, self.upper + bin_size, self.bins + 3),
            weights=signal_weights,
            histtype='step', 
            label=signal_labels, 
            color=['red']
            )
        
        hep.cms.label('', year=int(self.era), lumi=self.lumi/1000)

        axs.set_ylabel(f'Events / {bin_size:.2f}')
        axs.set_xlabel(self.xlabels[self.var])
        axs.set_xlim(edges[:-1][hist[-1,:] > 0][0], edges[1:][hist[-1,:] > 0][-1])

        handles, labels = axs.get_legend_handles_labels()
        handles.append(mpatches.Patch(color='white', label=self.region[self.scope], alpha=0))
        handles.append(mpatches.Patch(color='white', label=r'$70 < m_\mathrm{KKKK} < 200\, GeV$', alpha=0))
        handles.append(mpatches.Patch(color='white', label=r'$|\eta_\mathrm{KKKK}| < 2.4$', alpha=0))
        handles.append(mpatches.Patch(color='white', label=r'$75 < m_{\ell\ell} < 105\, GeV$', alpha=0))
        handles.append(mpatches.Patch(color='white', label=r'$p_{T}^{\ell\ell} > 30\, GeV$', alpha=0))
        handles.append(mpatches.Patch(color='white', label=r'Leading $K^{\pm}$ $p_{T} > 10\, GeV$', alpha=0))
        handles.append(mpatches.Patch(color='white', label=r'$m_{K^+K^-}^{lead} < 10\, GeV$, $m_{K^+K^-}^{sublead} < 10\, GeV$', alpha=0))
        handles.append(mpatches.Patch(color='white', label=r'$|\Delta m_{K^+K^-, K^+K^-}| < 1\, GeV$', alpha=0))
        handles.append(mpatches.Patch(color='white', label=r'Blinded Region:', alpha=0))
        handles.append(mpatches.Patch(color='white', label=r'$110 < m_\mathrm{KKKK} < 140\, GeV$', alpha=0))
        #handles.append(mpatches.Patch(color='white', label=f'Mean: {mean:.2f}', alpha=0))
        axs.legend(handles=handles, prop={'size': 12})
        #axs.set_ylim(bottom=0)

        plt.tight_layout()  

        if save_as is not None:
            plt.savefig(f"{save_as}.png")
            plt.savefig(f"{save_as}.pdf")

        return fig, axs
    
    def plot_reco_to_truth(self, save_as = None):
        fig, axs = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [3, 1]})

        signal_m_reco = [tree[self.var][self.signal_masks[key]].to_numpy() for key, tree in self.signal_branches.items()]
        signal_m_truth = [tree[f"truth_{self.var}"][self.signal_masks[key]].to_numpy() for key, tree in self.signal_branches.items()]
        signal_weights = list(self.signal_weights.values())

        bin_size = (self.upper - self.lower) / self.bins

        clipped_signal_m_reco = [np.clip(data, self.lower - bin_size/2, self.upper + bin_size/2) for data in signal_m_reco]
        clipped_signal_m_truth = [np.clip(data, self.lower - bin_size/2, self.upper + bin_size/2) for data in signal_m_truth]

        print(clipped_signal_m_reco[0].shape[0])

        reco_yield = [np.sum(weights) for weights in signal_weights]

        reco_hist, reco_edges, _ = axs[0].hist(
            clipped_signal_m_reco,
            bins=np.linspace(self.lower - bin_size, self.upper + bin_size, self.bins + 3),
            weights=signal_weights,
            label=f"Reco, yield: {reco_yield[0]:.2f}", 
            color=['seagreen']
            )
        
        truth_hist, truth_edges, _ = axs[0].hist(
            clipped_signal_m_truth,
            bins=np.linspace(self.lower - bin_size, self.upper + bin_size, self.bins + 3),
            weights=signal_weights,
            histtype='step', 
            label="Truth", 
            color=['red']
            )
        
        hep.cms.label('', year=int(self.era), lumi=self.lumi/1000, ax=axs[0])
        if "m" in self.var or "pt" in self.var:
            axs[0].set_ylabel(f'Events / {bin_size:.1f} GeV')
        else:
            axs[0].set_ylabel(f'Events / {bin_size:.1f}')
        handles, labels = axs[0].get_legend_handles_labels()
        handles.append(mpatches.Patch(color='white', label=self.region[self.scope], alpha=0))
        handles.append(mpatches.Patch(color='white', label=r'$70 < m_\mathrm{KKKK} < 200\, GeV$', alpha=0))
        handles.append(mpatches.Patch(color='white', label=r'$|\eta_\mathrm{KKKK}| < 2.4$', alpha=0))
        handles.append(mpatches.Patch(color='white', label=r'$75 < m_{\ell\ell} < 105\, GeV$', alpha=0))
        handles.append(mpatches.Patch(color='white', label=r'$p_{T}^{\ell\ell} > 30\, GeV$', alpha=0))
        handles.append(mpatches.Patch(color='white', label=r'Leading $K^{\pm}$ $p_{T} > 10\, GeV$', alpha=0))
        handles.append(mpatches.Patch(color='white', label=r'Blinded Region:', alpha=0))
        handles.append(mpatches.Patch(color='white', label=r'$122 < m_\mathrm{KKKK} < 128\, GeV$', alpha=0))
        axs[0].legend(handles=handles, prop={'size': 12})
        axs[0].set_ylim(bottom=0)

        bin_centers = 0.5 * (reco_edges[:-1] + reco_edges[1:])
        
        ratio = (reco_hist - truth_hist) / truth_hist

        axs[1].errorbar(bin_centers[reco_hist > 0], ratio[reco_hist > 0], linestyle='none', marker='.', color='black')
        axs[1].axhline(0, linestyle=':', color='black')
        axs[1].axhline(.25, linestyle=':', color='black')
        axs[1].axhline(-.25, linestyle=':', color='black')
        axs[1].set_ylim(-0.5, 0.5)
        axs[1].set_yticks([-0.25, 0, 0.25])

        axs[1].set_ylabel(r'$\frac{\mathrm{Reco} - \mathrm{Truth}}{\mathrm{Truth}}$')

        axs[1].set_xlim(truth_edges[:-1][truth_hist > 0][0], truth_edges[1:][truth_hist > 0][-1])

        axs[1].set_xlabel(self.xlabels[self.var])

        plt.tight_layout()
        fig.subplots_adjust(hspace=0)

        if save_as is not None:
            plt.savefig(f"{save_as}.png")
            plt.savefig(f"{save_as}.pdf")

        return fig, axs
    
    def plot_diff(self, save_as=None):
        fig, axs = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [3, 1]})

        bkg_var_1 = [tree[self.var_1][self.mc_masks[key]].to_numpy() for key, tree in self.mc_branches.items()]
        bkg_var_2 = [tree[self.var_2][self.mc_masks[key]].to_numpy() for key, tree in self.mc_branches.items()]
        bkg_weights = list(self.mc_weights.values())

        signal_var_1 = [tree[self.var_1][self.signal_masks[key]].to_numpy() for key, tree in self.signal_branches.items()]
        signal_var_2 = [tree[self.var_2][self.signal_masks[key]].to_numpy() for key, tree in self.signal_branches.items()]
        signal_weights = list(self.signal_weights.values())

        data_var_1 = [tree[self.var_1][self.data_masks[key]].to_numpy() for key, tree in self.data_branches.items()]
        data_var_2 = [tree[self.var_2][self.data_masks[key]].to_numpy() for key, tree in self.data_branches.items()]
        data_labels = list(self.data_branches.keys())
        
        bin_size = (self.upper - self.lower) / self.bins

        bkg_diff = [var_1 - var_2 for var_1, var_2 in zip(bkg_var_1, bkg_var_2)]
        bkg_diff = [np.clip(data, self.lower - bin_size/2, self.upper + bin_size/2) for data in bkg_diff]

        signal_diff = [var_1 - var_2 for var_1, var_2 in zip(signal_var_1, signal_var_2)]
        signal_diff = [np.clip(data, self.lower - bin_size/2, self.upper + bin_size/2) for data in signal_diff]

        data_diff = [var_1 - var_2 for var_1, var_2 in zip(data_var_1, data_var_2)]
        data_diff = np.array([np.clip(data, self.lower - bin_size/2, self.upper + bin_size/2) for data in data_diff])

        hist, edges, _ = axs[0].hist(
            bkg_diff,
            bins=np.linspace(self.lower - bin_size, self.upper + bin_size, self.bins + 3),
            weights=bkg_weights,
            stacked=True, 
            label=list(self.mc_branches.keys()), 
            color=['seagreen', 'gold', 'darkorange', 'orange', 'mediumspringgreen']
            )

        hist, edges, _ = axs[0].hist(
            signal_diff,
            bins=np.linspace(self.lower - bin_size, self.upper + bin_size, self.bins + 3),
            weights=signal_weights,
            histtype='step',
            label=list(self.signal_branches.keys()),
            color=['red']
            )
        
        data_hist, data_edges = np.histogram(
            data_diff[np.abs(data_diff) > 0.06],
            bins=self.bins + 2,
            range=(self.lower - bin_size, self.upper + bin_size)
            )
        
        bin_centers = 0.5 * (data_edges[:-1] + data_edges[1:])

        data_stat_error = np.sqrt(data_hist[data_hist > 0])

        bkg_hist, _ = np.histogram(
            np.clip(np.concatenate(bkg_diff), self.lower - bin_size/2, self.upper + bin_size/2),
            bins=self.bins + 2, 
            range=(self.lower - bin_size, self.upper + bin_size), 
            weights=np.concatenate(bkg_weights)
            )

        mc_stat_error=np.sqrt(
                np.histogram(
                    np.clip(np.concatenate(bkg_diff), self.lower - bin_size/2, self.upper + bin_size/2), 
                    bins=self.bins + 2, 
                    range=(self.lower - bin_size, self.upper + bin_size), 
                    weights=np.concatenate(bkg_weights)**2
                    )[0]
        )

        axs[0].errorbar(bin_centers[data_hist > 0], data_hist[data_hist > 0], yerr=data_stat_error, linestyle='none', marker='.', color='black', label=data_labels[0])
        #plt.errorbar(bin_centers[data_hist > 0], data_hist[data_hist > 0], yerr=mc_stat_error, linestyle='none', marker='none', color='black')
        
        axs[0].axvline(.06, linestyle='--', color='black')
        axs[0].axvline(-.06, linestyle='--', color='black')

        axs[0].set_title(r"$\it{" + "Private \: work" + "}$" +" "+ r"$\bf{" + "(CMS \: data/simulation)" + "}$",loc="left", fontsize=20)

        axs[0].set_title(f"{self.lumi/1000:.1f} fb$^{{-1}}$, {self.era} (13 $\mathrm{{{{TeV}}}}$)",loc="right", fontsize=20)

        if "m" in self.var or "pt" in self.var:
            axs[0].set_ylabel(f'Events / {bin_size:.3f} GeV', fontsize=20)
        else:
            axs[0].set_ylabel(f'Events / {bin_size:.1f}', fontsize=20)

        axs[0].set_xlabel(self.xlabels[self.var])
        axs[0].set_ylim(bottom=0, top=1.5*max(data_hist[data_hist > 0]))
        axs[0].set_xlim(edges[:-1][hist > 0][0], edges[1:][hist > 0][-1])

        handles, labels = axs[0].get_legend_handles_labels()
        handles.append(mpatches.Patch(color='white', label=self.region[self.scope], alpha=0))
        #handles.append(mpatches.Patch(color='white', label=r'$60 < m_\mathrm{KKKK} < 220\, GeV$', alpha=0))
        #handles.append(mpatches.Patch(color='white', label=r'$|\eta_\mathrm{KKKK}| < 3$', alpha=0))
        #handles.append(mpatches.Patch(color='white', label=r'$70 < m_{\ell\ell} < 110\, GeV$', alpha=0))
        #handles.append(mpatches.Patch(color='white', label=r'$p_{T}^{\ell\ell} > 20\, GeV$', alpha=0))
        #handles.append(mpatches.Patch(color='white', label=r'Leading $K^{\pm}$ $p_{T} > 5\, GeV$', alpha=0))
        #if self.scope == "mm" or self.scope == "ee":
        #    handles.append(mpatches.Patch(color='white', label=r'$m_{K^+K^-} < 5\, GeV$', alpha=0))
        #    handles.append(mpatches.Patch(color='white', label=r'$|\Delta m_{K^+K^-, K^+K^-}| < 1.5\, GeV$', alpha=0))
        #handles.append(mpatches.Patch(color='white', label=r'Blinded Region:', alpha=0))
        #handles.append(mpatches.Patch(color='white', label=r'$110 < m_\mathrm{KKKK} < 140\, GeV$', alpha=0))

        axs[0].legend(handles=handles, fontsize=12)

        ratio = (data_hist[data_hist > 0] - bkg_hist[data_hist > 0]) / bkg_hist[data_hist > 0]
        stat_ratio_errors = np.abs(data_stat_error/bkg_hist[data_hist > 0])
        mc_ratio_errors = np.abs(data_hist[data_hist > 0]*mc_stat_error[data_hist > 0]/bkg_hist[data_hist > 0]**2)
        axs[1].errorbar(bin_centers[data_hist > 0], ratio, yerr=stat_ratio_errors, linestyle='none', marker='.', color='black')
        axs[1].fill_between(bin_centers[data_hist > 0], - mc_ratio_errors, mc_ratio_errors, step='mid', facecolor='gray', edgecolor='white', alpha=0.2, label='MC stat. Error')
        axs[1].legend(prop={'size': 10})

        axs[1].axhline(0, linestyle=':', color='black')
        axs[1].axhline(.25, linestyle=':', color='black')
        axs[1].axhline(-.25, linestyle=':', color='black')
        axs[1].set_ylim(-0.5, 0.5)
        axs[1].set_yticks([-0.25, 0, 0.25])

        axs[1].set_ylabel(r'$\frac{\mathrm{Data} - \mathrm{MC}}{\mathrm{MC}}$', fontsize=20)

        axs[1].set_xlim(data_edges[:-1][data_hist > 0][0], data_edges[1:][data_hist > 0][-1])
        axs[1].set_xlabel(self.xlabels[self.var], fontsize=20)

        plt.tight_layout()
        fig.subplots_adjust(hspace=0)

        if save_as is not None:
            plt.savefig(f"{save_as}.png")
            plt.savefig(f"{save_as}.pdf")

        return fig, axs
    
    def plot_heatmap(self, xvar, yvar, save_as=None):
        fig, axs = plt.subplots()

        signal_x_reco = [tree[xvar][self.signal_masks[key]].to_numpy() for key, tree in self.signal_branches.items()]
        signal_x_truth = [tree[f"truth_{xvar}"][self.signal_masks[key]].to_numpy() for key, tree in self.signal_branches.items()]

        signal_x_diff = [reco - truth for reco, truth in zip(signal_x_reco, signal_x_truth)]

        signal_y_reco = [tree[yvar][self.signal_masks[key]].to_numpy() for key, tree in self.signal_branches.items()]
        signal_y_truth = [tree[f"truth_{yvar}"][self.signal_masks[key]].to_numpy() for key, tree in self.signal_branches.items()]

        signal_y_diff = [reco - truth for reco, truth in zip(signal_y_reco, signal_y_truth)]

        signal_weights = list(self.signal_weights.values())

        hist, x_edges, y_edges, im = axs.hist2d(
            np.concatenate(signal_x_diff),
            np.concatenate(signal_y_diff),
            bins=[np.linspace(self.xlower, self.xupper, self.bins), np.linspace(self.ylower, self.yupper, self.bins)],
            #weights=np.concatenate(signal_weights),
            cmap='viridis'
            )
        
        hep.cms.label('', year=int(self.era), lumi=self.lumi/1000, ax=axs)
        axs.set_xlabel(self.xlabels[xvar])
        axs.set_ylabel(self.xlabels[yvar])

        # Flatten the histogram and find the indices of non-zero elements
        non_zero_indices = np.flatnonzero(hist)

        # Get the first and last non-zero indices
        first_non_zero_idx = non_zero_indices[0]
        last_non_zero_idx = non_zero_indices[-1]

        # Calculate the bin edges corresponding to these indices
        first_bin_edge = x_edges[first_non_zero_idx // hist.shape[1]]
        last_bin_edge = x_edges[last_non_zero_idx // hist.shape[1] + 1]

        first_bin_edge_y = y_edges[first_non_zero_idx % hist.shape[1]]
        last_bin_edge_y = y_edges[last_non_zero_idx % hist.shape[1] + 1]


        # Set the x limits of the plot
        axs.set_xlim(first_bin_edge, last_bin_edge) #axs.set_xlim(x_edges[:-1][hist > 0][0], x_edges[1:][hist > 0][-1])
        axs.set_ylim(first_bin_edge_y, last_bin_edge_y) #axs.set_ylim(y_edges[:-1][hist > 0][0], y_edges[1:][hist > 0][-1])

        plt.colorbar(im, ax=axs)
        plt.tight_layout()

        if save_as is not None:
            plt.savefig(f"{save_as}.png")
            plt.savefig(f"{save_as}.pdf")

        return fig, axs
    
    def plot_2d(self, save_as=None):
        fig, ax = plt.subplots()

        signal_x_reco = [tree[self.var][self.signal_masks[key]].to_numpy() for key, tree in self.signal_branches.items()][0]
        signal_x_truth = [tree[f"truth_{self.var}"][self.signal_masks[key]].to_numpy() for key, tree in self.signal_branches.items()][0]
        signal_y = [tree[self.var2][self.signal_masks[key]].to_numpy() for key, tree in self.signal_branches.items()][0]
        signal_y_truth = [tree[f"truth_{self.var2}"][self.signal_masks[key]].to_numpy() for key, tree in self.signal_branches.items()][0]

        DeltaR = np.sqrt((signal_x_reco - signal_x_truth)**2 + (signal_y - signal_y_truth)**2)

        mean = np.mean(DeltaR[DeltaR > 0.4])
        frac = DeltaR[DeltaR > 0.4].shape[0]/DeltaR.shape[0]

        plt.plot(signal_x_truth[DeltaR > 0.4], signal_y_truth[DeltaR > 0.4], 'o', color="red", label='Truth')
        plt.plot(signal_x_reco[DeltaR > 0.4], signal_y[DeltaR > 0.4], 'o', color="blue", label='Reco')    

        handles, labels = ax.get_legend_handles_labels()
        handles.append(mpatches.Patch(color='white', label=self.region[self.scope], alpha=0))
        handles.append(mpatches.Patch(color='white', label=r'$70 < m_\mathrm{KKKK} < 200\, GeV$', alpha=0))
        handles.append(mpatches.Patch(color='white', label=r'$|\eta_\mathrm{KKKK}| < 2.4$', alpha=0))
        handles.append(mpatches.Patch(color='white', label=r'$75 < m_{\ell\ell} < 105\, GeV$', alpha=0))
        handles.append(mpatches.Patch(color='white', label=r'$p_{T}^{\ell\ell} > 30\, GeV$', alpha=0))
        handles.append(mpatches.Patch(color='white', label=r'$\Delta R$ > 0.4', alpha=0))
        handles.append(mpatches.Patch(color='white', label=fr'Mean $\Delta R$: {mean:.2f}', alpha=0))
        #ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.04, 1), borderaxespad=0, prop={'size': 12})

        hep.cms.label('', year=int(self.era), lumi=self.lumi/1000, ax=ax)
        ax.set_xlabel(self.xlabels[self.var])
        ax.set_ylabel(self.xlabels[self.var2])
        ax.set_xlim(self.xlower, self.xupper)
        ax.set_ylim(self.ylower, self.yupper)

        if save_as is not None:
            plt.savefig(f"{save_as}.png")
            plt.savefig(f"{save_as}.pdf")

        fig, ax = plt.subplots()

        DeltaR = np.clip(DeltaR, 0, 7.05)

        hist, edges, _ = ax.hist(
            DeltaR[DeltaR > 0.4],
            bins=np.linspace(0, 7.1, 72),
            color='seagreen'
            )
        print(edges)
        hep.cms.label('', year=int(self.era), lumi=self.lumi/1000, ax=ax)
        
        handles, labels = ax.get_legend_handles_labels()
        handles.append(mpatches.Patch(color='white', label=self.region[self.scope], alpha=0))
        handles.append(mpatches.Patch(color='white', label=r'$70 < m_\mathrm{KKKK} < 200\, GeV$', alpha=0))
        handles.append(mpatches.Patch(color='white', label=r'$|\eta_\mathrm{KKKK}| < 2.4$', alpha=0))
        handles.append(mpatches.Patch(color='white', label=r'$75 < m_{\ell\ell} < 105\, GeV$', alpha=0))
        handles.append(mpatches.Patch(color='white', label=r'$p_{T}^{\ell\ell} > 30\, GeV$', alpha=0))
        handles.append(mpatches.Patch(color='white', label=r'Leading $K^{\pm}$ $p_{T} > 10\, GeV$', alpha=0))
        handles.append(mpatches.Patch(color='white', label=r'$\Delta R$ > 0.4', alpha=0))
        handles.append(mpatches.Patch(color='white', label=fr'Mean $\Delta R$: {mean:.2f}', alpha=0))
        handles.append(mpatches.Patch(color='white', label=fr'Fraction $\Delta R$ > 0.4: {frac:.2f}', alpha=0))
        ax.legend(handles=handles, loc="upper right", prop={'size': 12})

        ax.set_xlabel(r'$\Delta R_{leading\  K^+}^{reco - truth}$')
        ax.set_xlim(0, 7.1)
        ax.set_ylabel('Events')

        return fig, ax
    
    def plot_mean(self, save_as=None):
        fig, ax = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [3, 1]})

        combinations = ["12", "14", "23", "34"]

        bkg_vars = {}
        for key, tree in self.mc_branches.items():
            arrays = [tree[self.var + "_" + comb][self.mc_masks[key]].to_numpy() for comb in combinations]
            bkg_vars[key] = np.stack(arrays)

        bkg_means = [np.mean(bkg_vars[keys], axis=0) for keys in bkg_vars.keys()]
        bkg_weights = list(self.mc_weights.values())
        bkg_labels = list(self.mc_branches.keys())

        signal_vars = {}
        for key, tree in self.signal_branches.items():
            arrays = [tree[self.var + "_" + comb][self.signal_masks[key]].to_numpy() for comb in combinations]
            signal_vars[key] = np.stack(arrays)

        signal_means = [np.mean(signal_vars[key], axis=0) for key in signal_vars.keys()]
        signal_weights = list(self.signal_weights.values())
        signal_labels = list(self.signal_branches.keys())

        data_vars = {}
        for key, tree in self.data_branches.items():
            arrays = [tree[self.var + "_" + comb][self.data_masks[key]].to_numpy() for comb in combinations]
            data_vars[key] = np.stack(arrays)

        data_means = [np.mean(data_vars[key], axis=0) for key in data_vars.keys()]
        data_labels = list(self.data_branches.keys())

        bin_size = (self.upper - self.lower) / self.bins

        for i, (data, weights) in enumerate(zip(bkg_means, bkg_weights)):
            hist, _ = np.histogram(
                np.clip(data, self.lower - bin_size/2, self.upper + bin_size/2), 
                bins=self.bins + 2, 
                range=(self.lower - bin_size, self.upper + bin_size), 
                weights=weights
                )
            event_count = np.sum(hist)
            bkg_labels[i] = f"{bkg_labels[i]}, yield: {event_count:.2f}"

        clipped_bkg_means = [np.clip(data, self.lower - bin_size/2, self.upper + bin_size/2) for data in bkg_means]
        clipped_signal_means = [np.clip(data, self.lower - bin_size/2, self.upper + bin_size/2) for data in signal_means]
        clipped_data_means = [np.clip(data, self.lower - bin_size/2, self.upper + bin_size/2) for data in data_means]

        ax[0].hist(
            clipped_bkg_means,
            bins=np.linspace(self.lower - bin_size, self.upper + bin_size, self.bins + 3),
            weights=bkg_weights,
            stacked=True, 
            label=bkg_labels, 
            color=['seagreen', 'gold', 'darkorange', 'orange', 'mediumspringgreen'] #, 'lightcoral']
            )
        
        ax[0].hist(
            clipped_signal_means,
            bins=np.linspace(self.lower - bin_size, self.upper + bin_size, self.bins + 3),
            weights=signal_weights,
            histtype='step', 
            label=signal_labels, 
            color=['red']
            )
        
        data_hist, data_edges = np.histogram(
            clipped_data_means, 
            bins=self.bins +2,
            range=(self.lower - bin_size, self.upper + bin_size)
            )

        bkg_hist, _ = np.histogram(np.concatenate(bkg_means), bins=self.bins, range=(self.lower, self.upper), weights=np.concatenate(bkg_weights))
        bkg_stat_error = np.sqrt(
            np.histogram(
                np.concatenate(bkg_means), 
                bins=self.bins, 
                range=(self.lower, self.upper), 
                weights=np.concatenate(bkg_weights)**2
                )[0]
        )

        bin_centers = 0.5 * (data_edges[:-1] + data_edges[1:])

        stat_error = np.sqrt(data_hist[data_hist > 0])

        data_int = np.sum(data_hist)

        data_labels = f"{data_labels[0]}, \nyield: {data_int:.2f}"

        ax[0].errorbar(bin_centers[data_hist > 0], data_hist[data_hist > 0], yerr=stat_error, linestyle='none', marker='.', color='black', label=data_labels)

        hep.cms.label('Work in Progress', data=True, year=int(self.era), lumi=self.lumi/1000, ax=ax[0])
        ax[0].set_ylabel(f'Events / {bin_size:.1f}')
        ax[0].set_ylim(bottom=0)
        handles, labels = ax[0].get_legend_handles_labels()
        handles.append(mpatches.Patch(color='white', label=self.region[self.scope], alpha=0))
        handles.append(mpatches.Patch(color='white', label=r'$70 < m_\mathrm{KKKK} < 200\, GeV$', alpha=0))
        handles.append(mpatches.Patch(color='white', label=r'$|\eta_\mathrm{KKKK}| < 2.4$', alpha=0))
        handles.append(mpatches.Patch(color='white', label=r'$75 < m_{\ell\ell} < 105\, GeV$', alpha=0))
        handles.append(mpatches.Patch(color='white', label=r'$p_{T}^{\ell\ell} > 30\, GeV$', alpha=0))
        handles.append(mpatches.Patch(color='white', label=r'Leading $K^{\pm}$ $p_{T} > 10\, GeV$', alpha=0))
        handles.append(mpatches.Patch(color='white', label=r'Two $a$ candidates needed with:', alpha=0))
        handles.append(mpatches.Patch(color='white', label=r'$m_{K^+K^-} < 10\, GeV$', alpha=0))
        handles.append(mpatches.Patch(color='white', label=r'$|\Delta m_{K^+K^-, K^+K^-}| < 1\, GeV$', alpha=0))
        handles.append(mpatches.Patch(color='white', label=r'Blinded Region:', alpha=0))
        handles.append(mpatches.Patch(color='white', label=r'$110 < m_\mathrm{KKKK} < 140\, GeV$', alpha=0))
        ax[0].legend(handles=handles, prop={'size': 12})

        ratio = (data_hist[data_hist > 0] - bkg_hist[data_hist > 0]) / bkg_hist[data_hist > 0]
        stat_ratio_errors = stat_error/bkg_hist[data_hist > 0]
        mc_ratio_errors = data_hist[data_hist > 0]*bkg_stat_error[data_hist > 0]/bkg_hist[data_hist > 0]**2
        ax[1].errorbar(bin_centers[data_hist > 0], ratio, yerr=stat_ratio_errors, linestyle='none', marker='.', color='black')
        ax[1].fill_between(bin_centers[data_hist > 0], - mc_ratio_errors, mc_ratio_errors, step='mid', facecolor='gray', edgecolor='white', alpha=0.2, label='MC Error')

        ax[1].axhline(0, linestyle=':', color='black')
        ax[1].axhline(.25, linestyle=':', color='black')
        ax[1].axhline(-.25, linestyle=':', color='black')

        ax[1].set_xlim(data_edges[:-1][data_hist > 0][0], data_edges[1:][data_hist > 0][-1])
        ax[1].set_ylim(-0.5, 0.5)
        ax[1].set_yticks([-0.25, 0, 0.25])

        ax[1].set_ylabel(r'$\frac{\mathrm{Data} - \mathrm{MC}}{\mathrm{MC}}$')
        ax[1].set_xlabel(self.xlabels[self.var])

        plt.tight_layout()
        fig.subplots_adjust(hspace=0)

        if save_as is not None:
            plt.savefig(f"{save_as}.png")
            plt.savefig(f"{save_as}.pdf")

        return fig, ax
    
    def plot_data_to_mc_wo_signal(self, save_as=None):

        fig, axs = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [3, 1]})

        bkg_m = [tree[self.var][self.mc_masks[key]].to_numpy() for key, tree in self.mc_branches.items()]
        bkg_weights = list(self.mc_weights.values())
        bkg_labels = list(self.mc_branches.keys())

        data_m = [tree[self.var][self.data_masks[key]].to_numpy() for key, tree in self.data_branches.items()]
        data_labels = list(self.data_branches.keys())

        bin_size = (self.upper - self.lower) / self.bins

        if self.var == "H_mass":

            for i, (data, weights) in enumerate(zip(bkg_m, bkg_weights)):
                hist, _ = np.histogram(data, bins=self.bins, range=(self.lower, self.upper), weights=weights)
                event_count = np.sum(hist)
                bkg_labels[i] = f"{bkg_labels[i]}, yield: {event_count:.2f}"

            axs[0].hist(
                bkg_m, 
                bins=np.linspace(self.lower, self.upper, self.bins + 1),
                weights=bkg_weights,
                stacked=True, 
                label=bkg_labels, 
                color=['seagreen', 'gold', 'darkorange', 'orange', 'mediumspringgreen'] #, 'lightcoral']
                )

            data_hist, data_edges = np.histogram(
                data_m, 
                bins=self.bins,
                range=(self.lower, self.upper)
                )

            bkg_hist, _ = np.histogram(np.concatenate(bkg_m), bins=self.bins, range=(self.lower, self.upper), weights=np.concatenate(bkg_weights))
            bkg_stat_error = np.sqrt(
                np.histogram(
                    np.concatenate(bkg_m), 
                    bins=self.bins, 
                    range=(self.lower, self.upper), 
                    weights=np.concatenate(bkg_weights)**2
                    )[0]
            )
        else:
            for i, (data, weights) in enumerate(zip(bkg_m, bkg_weights)):
                hist, _ = np.histogram(np.clip(data, self.lower - bin_size/2, self.upper + bin_size/2), 
                                       bins=self.bins +2, 
                                       range=(self.lower - bin_size, self.upper + bin_size), 
                                       weights=weights
                                       )
                event_count = np.sum(hist)
                bkg_labels[i] = f"{bkg_labels[i]}, yield: {event_count:.2f}"

            clipped_bkg_m = [np.clip(data, self.lower - bin_size/2, self.upper + bin_size/2) for data in bkg_m]

            axs[0].hist(
                clipped_bkg_m, 
                bins=np.linspace(self.lower - bin_size, self.upper + bin_size, self.bins + 3),
                weights=bkg_weights,
                stacked=True, 
                label=bkg_labels, 
                color=['seagreen', 'gold', 'darkorange', 'orange', 'mediumspringgreen'] #, 'lightcoral']
                )

            clipped_data_m = [np.clip(data, self.lower - bin_size/2, self.upper + bin_size/2) for data in data_m]

            data_hist, data_edges = np.histogram(
                clipped_data_m,
                bins=self.bins + 2,
                range=(self.lower - bin_size, self.upper + bin_size)
                )

            bkg_hist, _ = np.histogram(
                np.clip(np.concatenate(bkg_m), self.lower - bin_size/2, self.upper + bin_size/2),
                bins=self.bins + 2, 
                range=(self.lower - bin_size, self.upper + bin_size), 
                weights=np.concatenate(bkg_weights))
            bkg_stat_error = np.sqrt(
                np.histogram(
                    np.clip(np.concatenate(bkg_m), self.lower - bin_size/2, self.upper + bin_size/2), 
                    bins=self.bins + 2, 
                    range=(self.lower - bin_size, self.upper + bin_size), 
                    weights=np.concatenate(bkg_weights)**2
                    )[0]
            )

        bin_centers = 0.5 * (data_edges[:-1] + data_edges[1:])

        stat_error = np.sqrt(data_hist[data_hist > 0])
        
        data_int = np.sum(data_hist)
        data_labels = f"{data_labels[0]}, \nyield: {data_int:.2f}"

        axs[0].errorbar(bin_centers[data_hist > 0], data_hist[data_hist > 0], yerr=stat_error, linestyle='none', marker='.', color='black', label=data_labels)

        hep.cms.label('Work in Progress', data=True, year=int(self.era), lumi=self.lumi/1000, ax=axs[0])
        if "m" in self.var or "pt" in self.var:
            axs[0].set_ylabel(f'Events / {bin_size:.1f} GeV')
        else:
            axs[0].set_ylabel(f'Events / {bin_size:.1f}')
        handles, labels = axs[0].get_legend_handles_labels()
        handles.append(mpatches.Patch(color='white', label=self.region[self.scope], alpha=0))
        handles.append(mpatches.Patch(color='white', label=r'$70 < m_\mathrm{KKKK} < 200\, GeV$', alpha=0))
        handles.append(mpatches.Patch(color='white', label=r'$|\eta_\mathrm{KKKK}| < 2.4$', alpha=0))
        handles.append(mpatches.Patch(color='white', label=r'$75 < m_{\ell\ell} < 105\, GeV$', alpha=0))
        handles.append(mpatches.Patch(color='white', label=r'$p_{T}^{\ell\ell} > 30\, GeV$', alpha=0))
        handles.append(mpatches.Patch(color='white', label=r'Leading $K^{\pm}$ $p_{T} > 10\, GeV$', alpha=0))
        handles.append(mpatches.Patch(color='white', label=r'Two $a$ candidates needed with:', alpha=0))
        handles.append(mpatches.Patch(color='white', label=r'$m_{K^+K^-} < 10\, GeV$', alpha=0))
        handles.append(mpatches.Patch(color='white', label=r'$|\Delta m_{K^+K^-, K^+K^-}| < 1\, GeV$', alpha=0))
        handles.append(mpatches.Patch(color='white', label=r'Blinded Region:', alpha=0))
        handles.append(mpatches.Patch(color='white', label=r'$110 < m_\mathrm{KKKK} < 140\, GeV$', alpha=0))
        axs[0].legend(handles=handles, prop={'size': 12})
        axs[0].set_ylim(bottom=0)

        ratio = (data_hist[data_hist > 0] - bkg_hist[data_hist > 0]) / bkg_hist[data_hist > 0]
        stat_ratio_errors = stat_error/bkg_hist[data_hist > 0] 
        mc_ratio_errors = data_hist[data_hist > 0]*bkg_stat_error[data_hist > 0]/bkg_hist[data_hist > 0]**2
        axs[1].errorbar(bin_centers[data_hist > 0], ratio, yerr=stat_ratio_errors, linestyle='none', marker='.', color='black')
        axs[1].fill_between(bin_centers[data_hist > 0], - mc_ratio_errors, mc_ratio_errors, step='mid', facecolor='gray', edgecolor='white', alpha=0.2, label='MC Error')

        axs[1].axhline(0, linestyle=':', color='black')
        axs[1].axhline(.25, linestyle=':', color='black')
        axs[1].axhline(-.25, linestyle=':', color='black')
        axs[1].set_ylim(-0.5, 0.5)
        axs[1].set_yticks([-0.25, 0, 0.25])

        axs[1].set_ylabel(r'$\frac{\mathrm{Data} - \mathrm{MC}}{\mathrm{MC}}$')

        axs[1].set_xlim(data_edges[:-1][data_hist > 0][0], data_edges[1:][data_hist > 0][-1])
        axs[1].set_xlabel(self.xlabels[self.var])
        
        plt.tight_layout()
        fig.subplots_adjust(hspace=0)
        
        if save_as is not None:
            plt.savefig(f"{save_as}.png")
            plt.savefig(f"{save_as}.pdf")

        return fig, axs

    def plot_both(self, save_as=None):
        fig, axs = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [3, 1]})

        bkg_var_1 = [tree[self.var_1][self.mc_masks[key]].to_numpy() for key, tree in self.mc_branches.items()]
        bkg_var_2 = [tree[self.var_2][self.mc_masks[key]].to_numpy() for key, tree in self.mc_branches.items()]
        bkg_weights = list(self.mc_weights.values())

        signal_var_1 = [tree[self.var_1][self.signal_masks[key]].to_numpy() for key, tree in self.signal_branches.items()]
        signal_var_2 = [tree[self.var_2][self.signal_masks[key]].to_numpy() for key, tree in self.signal_branches.items()]
        signal_weights = list(self.signal_weights.values())

        data_var_1 = [tree[self.var_1][self.data_masks[key] & ((tree[self.var_1] < 1.4) | (tree[self.var_1] > 1.6))].to_numpy() for key, tree in self.data_branches.items()]
        data_var_2 = [tree[self.var_2][self.data_masks[key] & ((tree[self.var_2] < 1.4) | (tree[self.var_2] > 1.6))].to_numpy() for key, tree in self.data_branches.items()]
        data_labels = list(self.data_branches.keys())
        
        bin_size = (self.upper - self.lower) / self.bins

        bkg_diff = [np.append(var_1,var_2) for var_1, var_2 in zip(bkg_var_1, bkg_var_2)]
        bkg_diff = [np.clip(data, self.lower - bin_size/2, self.upper + bin_size/2) for data in bkg_diff]
        bkg_weights = [np.append(weight, weight) for weight in bkg_weights]

        signal_diff = [np.append(var_1,var_2) for var_1, var_2 in zip(signal_var_1, signal_var_2)]
        signal_diff = [np.clip(data, self.lower - bin_size/2, self.upper + bin_size/2) for data in signal_diff]
        signal_weights = [np.append(weight, weight) for weight in signal_weights]

        data_diff = [np.append(var_1,var_2) for var_1, var_2 in zip(data_var_1, data_var_2)]
        data_diff = [np.clip(data, self.lower - bin_size/2, self.upper + bin_size/2) for data in data_diff]

        hist, edges, _ = axs[0].hist(
            bkg_diff,
            bins=np.linspace(self.lower - bin_size, self.upper + bin_size, self.bins + 3),
            weights=bkg_weights,
            stacked=True, 
            label=list(self.mc_branches.keys()), 
            color=['seagreen', 'gold', 'darkorange', 'orange', 'mediumspringgreen']
            )

        hist, edges, _ = axs[0].hist(
            signal_diff,
            bins=np.linspace(self.lower - bin_size, self.upper + bin_size, self.bins + 3),
            weights=signal_weights,
            histtype='step',
            label=list(self.signal_branches.keys()),
            color=['red']
            )
        
        data_hist, data_edges = np.histogram(
            data_diff,
            bins=self.bins + 2,
            range=(self.lower - bin_size, self.upper + bin_size)
            )
        
        bin_centers = 0.5 * (data_edges[:-1] + data_edges[1:])

        data_stat_error = np.sqrt(data_hist[data_hist > 0])

        bkg_hist, _ = np.histogram(
            np.clip(np.concatenate(bkg_diff), self.lower - bin_size/2, self.upper + bin_size/2),
            bins=self.bins + 2, 
            range=(self.lower - bin_size, self.upper + bin_size), 
            weights=np.concatenate(bkg_weights)
            )

        mc_stat_error=np.sqrt(
                np.histogram(
                    np.clip(np.concatenate(bkg_diff), self.lower - bin_size/2, self.upper + bin_size/2), 
                    bins=self.bins + 2, 
                    range=(self.lower - bin_size, self.upper + bin_size), 
                    weights=np.concatenate(bkg_weights)**2
                    )[0]
        )

        axs[0].errorbar(bin_centers[data_hist > 0], data_hist[data_hist > 0], yerr=data_stat_error, linestyle='none', marker='.', color='black', label=data_labels[0])
        #plt.errorbar(bin_centers[data_hist > 0], data_hist[data_hist > 0], yerr=mc_stat_error, linestyle='none', marker='none', color='black')
        
        axs[0].set_title(r"$\it{" + "Private \: work" + "}$" +" "+ r"$\bf{" + "(CMS \: data/simulation)" + "}$",loc="left", fontsize=20)

        axs[0].set_title(f"{self.lumi/1000:.1f} fb$^{{-1}}$, {self.era} (13 $\mathrm{{{{TeV}}}}$)",loc="right", fontsize=20)

        if "m" in self.var or "pt" in self.var:
            axs[0].set_ylabel(f'Events / {bin_size:.1f} GeV', fontsize=20)
        else:
            axs[0].set_ylabel(f'Events / {bin_size:.1f}', fontsize=20)

        axs[0].set_xlabel(self.xlabels[self.var])
        axs[0].set_ylim(bottom=0)
        axs[0].set_xlim(edges[:-1][hist > 0][0], edges[1:][hist > 0][-1])

        handles, labels = axs[0].get_legend_handles_labels()
        handles.append(mpatches.Patch(color='white', label=self.region[self.scope], alpha=0))
        #handles.append(mpatches.Patch(color='white', label=r'$60 < m_\mathrm{KKKK} < 220\, GeV$', alpha=0))
        #handles.append(mpatches.Patch(color='white', label=r'$|\eta_\mathrm{KKKK}| < 3$', alpha=0))
        #handles.append(mpatches.Patch(color='white', label=r'$70 < m_{\ell\ell} < 110\, GeV$', alpha=0))
        #handles.append(mpatches.Patch(color='white', label=r'$p_{T}^{\ell\ell} > 20\, GeV$', alpha=0))
        #handles.append(mpatches.Patch(color='white', label=r'Leading $K^{\pm}$ $p_{T} > 5\, GeV$', alpha=0))
        #if self.scope == "mm" or self.scope == "ee":
        #    handles.append(mpatches.Patch(color='white', label=r'$m_{K^+K^-} < 5\, GeV$', alpha=0))
        #    handles.append(mpatches.Patch(color='white', label=r'$|\Delta m_{K^+K^-, K^+K^-}| < 1.5\, GeV$', alpha=0))
        #handles.append(mpatches.Patch(color='white', label=r'Blinded Region:', alpha=0))
        #handles.append(mpatches.Patch(color='white', label=r'$110 < m_\mathrm{KKKK} < 140\, GeV$', alpha=0))

        axs[0].legend(handles=handles,prop={'size': 20})

        ratio = (data_hist[data_hist > 0] - bkg_hist[data_hist > 0]) / bkg_hist[data_hist > 0]
        stat_ratio_errors = np.abs(data_stat_error/bkg_hist[data_hist > 0])
        mc_ratio_errors = np.abs(data_hist[data_hist > 0]*mc_stat_error[data_hist > 0]/bkg_hist[data_hist > 0]**2)
        axs[1].errorbar(bin_centers[data_hist > 0], ratio, yerr=stat_ratio_errors, linestyle='none', marker='.', color='black')
        axs[1].fill_between(bin_centers[data_hist > 0], - mc_ratio_errors, mc_ratio_errors, step='mid', facecolor='gray', edgecolor='white', alpha=0.2, label='MC stat. Error')
        axs[1].legend(prop={'size': 10})

        axs[1].axhline(0, linestyle=':', color='black')
        axs[1].axhline(.25, linestyle=':', color='black')
        axs[1].axhline(-.25, linestyle=':', color='black')
        axs[1].set_ylim(-0.5, 0.5)
        axs[1].set_yticks([-0.25, 0, 0.25])

        axs[1].set_ylabel(r'$\frac{\mathrm{Data} - \mathrm{MC}}{\mathrm{MC}}$', fontsize=20)

        axs[1].set_xlim(data_edges[:-1][data_hist > 0][0], data_edges[1:][data_hist > 0][-1])
        axs[1].set_xlabel(self.xlabels[self.var_1], fontsize=20)

        plt.tight_layout()
        fig.subplots_adjust(hspace=0)

        if save_as is not None:
            plt.savefig(f"{save_as}.png")
            plt.savefig(f"{save_as}.pdf")

        return fig, axs
    def check_sfs(self):
        systematics = {
            "pu": ["puweight", "puweight_up", "puweight_down"],
            "id_1": ["id_wgt_mu_1", "id_wgt_mu_1__MuonIDUp", "id_wgt_mu_1__MuonIDDown"] if self.scope == "mm" else ["id_wgt_ele_1", "id_wgt_ele_1__ElectronIDUp", "id_wgt_ele_1__ElectronIDDown"],
            "id_2": ["id_wgt_mu_2", "id_wgt_mu_2__MuonIDUp", "id_wgt_mu_2__MuonIDDown"] if self.scope == "mm" else ["id_wgt_ele_2", "id_wgt_ele_2__ElectronIDUp", "id_wgt_ele_2__ElectronIDDown"],
            "iso_1": ["iso_wgt_mu_1", "iso_wgt_mu_1__MuonIsoUp", "iso_wgt_mu_1__MuonIsoDown"] if self.scope == "mm" else ["1", "1"],
            "iso_2": ["iso_wgt_mu_2", "iso_wgt_mu_2__MuonIsoUp", "iso_wgt_mu_2__MuonIsoDown"] if self.scope == "mm" else ["1", "1"],
            "trigger_1": ["trigger_wgt_mu_1", "trigger_wgt_mu_1__MuonTriggerUp", "trigger_wgt_mu_1__MuonTriggerDown"] if self.scope == "mm" else ["trigger_wgt_ele_1", "trigger_wgt_ele_1_up", "trigger_wgt_ele_1_down"],
            "trigger_2": ["trigger_wgt_mu_2", "trigger_wgt_mu_2__MuonTriggerUp", "trigger_wgt_mu_2__MuonTriggerDown"] if self.scope == "mm" else ["trigger_wgt_ele_2", "trigger_wgt_ele_2_up", "trigger_wgt_ele_2_down"]
        }

        signal_m = [tree[self.var][self.signal_masks[key]].to_numpy() for key, tree in self.signal_branches.items()]
        signal_weights = list(self.signal_weights.values())

        nominal_histogram, edges = np.histogram(
            signal_m,
            bins=self.bins,
            range=(self.lower, self.upper),
            weights=signal_weights
        )

        for key, values in systematics.items():
            if self.scope == 'ee' and "iso" in key:
                continue
            
            sys_nominal = [tree[values[0]][self.signal_masks[label]].to_numpy() for label, tree in self.signal_branches.items()]
            sys_up = [tree[values[1]][self.signal_masks[label]].to_numpy() for label, tree in self.signal_branches.items()]
            sys_down = [tree[values[2]][self.signal_masks[label]].to_numpy() for label, tree in self.signal_branches.items()]

            print(f"{key} Nominal: {sys_nominal[0]}")
            print(f"{key} Up: {sys_up[0]}")
            print(f"{key} Down: {sys_down[0]}")

            weights_sys_up = signal_weights[0]/sys_nominal[0]*sys_up[0]
            weights_sys_down = signal_weights[0]/sys_nominal[0]*sys_down[0]

            fig, axs = plt.subplots()   
            hist_nom, edges, _ = plt.hist(
                signal_m,
                bins=np.linspace(self.lower, self.upper, self.bins + 1),
                weights=signal_weights,
                histtype='step',
                label="Nominal"
            )

            hist_up, edges, _ = plt.hist(
                signal_m[0],
                bins=np.linspace(self.lower, self.upper, self.bins + 1),
                weights=weights_sys_up,
                histtype='step',
                label=f"{key} Up"
            )
            hist_down, edges, _ = plt.hist(
                signal_m[0],
                bins=np.linspace(self.lower, self.upper, self.bins + 1),
                weights=weights_sys_down,
                histtype='step',
                label=f"{key} Down"
            )

            plt.xlim(edges[:-1][nominal_histogram > 0][0], edges[1:][nominal_histogram > 0][-1])
            plt.xlabel(self.xlabels[self.var])
            plt.ylabel('Events')
            #plt.ylim(bottom=0)
            plt.legend()

            max_nominal = np.max(hist_nom)
            max_up = np.max(hist_up)
            max_down = np.max(hist_down)

            plt.figtext(0.15, 0.8, f"Nominal: {max_nominal:.2f}", fontsize=12)
            plt.figtext(0.15, 0.75, f"{key} Up: {max_up:.2f}", fontsize=12)
            plt.figtext(0.15, 0.7, f"{key} Down: {max_down:.2f}", fontsize=12)

            axs.set_title(r"$\it{" + "Private \: work" + "}$" +" "+ r"$\bf{" + "(CMS \: data/simulation)" + "}$",loc="left", fontsize=20)
            axs.set_title(f"{self.lumi/1000:.1f} fb$^{{-1}}$, {self.era} (13 $\mathrm{{{{TeV}}}}$)",loc="right", fontsize=20)

            if max_up < max_nominal or max_down > max_nominal:
                plt.savefig(f"/web/jhornung/public_html/analysis_plots/systematics_{self.scope}_{key}.png")
                plt.savefig(f"/web/jhornung/public_html/analysis_plots/systematics_{self.scope}_{key}.pdf")

        plt.show()

        return fig, axs

if __name__ == "__main__":

    mc_loc = "/ceph/jhornung/MC_2018/2018"
    scope = "mm"
    bkgs = ["DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8+RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2+MINIAODSIM",
            "ttbar",
            "diboson",
            "wh",
            "single_top"
            ]
    bkg_labels = ['Drell-Yan', 
                  r'$t\,\bar{t}$', 
                  r'Diboson', 
                  r'$W\,H\rightarrow \ell\,\nu\, b\,\bar{b}$', 
                  r'Single Top'
                  ]
    mc_bkg = data(mc_loc, bkgs, bkg_labels, scope)
    bkg_branches = mc_bkg.branches_dict
    
    signal = ["HaaKKKK_M125_13TeV_powheg_pythia8-RunIISummer20UL18MiniAODv2-106X"]
    signal_label = [r'$H\rightarrow aa\rightarrow 4\,K$']
    mc_signal = data(mc_loc, signal, signal_label, scope)
    signal_branches = mc_signal.branches_dict
    
    data_loc = "/ceph/jhornung/Data_2018"
    datasamples = ["single_muon_data"]
    data_label = ["Single Muon Data"]
    muons = data(data_loc, datasamples, data_label, scope)
    muons_branches = muons.branches_dict
    
    config = {
        "mc": mc_bkg.branches_dict,
        "signal": mc_signal.branches_dict,
        "data": muons.branches_dict,
        "scope": scope,
        "z_pt_weights": mc_bkg.get_Z_pt_weights(),
        "var": "pt_vis",
        "lower": 30,
        "upper": 300,
        "bins": 45
    }
    
    mass = hist(config)
    fig, axs = mass.plot_data_to_mc()
    plt.show()