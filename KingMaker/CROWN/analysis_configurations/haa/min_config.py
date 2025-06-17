from __future__ import annotations  # needed for type annotations in > python 3.7

from typing import List

from .producers import event as event
from .producers import genparticles as genparticles
from .producers import muons as muons
from .producers import electrons as electrons
from .producers import photons as photons
from .producers import pairquantities as pairquantities
from .producers import pairselection as pairselection
from .producers import higgs as higgs
from .producers import higgsdaughters as higgsdaughters
from .producers import jets as jets
from .producers import fatjets as fatjets
from .producers import pfcands as pfcands
from .producers import scalefactors as scalefactors
from .producers import fromnano as fromnano
from .quantities import nanoAOD as nanoAOD
from .quantities import output as q
from code_generation.configuration import Configuration
from code_generation.modifiers import EraModifier
from code_generation.rules import RemoveProducer
from code_generation.systematics import SystematicShift


def build_config(
    era: str,
    sample: str,
    scopes: List[str],
    shifts: List[str],
    available_sample_types: List[str],
    available_eras: List[str],
    available_scopes: List[str],
):
    configuration = Configuration(
        era,
        sample,
        scopes,
        shifts,
        available_sample_types,
        available_eras,
        available_scopes,
    )

    configuration.add_config_parameters(
        "global",
        {
            "PU_reweighting_file": EraModifier(
                {
                    "2016preVFP": "data/jsonpog-integration/POG/LUM/2016preVFP_UL/puWeights.json.gz",
                    "2016postVFP": "data/jsonpog-integration/POG/LUM/2016postVFP_UL/puWeights.json.gz",
                    "2017": "data/jsonpog-integration/POG/LUM/2017_UL/puWeights.json.gz",
                    "2018": "data/jsonpog-integration/POG/LUM/2018_UL/puWeights.json.gz",
                }
            ),
            "PU_reweighting_era": EraModifier(
                {
                    "2016preVFP": "Collisions16_UltraLegacy_goldenJSON",
                    "2016postVFP": "Collisions16_UltraLegacy_goldenJSON",
                    "2017": "Collisions17_UltraLegacy_goldenJSON",
                    "2018": "Collisions18_UltraLegacy_goldenJSON",
                }
            ),
            "golden_json_file": EraModifier(
                {
                    "2016preVFP": "data/golden_json/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt",
                    "2016postVFP": "data/golden_json/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt",
                    "2017": "data/golden_json/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt",
                    "2018": "data/golden_json/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt",
                }
            ),
            "met_filters": EraModifier(
                {
                    "2016preVFP": [
                        "Flag_goodVertices",
                        "Flag_globalSuperTightHalo2016Filter",
                        "Flag_HBHENoiseFilter",
                        "Flag_HBHENoiseIsoFilter",
                        "Flag_EcalDeadCellTriggerPrimitiveFilter",
                        "Flag_BadPFMuonFilter",
                        # "Flag_BadPFMuonDzFilter", # only since nanoAODv9 available
                        "Flag_eeBadScFilter",
                    ],
                    "2016postVFP": [
                        "Flag_goodVertices",
                        "Flag_globalSuperTightHalo2016Filter",
                        "Flag_HBHENoiseFilter",
                        "Flag_HBHENoiseIsoFilter",
                        "Flag_EcalDeadCellTriggerPrimitiveFilter",
                        "Flag_BadPFMuonFilter",
                        # "Flag_BadPFMuonDzFilter", # only since nanoAODv9 available
                        "Flag_eeBadScFilter",
                    ],
                    "2017": [
                        "Flag_goodVertices",
                        "Flag_globalSuperTightHalo2016Filter",
                        "Flag_HBHENoiseFilter",
                        "Flag_HBHENoiseIsoFilter",
                        "Flag_EcalDeadCellTriggerPrimitiveFilter",
                        "Flag_BadPFMuonFilter",
                        # "Flag_BadPFMuonDzFilter", # only since nanoAODv9 available
                        "Flag_eeBadScFilter",
                        "Flag_ecalBadCalibFilter",
                    ],
                    "2018": [
                        "Flag_goodVertices",
                        "Flag_globalSuperTightHalo2016Filter",
                        "Flag_HBHENoiseFilter",
                        "Flag_HBHENoiseIsoFilter",
                        "Flag_EcalDeadCellTriggerPrimitiveFilter",
                        "Flag_BadPFMuonFilter",
                        # "Flag_BadPFMuonDzFilter", # only since nanoAODv9 available
                        "Flag_eeBadScFilter",
                        "Flag_ecalBadCalibFilter",
                    ],
                }
            ),
        },
    )

    # base selection:
    configuration.add_config_parameters(
        "global",
        {
            "min_muon_pt": 0.0,
#            "max_muon_eta": 2.4,
#            "max_muon_dxy": 0.045,
#            "max_muon_dz": 0.2,
            "muon_id": "Muon_looseId",
#            "muon_iso_cut": 0.3,
            "min_electron_pt": 0.0,
#            "max_electron_eta": 2.4,
#            "max_electron_dxy": 0.045,
#            "max_electron_dz": 0.2,
            "electron_id": "Electron_cutBased",
            "electron_id_wp": 2,
#            "electron_iso_cut": 0.3,
            "photon_id": "Photon_cutBased",
            "photon_id_wp": 2,
            
            "pfcands_pdgid": "211,-211",
            "min_pfcands_pt": 1.,
            "max_pfcands_eta": 2.5,
            "min_pfcands_mass": 0.139,
            "max_pfcands_mass": 0.14,
            #"pfcands_fromPV_flag": "3",
#            "higgs_reco_pt_cut": 10.0,

            "charged_pfcands_pdgid": "211,-211",
            "neutral_pfcands_pdgid": "130",
            "photon_pfcands_pdgid": "22",
            "fromPV": "3",

            #"kaon_mass": 0.493677,

            "genpart_pdgid": "23",
            "genpart_status": "62",

#            "min_jet_pt": 0.0,#20.0,
#            "max_jet_eta": 2.4,
            "btag_cut": 0.6,
            "deltaR_jet_veto": 0.4,

            "min_fatjet_pt": 0.0,
        },
    )
    # MM scope Muon selection
    configuration.add_config_parameters(
        ["mm"],
        {
            "muon_index_in_pair": 0,
            "second_muon_index_in_pair": 1,
            "min_muon_pt": 23.0,
            "max_muon_eta": 2.1,
            "max_muon_dxy": 0.045,
            "max_muon_dz": 0.2,
            "muon_id": "Muon_tightId",
            "muon_iso_cut": 0.15,
            "truegen_mother_pdgid": 23,
            "truegen_daughter_1_pdgid": 13,
            "truegen_daughter_2_pdgid": 13,
        },
    )
    
    # EE scope Electron selection
    configuration.add_config_parameters(
        ["ee"],
        {
            "electron_index_in_pair": 0,
            "second_electron_index_in_pair": 1,
            "min_electron_pt": 23.0,
            "max_electron_eta": 2.1,
            "max_electron_dxy": 0.045,
            "max_electron_dz": 0.2,
            "electron_id": "Electron_cutBased",
            "electron_id_wp": 3,
            "electron_iso_cut": 0.15,
            "truegen_mother_pdgid": 23,
            "truegen_daughter_1_pdgid": 11,
            "truegen_daughter_2_pdgid": 11,
        },
    )

    # EM scope Electron selection
    configuration.add_config_parameters(
        ["em"],
        {
            "min_electron_pt": 23.0,
            "max_electron_eta": 2.1,
            "max_electron_dxy": 0.045,
            "max_electron_dz": 0.2,
            "electron_id": "Electron_cutBased",
            "electron_id_wp": 3,
            "electron_iso_cut": 0.15,
            "min_muon_pt": 23.0,
            "max_muon_eta": 2.1,
            "max_muon_dxy": 0.045,
            "max_muon_dz": 0.2,
            "muon_id": "Muon_mediumId",
            "muon_iso_cut": 0.15,
        },
    )
    
    # Muon scale factors configuration
    configuration.add_config_parameters(
        ["mm","em"],
        {
            "muon_sf_file": EraModifier(
                {
                    "2016preVFP": "data/jsonpog-integration/POG/MUO/2016preVFP_UL/muon_Z.json.gz",
                    "2016postVFP": "data/jsonpog-integration/POG/MUO/2016postVFP_UL/muon_Z.json.gz",
                    "2017": "data/jsonpog-integration/POG/MUO/2017_UL/muon_Z.json.gz",
                    "2018": "data/jsonpog-integration/POG/MUO/2018_UL/muon_Z.json.gz",
                }
            ),
            "muon_id_sf_name": "NUM_TightID_DEN_TrackerMuons",
            "muon_iso_sf_name": "NUM_TightRelIso_DEN_TightIDandIPCut",
            "muon_trigger_sf_name": EraModifier(
                {
                    "2016preVFP": "NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdTight_and_PFIsoTight",
                    "2016postVFP": "NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdTight_and_PFIsoTight",
                    "2018": "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight",
                    "2017": "NUM_IsoMu27_DEN_CutBasedIdTight_and_PFIsoTight",
                }
            ),
            "muon_sf_year_id": EraModifier(
                {
                    "2016preVFP": "2016preVFP_UL",
                    "2016postVFP": "2016postVFP_UL",
                    "2017": "2017_UL",
                    "2018": "2018_UL",
                }
            ),
            "muon_sf_varation": "sf",  # "sf" is nominal, "systup"/"systdown" are up/down variations
        },
    )    
    
    # electron scale factors configuration
    configuration.add_config_parameters(
        ["ee", "em"],
        {
            "ele_sf_file": EraModifier(
                {
                    "2016preVFP": "data/jsonpog-integration/POG/EGM/2016preVFP_UL/electron.json.gz",
                    "2016postVFP": "data/jsonpog-integration/POG/EGM/2016postVFP_UL/electron.json.gz",
                    "2017": "data/jsonpog-integration/POG/EGM/2017_UL/electron.json.gz",
                    "2018": "data/jsonpog-integration/POG/EGM/2018_UL/electron.json.gz",
                }
            ),
            "ele_id_sf_name": "UL-Electron-ID-SF",
            "ele_sf_year_id": EraModifier(
                {
                    "2016preVFP": "2016preVFP",
                    "2016postVFP": "2016postVFP",
                    "2017": "2017",
                    "2018": "2018",
                }
            ),
            "ele_sf_varation": "sf",  # "sf" is nominal, "sfup"/"sfdown" are up/down variations
            "nom_ele_trigger_sf_file": EraModifier(
                {
                    "2016preVFP": "data/custom_top_sf/electron/2016preVFP_UL/trigger_2016preVFP.json.gz",
                    "2016postVFP": "data/custom_top_sf/electron/2016postVFP_UL/trigger_2016postVFP.json.gz",
                    "2017": "data/custom_top_sf/electron/2017_UL/trigger_2017.json.gz",
                    "2018": "data/custom_top_sf/electron/2018_UL/trigger_2018.json.gz",
                }
            ),
            "syst_ele_trigger_sf_file": EraModifier(
                {
                    "2016preVFP": "data/custom_top_sf/electron/2016preVFP_UL/trigger_2016preVFP_syststat.json.gz",
                    "2016postVFP": "data/custom_top_sf/electron/2016postVFP_UL/trigger_2016postVFP_syststat.json.gz",
                    "2017": "data/custom_top_sf/electron/2017_UL/trigger_2017_syststat.json.gz",
                    "2018": "data/custom_top_sf/electron/2018_UL/trigger_2018_syststat.json.gz",
                }
            ),
            "nom_ele_trigger_sf_name": "h2_scaleFactorsEGamma",
            "syst_ele_trigger_sf_name": "h2_uncertaintiesEGamma",
        },
    )

    ## all scopes misc settings
    configuration.add_config_parameters(
        scopes,
        {
            "deltaR_jet_veto": 0.4,
            "pairselection_min_dR": 0.5,
        },
    )

    configuration.add_producers(
        "global",
        [
            event.SampleFlags,
            event.PUweights,
            event.PUweights_up,
            event.PUweights_down,
            event.Lumi,
            event.MetFilter,
            muons.BaseMuons,
            electrons.BaseElectrons,
            photons.BasePhotons,
            jets.GoodJets,
            jets.JetCollection,
            jets.BasicJetQuantities,
            fatjets.FatJetCollection,
            fatjets.BasicFatJetQuantities,
            pfcands.BasePFCands,
            pfcands.ChargedPFCands,
            pfcands.NeutralPFCands,
            pfcands.PhotonPFCands,
            pfcands.fromPV
            #higgsdaughters.GetTrueDaughterP4s,
            #higgsdaughters.truth_dQuantities,
#            fromnano.copy,
            #simplejets.GoodJets,
            #simplejets.NumberOfJets,
            #simplepfcands.GoodPFCands,
            #simplepfcands.NumberOfPFCands,
        ],
    )
    if sample == "data":
        configuration.add_producers(
            "global",
            [
                event.JSONFilter,
            ],
        )
    
    if era == "2016preVFP" or era == "2016postVFP":
        configuration.add_modification_rule(
            "global",
            RemoveProducer(
                producers=[
                    jets.BasicJetQuantities,
                    fatjets.BasicFatJetQuantities,
                ],
                samples=["data", "Haa", "bkg"],
            ),
        )

    configuration.add_producers(
        scopes,
        [
            higgs.ChargePairs,
            higgs.GoodHiggsDaughtersFilter,
            higgs.HiggsQuantities,
            higgsdaughters.dQuantities,
            higgs.ps_quanities,
            higgs.GetPS,
            higgsdaughters.ps_daughters,
            higgs.ps_masses,
        ],
    )
    configuration.add_producers(
        "mm",
        [
            muons.GoodMuons,
            muons.NumberOfGoodMuons,
            pairselection.ZMMPairSelection,
            pairselection.GoodMMPairFilter,
            pairselection.LVMu1,
            pairselection.LVMu2,
            pairquantities.MMDiTauPairQuantities,
            genparticles.MMGenDiTauPairQuantities,
            scalefactors.MuonIDIsoTrigger_SF,
        ],
    )
    
    configuration.add_producers(
        "ee",
        [
            electrons.GoodElectrons,
            electrons.NumberOfGoodElectrons,
            pairselection.ZEEPairSelection,
            pairselection.GoodEEPairFilter,
            pairselection.LVEl1,
            pairselection.LVEl2,
            pairquantities.EEDiTauPairQuantities,
            genparticles.EEGenDiTauPairQuantities,
            scalefactors.ElectronIDTrigger_SF
        ],
    )

    configuration.add_producers(
        "em",
        [
            electrons.GoodElectrons,
            electrons.NumberOfGoodElectrons,
            muons.GoodMuons,
            muons.NumberOfGoodMuons,
            pairselection.EMPairSelection,
            pairselection.GoodEMPairFilter,
            pairselection.LVEl1,
            pairselection.LVMu2,
            pairquantities.EMDiTauPairQuantities,
            genparticles.EMGenDiTauPairQuantities,
            scalefactors.ElectronIDTrigger_SF,
            scalefactors.MuonIDIsoTrigger_SF
        ],
    )
    configuration.add_modification_rule(
        "global",
        RemoveProducer(
            producers=[
                event.PUweights,
                event.PUweights_up,
                event.PUweights_down,
                ],
            samples=["data"],
        ),
    )
    configuration.add_modification_rule(
        "mm",
        RemoveProducer(
            producers=[
                genparticles.MMGenDiTauPairQuantities,
                scalefactors.MuonIDIsoTrigger_SF,
            ],
            samples=["data"],
        ),
    )
    configuration.add_modification_rule(
        "em",
        RemoveProducer(
            producers=[
                genparticles.EMGenDiTauPairQuantities,
                scalefactors.ElectronIDTrigger_SF,
                scalefactors.MuonIDIsoTrigger_SF   
            ],
            samples=["data"],
        ),
    )
    configuration.add_modification_rule(
        "ee",
        RemoveProducer(
            producers=[
                genparticles.EEGenDiTauPairQuantities,
                scalefactors.ElectronIDTrigger_SF
            ],
            samples=["data"],
        ),
    )
    configuration.add_outputs(
        "global",
	[
        nanoAOD.PV_npvs,
        nanoAOD.PV_npvsGood,
        q.H_pt,
	    q.H_eta,
	    q.H_phi,
	    q.H_mass,
        #q.highPtPair,
        #q.higgsdaughters,
        q.ps_1_mass,
        q.ps_2_mass,
        #q.ps1Pair,
        #q.ps2Pair,
        #q.higgsdaughters,
        #q.ps_mass_12,
        #q.ps_mass_14,
        #q.ps_mass_23,
        #q.ps_mass_34,
        #q.truth_H_pt,
        #q.truth_H_eta,
        #q.pfcands_iso,
        #q.pfcands_pdgid,
        #q.pfcands_pt,
        #q.pfcands_eta,
        #q.pfcands_phi,
        #q.pfcands_mass,
        #q.genpart_pdgid,
        #q.genpart_status,
        #q.genpart_eta,
        #q.genpart_phi,
        #q.daughter_isos,
	    q.d1_pt,
	    q.d1_eta,
	    q.d1_phi,
        q.d1_prompt,
        q.d1_iso,
	    q.d2_pt,
	    q.d2_eta,
	    q.d2_phi,
        q.d2_prompt,
        q.d2_iso,
	    q.d3_pt,
	    q.d3_eta,
	    q.d3_phi,
        q.d3_prompt,
        q.d3_iso,
	    q.d4_pt,
	    q.d4_eta,
	    q.d4_phi,
        q.d4_prompt,
        q.d4_iso,
	    #q.truth_d1_pt,
	    #q.truth_d1_eta,
	    #q.truth_d1_phi,
	    #q.truth_d2_pt,
	    #q.truth_d2_eta,
	    #q.truth_d2_phi,
	    #q.truth_d3_pt,
	    #q.truth_d3_eta,
	    #q.truth_d3_phi,
	    #q.truth_d4_pt,
	    #q.truth_d4_eta,
	    #q.truth_d4_phi,
        q.njets,
        q.jpt_1,
        q.jeta_1,
        q.jphi_1,
        q.jmass_1,
        q.jbtagDeepB_1,
        q.jbtagDeepFlavQG_1,
        q.jchEmEF_1,
        q.jchHEF_1,
        q.jneEmEF_1,
        q.jneHEF_1,
        q.jpt_2,
        q.jeta_2,
        q.jphi_2,
        q.jmass_2,
        q.jbtagDeepB_2,
        q.jbtagDeepFlavQG_2,
        q.jchEmEF_2,
        q.jchHEF_2,
        q.jneEmEF_2,
        q.jneHEF_2,
        q.jpt_3,
        q.jeta_3,
        q.jphi_3,
        q.jmass_3,
        q.jbtagDeepB_3,
        q.jbtagDeepFlavQG_3,
        q.jchEmEF_3,
        q.jchHEF_3,
        q.jneEmEF_3,
        q.jneHEF_3,
        q.jpt_4,
        q.jeta_4,
        q.jphi_4,
        q.jmass_4,
        q.jbtagDeepB_4,
        q.jbtagDeepFlavQG_4,
        q.jchEmEF_4,
        q.jchHEF_4,
        q.jneEmEF_4,
        q.jneHEF_4,
        q.jpt_5,
        q.jeta_5,
        q.jphi_5,
        q.jmass_5,
        q.jbtagDeepB_5,
        q.jbtagDeepFlavQG_5,
        q.jchEmEF_5,
        q.jchHEF_5,
        q.jneEmEF_5,
        q.jneHEF_5, 
	    q.nfatjets,
        q.fjpt_1,
        q.fjeta_1,
        q.fjphi_1,
        q.fjmass_1,
        q.fjmsoftdrop_1,
        q.fjn2b1_1,
        q.fjn3b1_1,
        q.fjparticleNet_H4qvsQCD_1,
        q.fjparticleNet_QCD_1,
        q.fjparticleNet_mass_1,
        q.fjpt_2,
        q.fjeta_2,
        q.fjphi_2,
        q.fjmass_2,
        q.fjmsoftdrop_2,
        q.fjn2b1_2,
        q.fjn3b1_2,
        q.fjparticleNet_H4qvsQCD_2,
        q.fjparticleNet_QCD_2,
        q.fjparticleNet_mass_2,
        q.fjpt_3,
        q.fjeta_3,
        q.fjphi_3,
        q.fjmass_3,
        q.fjmsoftdrop_3,
        q.fjn2b1_3,
        q.fjn3b1_3,
        q.fjparticleNet_H4qvsQCD_3,
        q.fjparticleNet_QCD_3,
        q.fjparticleNet_mass_3,
        ],
    )

    if sample != "data":
        configuration.add_outputs(
            "global",
            [
                nanoAOD.genWeight,
            ],
        )
	 	
    configuration.add_outputs(
       "mm",
        [
            q.is_data,
            q.is_embedding,
            q.is_ttbar,
            q.is_dyjets,
            q.is_wjets,
            q.is_diboson,
            nanoAOD.run,
            q.lumi,
            nanoAOD.event,
            q.puweight,
            q.puweight_up,
            q.puweight_down,
            q.pt_1,
            q.pt_2,
            q.eta_1,
            q.eta_2,
            q.phi_1,
            q.phi_2,
            q.m_vis,
            q.pt_vis,
            q.eta_vis,
            q.phi_vis,
            q.gen_pt_1,
            q.gen_eta_1,
            q.gen_phi_1,
            q.gen_mass_1,
            #q.gen_pdgid_1,
            q.gen_pt_2,
            q.gen_eta_2,
            q.gen_phi_2,
            q.gen_mass_2,
            #q.gen_pdgid_2,
            q.gen_m_vis,
            q.gen_pt_vis,
            q.id_wgt_mu_1,
            q.id_wgt_mu_2,
            q.iso_wgt_mu_1,
            q.iso_wgt_mu_2,
            q.trigger_wgt_mu_1,
            q.trigger_wgt_mu_2,
        ],
    )
    configuration.add_outputs(
        "ee",
        [
            q.is_data,
            q.is_embedding,
            q.is_ttbar,
            q.is_dyjets,
            q.is_wjets,
            q.is_diboson,
            nanoAOD.run,
            q.lumi,
            nanoAOD.event,
            q.puweight,
            q.puweight_up,
            q.puweight_down,
            q.pt_1,
            q.pt_2,
            q.eta_1,
            q.eta_2,
            q.phi_1,
            q.phi_2,
            q.m_vis,
            q.pt_vis,
            q.eta_vis,
            q.phi_vis,
            q.gen_pt_1,
            q.gen_eta_1,
            q.gen_phi_1,
            q.gen_mass_1,
            #q.gen_pdgid_1,
            q.gen_pt_2,
            q.gen_eta_2,
            q.gen_phi_2,
            q.gen_mass_2,
            #q.gen_pdgid_2,
            q.gen_m_vis,
            q.gen_pt_vis,
            q.id_wgt_ele_1,
            q.id_wgt_ele_2,
            q.trigger_wgt_ele_1,
            q.trigger_wgt_ele_1_up,
            q.trigger_wgt_ele_1_down,
            q.trigger_wgt_ele_2,
            q.trigger_wgt_ele_2_up,
            q.trigger_wgt_ele_2_down
            ],
    )
    configuration.add_outputs(
        "em",
        [
            q.is_data,
            q.is_embedding,
            q.is_ttbar,
            q.is_dyjets,
            q.is_wjets,
            q.is_diboson,
            nanoAOD.run,
            q.lumi,
            nanoAOD.event,
            q.puweight,
            q.puweight_up,
            q.puweight_down,
            q.pt_1,
            q.pt_2,
            q.eta_1,
            q.eta_2,
            q.phi_1,
            q.phi_2,
            q.m_vis,
            q.pt_vis,
            q.eta_vis,
            q.phi_vis,
            q.gen_pt_1,
            q.gen_eta_1,
            q.gen_phi_1,
            q.gen_mass_1,
            #q.gen_pdgid_1,
            q.gen_pt_2,
            q.gen_eta_2,
            q.gen_phi_2,
            q.gen_mass_2,
            #q.gen_pdgid_2,
            q.gen_m_vis,
            q.gen_pt_vis,
            q.gen_eta_vis,
            q.id_wgt_ele_1,
            q.id_wgt_mu_2,
            q.iso_wgt_mu_2,
            q.trigger_wgt_ele_1,
            q.trigger_wgt_ele_1_up,
            q.trigger_wgt_ele_1_down,
            q.trigger_wgt_mu_2
            ],
    )

    configuration.add_shift(
        SystematicShift(
            name="MuonIDUp",
            shift_config={"mm": {"muon_sf_varation": "systup"}},
            producers={
                "mm": [
                    scalefactors.Muon_1_ID_SF,
                    scalefactors.Muon_2_ID_SF,
                ],
                "em": [
                    scalefactors.Muon_2_ID_SF,
                ],
            },
        )
    )
    configuration.add_shift(
        SystematicShift(
            name="MuonIDDown",
            shift_config={"mm": {"muon_sf_varation": "systdown"}},
            producers={
                "mm": [
                    scalefactors.Muon_1_ID_SF,
                    scalefactors.Muon_2_ID_SF,
                ],
                "em": [
                    scalefactors.Muon_2_ID_SF,
                ],
            },
        )
    )

    configuration.add_shift(
        SystematicShift(
            name="MuonIsoUp",
            shift_config={"mm": {"muon_sf_varation": "systup"}},
            producers={
                "mm": [
                    scalefactors.Muon_1_Iso_SF,
                    scalefactors.Muon_2_Iso_SF,
                ],
                "em": [
                    scalefactors.Muon_2_Iso_SF,
                ],
            },
        )
    )

    configuration.add_shift(
        SystematicShift(
            name="MuonIsoDown",
            shift_config={"mm": {"muon_sf_varation": "systdown"}},
            producers={
                "mm": [
                    scalefactors.Muon_1_Iso_SF,
                    scalefactors.Muon_2_Iso_SF,
                ],
                "em": [
                    scalefactors.Muon_2_Iso_SF,
                ],
            },
        )
    )

    configuration.add_shift(
        SystematicShift(
            name="MuonTriggerUp",
            shift_config={"mm": {"muon_sf_varation": "systup"}},
            producers={
                "mm": [
                    scalefactors.Muon_1_Trigger_SF,
                    scalefactors.Muon_2_Trigger_SF,
                ],
                "em": [
                    scalefactors.Muon_2_Trigger_SF,
                ],
            },
        )
    )

    configuration.add_shift(
        SystematicShift(
            name="MuonTriggerDown",
            shift_config={"mm": {"muon_sf_varation": "systdown"}},
            producers={
                "mm": [
                    scalefactors.Muon_1_Trigger_SF,
                    scalefactors.Muon_2_Trigger_SF,
                ],
                "em": [
                    scalefactors.Muon_2_Trigger_SF,
                ],
            },
        )
    )

    configuration.add_shift(
        SystematicShift(
            name="ElectronIDUp",
            shift_config={"ee": {"ele_sf_varation": "sfup"}},
            producers={
                "ee": [
                    scalefactors.Ele_1_ID_SF,
                    scalefactors.Ele_2_ID_SF,
                ],
                "em": [
                    scalefactors.Ele_1_ID_SF,
                ],
            },
        )
    )

    configuration.add_shift(
        SystematicShift(
            name="ElectronIDDown",
            shift_config={"ee": {"ele_sf_varation": "sfdown"}},
            producers={
                "ee": [
                    scalefactors.Ele_1_ID_SF,
                    scalefactors.Ele_2_ID_SF,
                ],
                "em": [
                    scalefactors.Ele_1_ID_SF,
                ],
            },
        )
    )


    #########################
    # Finalize and validate the configuration
    #########################
    configuration.optimize()
    configuration.validate()
    configuration.report()
    return configuration.expanded_configuration()
