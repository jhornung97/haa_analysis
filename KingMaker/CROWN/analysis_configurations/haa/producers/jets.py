from ..quantities import output as q
from ..quantities import nanoAOD as nanoAOD
from code_generation.producer import Producer, ProducerGroup

####################
# Set of producers used for selection possible good jets
####################
'''
JetPtCorrection = Producer(
    name="JetPtCorrection",
    call="physicsobject::jet::JetPtCorrection({df}, {output}, {input}, {jet_reapplyJES}, {jet_jes_sources}, {jet_jes_shift}, {jet_jer_shift}, {jet_jec_file}, {jet_jer_tag}, {jet_jes_tag}, {jet_jec_algo})",
    input=[
        nanoAOD.Jet_pt,
        nanoAOD.Jet_eta,
        nanoAOD.Jet_phi,
        nanoAOD.Jet_area,
        nanoAOD.Jet_rawFactor,
        nanoAOD.Jet_ID,
        nanoAOD.GenJet_pt,
        nanoAOD.GenJet_eta,
        nanoAOD.GenJet_phi,
        nanoAOD.rho,
    ],
    output=[q.Jet_pt_corrected],
    scopes=["global"],
)

JetMassCorrection = Producer(
    name="JetMassCorrection",
    call="physicsobject::ObjectMassCorrectionWithPt({df}, {output}, {input})",
    input=[
        nanoAOD.Jet_mass,
        nanoAOD.Jet_pt,
        q.Jet_pt_corrected,
    ],
    output=[q.Jet_mass_corrected],
    scopes=["global"],
)
# in data and embdedded sample, we simply rename the nanoAOD jets to the jet_pt_corrected column
RenameJetPt = Producer(
    name="RenameJetPt",
    call="basefunctions::rename<ROOT::RVec<float>>({df}, {input}, {output})",
    input=[nanoAOD.Jet_pt],
    output=[q.Jet_pt_corrected],
    scopes=["global"],
)
RenameJetMass = Producer(
    name="RenameJetMass",
    call="basefunctions::rename<ROOT::RVec<float>>({df}, {input}, {output})",
    input=[nanoAOD.Jet_mass],
    output=[q.Jet_mass_corrected],
    scopes=["global"],
)
RenameJetsData = ProducerGroup(
    name="RenameJetsData",
    call=None,
    input=None,
    output=None,
    scopes=["global"],
    subproducers=[RenameJetPt, RenameJetMass],
)
JetEnergyCorrection = ProducerGroup(
    name="JetEnergyCorrection",
    call=None,
    input=None,
    output=None,
    scopes=["global"],
    subproducers=[JetPtCorrection, JetMassCorrection],
)
'''
JetPtCut = Producer(
    name="JetPtCut",
    call="physicsobject::CutPt({df}, {input}, {output}, {min_jet_pt})",
    input=[nanoAOD.Jet_pt],
    output=[],
    scopes=["global"],
)
BJetPtCut = Producer(
    name="BJetPtCut",
    call="physicsobject::CutPt({df}, {input}, {output}, {min_bjet_pt})",
    input=[nanoAOD.Jet_pt],
    output=[],
    scopes=["global"],
)
JetEtaCut = Producer(
    name="JetEtaCut",
    call="physicsobject::CutEta({df}, {input}, {output}, {max_jet_eta})",
    input=[nanoAOD.Jet_eta],
    output=[],
    scopes=["global"],
)
BJetEtaCut = Producer(
    name="BJetEtaCut",
    call="physicsobject::CutEta({df}, {input}, {output}, {max_bjet_eta})",
    input=[nanoAOD.Jet_eta],
    output=[],
    scopes=["global"],
)
'''
JetIDCut = Producer(
    name="JetIDCut",
    call="physicsobject::jet::CutID({df}, {output}, {input}, {jet_id})",
    input=[nanoAOD.Jet_ID],
    output=[q.jet_id_mask],
    scopes=["global"],
)

JetPUIDCut = Producer(
    name="JetPUIDCut",
    call="physicsobject::jet::CutPUID({df}, {output}, {input}, {jet_puid}, {jet_puid_max_pt})",
    input=[nanoAOD.Jet_PUID, q.Jet_pt_corrected],
    output=[q.jet_puid_mask],
    scopes=["global"],
)
'''
BTagCut = Producer(
    name="BTagCut",
    call="physicsobject::jet::AntiCutRawID({df}, {input}, {output}, {btag_cut})",
    input=[nanoAOD.BJet_discriminator],
    output=[],
    scopes=["global"],
)
GoodJets = ProducerGroup(
    name="GoodJets",
    call="physicsobject::CombineMasks({df}, {output}, {input})",
    input=[],
    output=[q.good_jets_mask],
    scopes=["global"],
    subproducers=[
        #JetPtCut,
        #JetEtaCut, 
        BTagCut
    ],
)
'''
GoodBJets = ProducerGroup(
    name="GoodBJets",
    call="physicsobject::CombineMasks({df}, {output}, {input})",
    input=[q.jet_id_mask, q.jet_puid_mask],
    output=[q.good_bjets_mask],
    scopes=["global"],
    subproducers=[BJetPtCut, BJetEtaCut, BTagCut],
)
'''
####################
# Set of producers to apply a veto of jets overlapping with ditaupair candidates and ordering jets by their pt
# 1. check all jets vs the two lepton candidates, if they are not within deltaR = 0.5, keep them --> mask
# 2. Combine mask with good_jets_mask
# 3. Generate JetCollection, an RVec containing all indices of good Jets in pt order
# 4. generate jet quantity outputs
####################
VetoOverlappingJetsLooseElectrons = Producer(
    name="VetoOverlappingJetsLooseElectrons",
    call="jet::VetoOverlappingJetsLooseLeptons({df}, {output}, {input}, {deltaR_jet_veto})",
    input=[nanoAOD.Jet_eta, nanoAOD.Jet_phi, q.base_electrons_mask, nanoAOD.Electron_eta, nanoAOD.Electron_phi],
    output=[q.jet_electron_overlap_veto_mask],
    scopes=["global"],
)
VetoOverlappingJetsLooseMuons = Producer(
    name="VetoOverlappingJetsLooseMuons",
    call="jet::VetoOverlappingJetsLooseLeptons({df}, {output}, {input}, {deltaR_jet_veto})",
    input=[nanoAOD.Jet_eta, nanoAOD.Jet_phi, q.base_muons_mask, nanoAOD.Muon_eta, nanoAOD.Muon_phi],
    output=[q.jet_muon_overlap_veto_mask],
    scopes=["global"],
)
VetoOverlappingJetsLoosePhotons = Producer(
    name="VetoOverlappingJetsLoosePhotons",
    call="jet::VetoOverlappingJetsLooseLeptons({df}, {output}, {input}, {deltaR_jet_veto})",
    input=[nanoAOD.Jet_eta, nanoAOD.Jet_phi, q.base_photons_mask, nanoAOD.Photon_eta, nanoAOD.Photon_phi],
    output=[q.jet_photon_overlap_veto_mask],
    scopes=["global"],
)

GoodJetsWithVeto = ProducerGroup(
    name="GoodJetsWithVeto",
    call="physicsobject::CombineMasks({df}, {output}, {input})",
    input=[q.good_jets_mask],
    output=[],
    scopes=["global"],
    subproducers=[VetoOverlappingJetsLooseElectrons, VetoOverlappingJetsLooseMuons, VetoOverlappingJetsLoosePhotons],
)
'''
GoodBJetsWithVeto = Producer(
    name="GoodBJetsWithVeto",
    call="physicsobject::CombineMasks({df}, {output}, {input})",
    input=[q.good_bjets_mask, q.jet_overlap_veto_mask],
    output=[],
    scopes=["mt", "et", "tt", "em", "mm", "ee"],
)
'''
JetCollection = ProducerGroup(
    name="JetCollection",
    call="jet::OrderJetsByPt({df}, {output}, {input})",
    input=[nanoAOD.Jet_pt],
    output=[q.good_jet_collection],
    scopes=["global"],
    subproducers=[GoodJetsWithVeto],
)
'''
BJetCollection = ProducerGroup(
    name="BJetCollection",
    call="jet::OrderJetsByPt({df}, {output}, {input})",
    input=[q.Jet_pt_corrected],
    output=[q.good_bjet_collection],
    scopes=["mt", "et", "tt", "em", "mm", "ee"],
    subproducers=[GoodBJetsWithVeto],
)
'''
##########################
# Basic Jet Quantities
# njets, pt, eta, phi, b-tag value
##########################

LVJet1 = Producer(
    name="LVJet1",
    call="lorentzvectors::build({df}, {input_vec}, 0, {output})",
    input=[
        q.good_jet_collection,
        nanoAOD.Jet_pt,
        nanoAOD.Jet_eta,
        nanoAOD.Jet_phi,
        nanoAOD.Jet_mass,
    ],
    output=[q.jet_p4_1],
    scopes=["global"],
)
LVJet2 = Producer(
    name="LVJet2",
    call="lorentzvectors::build({df}, {input_vec}, 1, {output})",
    input=[
        q.good_jet_collection,
        nanoAOD.Jet_pt,
        nanoAOD.Jet_eta,
        nanoAOD.Jet_phi,
        nanoAOD.Jet_mass,
    ],
    output=[q.jet_p4_2],
    scopes=["global"],
)
LVJet3 = Producer(
    name="LVJet3",
    call="lorentzvectors::build({df}, {input_vec}, 2, {output})",
    input=[
        q.good_jet_collection,
        nanoAOD.Jet_pt,
        nanoAOD.Jet_eta,
        nanoAOD.Jet_phi,
        nanoAOD.Jet_mass,
    ],
    output=[q.jet_p4_3],
    scopes=["global"],
)
LVJet4 = Producer(
    name="LVJet4",
    call="lorentzvectors::build({df}, {input_vec}, 3, {output})",
    input=[
        q.good_jet_collection,
        nanoAOD.Jet_pt,
        nanoAOD.Jet_eta,
        nanoAOD.Jet_phi,
        nanoAOD.Jet_mass,
    ],
    output=[q.jet_p4_4],
    scopes=["global"],
)
LVJet5 = Producer(
    name="LVJet5",
    call="lorentzvectors::build({df}, {input_vec}, 4, {output})",
    input=[
        q.good_jet_collection,
        nanoAOD.Jet_pt,
        nanoAOD.Jet_eta,
        nanoAOD.Jet_phi,
        nanoAOD.Jet_mass,
    ],
    output=[q.jet_p4_5],
    scopes=["global"],
)
NumberOfJets = Producer(
    name="NumberOfJets",
    call="quantities::jet::NumberOfJets({df}, {output}, {input})",
    input=[q.good_jet_collection],
    output=[q.njets],
    scopes=["global"],
)
jpt_1 = Producer(
    name="jpt_1",
    call="quantities::pt({df}, {output}, {input})",
    input=[q.jet_p4_1],
    output=[q.jpt_1],
    scopes=["global"],
)
jpt_2 = Producer(
    name="jpt_2",
    call="quantities::pt({df}, {output}, {input})",
    input=[q.jet_p4_2],
    output=[q.jpt_2],
    scopes=["global"],
)
jpt_3 = Producer(
    name="jpt_3",
    call="quantities::pt({df}, {output}, {input})",
    input=[q.jet_p4_3],
    output=[q.jpt_3],
    scopes=["global"],
)
jpt_4 = Producer(
    name="jpt_4",
    call="quantities::pt({df}, {output}, {input})",
    input=[q.jet_p4_4],
    output=[q.jpt_4],
    scopes=["global"],
)
jpt_5 = Producer(
    name="jpt_5",
    call="quantities::pt({df}, {output}, {input})",
    input=[q.jet_p4_5],
    output=[q.jpt_5],
    scopes=["global"],
)
jeta_1 = Producer(
    name="jeta_1",
    call="quantities::eta({df}, {output}, {input})",
    input=[q.jet_p4_1],
    output=[q.jeta_1],
    scopes=["global"],
)
jeta_2 = Producer(
    name="jeta_2",
    call="quantities::eta({df}, {output}, {input})",
    input=[q.jet_p4_2],
    output=[q.jeta_2],
    scopes=["global"],
)
jeta_3 = Producer(
    name="jeta_3",
    call="quantities::eta({df}, {output}, {input})",
    input=[q.jet_p4_3],
    output=[q.jeta_3],
    scopes=["global"],
)
jeta_4 = Producer(
    name="jeta_4",
    call="quantities::eta({df}, {output}, {input})",
    input=[q.jet_p4_4],
    output=[q.jeta_4],
    scopes=["global"],
)
jeta_5 = Producer(
    name="jeta_5",
    call="quantities::eta({df}, {output}, {input})",
    input=[q.jet_p4_5],
    output=[q.jeta_5],
    scopes=["global"],
)
jphi_1 = Producer(
    name="jphi_1",
    call="quantities::phi({df}, {output}, {input})",
    input=[q.jet_p4_1],
    output=[q.jphi_1],
    scopes=["global"],
)
jphi_2 = Producer(
    name="jphi_2",
    call="quantities::phi({df}, {output}, {input})",
    input=[q.jet_p4_2],
    output=[q.jphi_2],
    scopes=["global"],
)
jphi_3 = Producer(
    name="jphi_3",
    call="quantities::phi({df}, {output}, {input})",
    input=[q.jet_p4_3],
    output=[q.jphi_3],
    scopes=["global"],
)
jphi_4 = Producer(
    name="jphi_4",
    call="quantities::phi({df}, {output}, {input})",
    input=[q.jet_p4_4],
    output=[q.jphi_4],
    scopes=["global"],
)
jphi_5 = Producer(
    name="jphi_5",
    call="quantities::phi({df}, {output}, {input})",
    input=[q.jet_p4_5],
    output=[q.jphi_5],
    scopes=["global"],
)
jmass_1 = Producer(
    name="jmass_1",
    call="quantities::mass({df}, {output}, {input})",
    input=[q.jet_p4_1],
    output=[q.jmass_1],
    scopes=["global"],
)
jmass_2 = Producer(
    name="jmass_2",
    call="quantities::mass({df}, {output}, {input})",
    input=[q.jet_p4_2],
    output=[q.jmass_2],
    scopes=["global"],
)
jmass_3 = Producer(
    name="jmass_3",
    call="quantities::mass({df}, {output}, {input})",
    input=[q.jet_p4_3],
    output=[q.jmass_3],
    scopes=["global"],
)
jmass_4 = Producer(
    name="jmass_4",
    call="quantities::mass({df}, {output}, {input})",
    input=[q.jet_p4_4],
    output=[q.jmass_4],
    scopes=["global"],
)
jmass_5 = Producer(
    name="jmass_5",
    call="quantities::mass({df}, {output}, {input})",
    input=[q.jet_p4_5],
    output=[q.jmass_5],
    scopes=["global"],
)
jbtagDeepB_1 = Producer(
    name="jbtagDeepB_1",
    call="quantities::jet::btagValue({df}, {output}, {input}, 0)",
    input=[nanoAOD.Jet_btagDeepB, q.good_jet_collection],
    output=[q.jbtagDeepB_1],
    scopes=["global"],
)
jbtagDeepB_2 = Producer(
    name="jbtagDeepB_2",
    call="quantities::jet::btagValue({df}, {output}, {input}, 1)",
    input=[nanoAOD.Jet_btagDeepB, q.good_jet_collection],
    output=[q.jbtagDeepB_2],
    scopes=["global"],
)
jbtagDeepB_3 = Producer(
    name="jbtagDeepB_3",
    call="quantities::jet::btagValue({df}, {output}, {input}, 2)",
    input=[nanoAOD.Jet_btagDeepB, q.good_jet_collection],
    output=[q.jbtagDeepB_3],
    scopes=["global"],
)
jbtagDeepB_4 = Producer(
    name="jbtagDeepB_4",
    call="quantities::jet::btagValue({df}, {output}, {input}, 3)",
    input=[nanoAOD.Jet_btagDeepB, q.good_jet_collection],
    output=[q.jbtagDeepB_4],
    scopes=["global"],
)
jbtagDeepB_5 = Producer(
    name="jbtagDeepB_5",
    call="quantities::jet::btagValue({df}, {output}, {input}, 4)",
    input=[nanoAOD.Jet_btagDeepB, q.good_jet_collection],
    output=[q.jbtagDeepB_5],
    scopes=["global"],
)
jbtagDeepFlavQG_1 = Producer(
    name="jbtagDeepFlavQG_1",
    call="quantities::jet::btagValue({df}, {output}, {input}, 0)",
    input=[nanoAOD.Jet_btagDeepFlavQG, q.good_jet_collection],
    output=[q.jbtagDeepFlavQG_1],
    scopes=["global"],
)
jbtagDeepFlavQG_2 = Producer(
    name="jbtagDeepFlavQG_2",
    call="quantities::jet::btagValue({df}, {output}, {input}, 1)",
    input=[nanoAOD.Jet_btagDeepFlavQG, q.good_jet_collection],
    output=[q.jbtagDeepFlavQG_2],
    scopes=["global"],
)
jbtagDeepFlavQG_3 = Producer(
    name="jbtagDeepFlavQG_3",
    call="quantities::jet::btagValue({df}, {output}, {input}, 2)",
    input=[nanoAOD.Jet_btagDeepFlavQG, q.good_jet_collection],
    output=[q.jbtagDeepFlavQG_3],
    scopes=["global"],
)
jbtagDeepFlavQG_4 = Producer(
    name="jbtagDeepFlavQG_4",
    call="quantities::jet::btagValue({df}, {output}, {input}, 3)",
    input=[nanoAOD.Jet_btagDeepFlavQG, q.good_jet_collection],
    output=[q.jbtagDeepFlavQG_4],
    scopes=["global"],
)
jbtagDeepFlavQG_5 = Producer(
    name="jbtagDeepFlavQG_5",
    call="quantities::jet::btagValue({df}, {output}, {input}, 4)",
    input=[nanoAOD.Jet_btagDeepFlavQG, q.good_jet_collection],
    output=[q.jbtagDeepFlavQG_5],
    scopes=["global"],
)
jchEmEF_1 = Producer(
    name="jchEmEF_1",
    call="quantities::jet::btagValue({df}, {output}, {input}, 0)",
    input=[nanoAOD.Jet_chEmEF, q.good_jet_collection],
    output=[q.jchEmEF_1],
    scopes=["global"],
)
jchEmEF_2 = Producer(
    name="jchEmEF_2",
    call="quantities::jet::btagValue({df}, {output}, {input}, 1)",
    input=[nanoAOD.Jet_chEmEF, q.good_jet_collection],
    output=[q.jchEmEF_2],
    scopes=["global"],
)
jchEmEF_3 = Producer(
    name="jchEmEF_3",
    call="quantities::jet::btagValue({df}, {output}, {input}, 2)",
    input=[nanoAOD.Jet_chEmEF, q.good_jet_collection],
    output=[q.jchEmEF_3],
    scopes=["global"],
)
jchEmEF_4 = Producer(
    name="jchEmEF_4",
    call="quantities::jet::btagValue({df}, {output}, {input}, 3)",
    input=[nanoAOD.Jet_chEmEF, q.good_jet_collection],
    output=[q.jchEmEF_4],
    scopes=["global"],
)
jchEmEF_5 = Producer(
    name="jchEmEF_5",
    call="quantities::jet::btagValue({df}, {output}, {input}, 4)",
    input=[nanoAOD.Jet_chEmEF, q.good_jet_collection],
    output=[q.jchEmEF_5],
    scopes=["global"],
)
jchHEF_1 = Producer(
    name="jchHEF_1",
    call="quantities::jet::btagValue({df}, {output}, {input}, 0)",
    input=[nanoAOD.Jet_chHEF, q.good_jet_collection],
    output=[q.jchHEF_1],
    scopes=["global"],
)
jchHEF_2 = Producer(
    name="jchHEF_2",
    call="quantities::jet::btagValue({df}, {output}, {input}, 1)",
    input=[nanoAOD.Jet_chHEF, q.good_jet_collection],
    output=[q.jchHEF_2],
    scopes=["global"],
)
jchHEF_3 = Producer(
    name="jchHEF_3",
    call="quantities::jet::btagValue({df}, {output}, {input}, 2)",
    input=[nanoAOD.Jet_chHEF, q.good_jet_collection],
    output=[q.jchHEF_3],
    scopes=["global"],
)
jchHEF_4 = Producer(
    name="jchHEF_4",
    call="quantities::jet::btagValue({df}, {output}, {input}, 3)",
    input=[nanoAOD.Jet_chHEF, q.good_jet_collection],
    output=[q.jchHEF_4],
    scopes=["global"],
)
jchHEF_5 = Producer(
    name="jchHEF_5",
    call="quantities::jet::btagValue({df}, {output}, {input}, 4)",
    input=[nanoAOD.Jet_chHEF, q.good_jet_collection],
    output=[q.jchHEF_5],
    scopes=["global"],
)
jneEmEF_1 = Producer(
    name="jneEmEF_1",
    call="quantities::jet::btagValue({df}, {output}, {input}, 0)",
    input=[nanoAOD.Jet_neEmEF, q.good_jet_collection],
    output=[q.jneEmEF_1],
    scopes=["global"],
)
jneEmEF_2 = Producer(
    name="jneEmEF_2",
    call="quantities::jet::btagValue({df}, {output}, {input}, 1)",
    input=[nanoAOD.Jet_neEmEF, q.good_jet_collection],
    output=[q.jneEmEF_2],
    scopes=["global"],
)
jneEmEF_3 = Producer(
    name="jneEmEF_3",
    call="quantities::jet::btagValue({df}, {output}, {input}, 2)",
    input=[nanoAOD.Jet_neEmEF, q.good_jet_collection],
    output=[q.jneEmEF_3],
    scopes=["global"],
)
jneEmEF_4 = Producer(
    name="jneEmEF_4",
    call="quantities::jet::btagValue({df}, {output}, {input}, 3)",
    input=[nanoAOD.Jet_neEmEF, q.good_jet_collection],
    output=[q.jneEmEF_4],
    scopes=["global"],
)
jneEmEF_5 = Producer(
    name="jneEmEF_5",
    call="quantities::jet::btagValue({df}, {output}, {input}, 4)",
    input=[nanoAOD.Jet_neEmEF, q.good_jet_collection],
    output=[q.jneEmEF_5],
    scopes=["global"],
)
jneHEF_1 = Producer(
    name="jneHEF_1",
    call="quantities::jet::btagValue({df}, {output}, {input}, 0)",
    input=[nanoAOD.Jet_neHEF, q.good_jet_collection],
    output=[q.jneHEF_1],
    scopes=["global"],
)
jneHEF_2 = Producer(
    name="jneHEF_2",
    call="quantities::jet::btagValue({df}, {output}, {input}, 1)",
    input=[nanoAOD.Jet_neHEF, q.good_jet_collection],
    output=[q.jneHEF_2],
    scopes=["global"],
)
jneHEF_3 = Producer(
    name="jneHEF_3",
    call="quantities::jet::btagValue({df}, {output}, {input}, 2)",
    input=[nanoAOD.Jet_neHEF, q.good_jet_collection],
    output=[q.jneHEF_3],
    scopes=["global"],
)
jneHEF_4 = Producer(
    name="jneHEF_4",
    call="quantities::jet::btagValue({df}, {output}, {input}, 3)",
    input=[nanoAOD.Jet_neHEF, q.good_jet_collection],
    output=[q.jneHEF_4],
    scopes=["global"],
)
jneHEF_5 = Producer(
    name="jneHEF_5",
    call="quantities::jet::btagValue({df}, {output}, {input}, 4)",
    input=[nanoAOD.Jet_neHEF, q.good_jet_collection],
    output=[q.jneHEF_5],
    scopes=["global"],
)
mjj = Producer(
    name="jphi_2",
    call="quantities::m_vis({df}, {output}, {input_vec})",
    input=[q.jet_p4_1, q.jet_p4_2],
    output=[q.mjj],
    scopes=["mt", "et", "tt", "em", "mm", "ee"],
)
BasicJetQuantities = ProducerGroup(
    name="BasicJetQuantities",
    call=None,
    input=None,
    output=None,
    scopes=["global"],
    subproducers=[
        LVJet1,
        LVJet2,
        LVJet3,
        LVJet4,
        LVJet5,
        NumberOfJets,
        jpt_1,
        jeta_1,
        jphi_1,
        jmass_1,
        jbtagDeepB_1,
        jbtagDeepFlavQG_1,
        jchEmEF_1,
        jchHEF_1,
        jneEmEF_1,
        jneHEF_1,
        jpt_2,
        jeta_2,
        jphi_2,
        jmass_2,
        jbtagDeepB_2,
        jbtagDeepFlavQG_2,
        jchEmEF_2,
        jchHEF_2,
        jneEmEF_2,
        jneHEF_2,
        jpt_3,
        jeta_3,
        jphi_3,
        jmass_3,
        jbtagDeepB_3,
        jbtagDeepFlavQG_3,
        jchEmEF_3,
        jchHEF_3,
        jneEmEF_3,
        jneHEF_3,
        jpt_4,
        jeta_4,
        jphi_4,
        jmass_4,
        jbtagDeepB_4,
        jbtagDeepFlavQG_4,
        jchEmEF_4,
        jchHEF_4,
        jneEmEF_4,
        jneHEF_4,
        jpt_5,
        jeta_5,
        jphi_5,
        jmass_5,
        jbtagDeepB_5,
        jbtagDeepFlavQG_5,
        jchEmEF_5,
        jchHEF_5,
        jneEmEF_5,
        jneHEF_5,
    ],
)

##########################
# Basic b-Jet Quantities
# nbtag, pt, eta, phi, b-tag value
##########################
'''
LVBJet1 = Producer(
    name="LVBJet1",
    call="lorentzvectors::build({df}, {input_vec}, 0, {output})",
    input=[
        q.good_bjet_collection,
        q.Jet_pt_corrected,
        nanoAOD.Jet_eta,
        nanoAOD.Jet_phi,
        q.Jet_mass_corrected,
    ],
    output=[q.bjet_p4_1],
    scopes=["mt", "et", "tt", "em", "mm", "ee"],
)
LVBJet2 = Producer(
    name="LVBJet2",
    call="lorentzvectors::build({df}, {input_vec}, 1, {output})",
    input=[
        q.good_bjet_collection,
        q.Jet_pt_corrected,
        nanoAOD.Jet_eta,
        nanoAOD.Jet_phi,
        q.Jet_mass_corrected,
    ],
    output=[q.bjet_p4_2],
    scopes=["mt", "et", "tt", "em", "mm", "ee"],
)
NumberOfBJets = Producer(
    name="NumberOfBJets",
    call="quantities::jet::NumberOfJets({df}, {output}, {input})",
    input=[q.good_bjet_collection],
    output=[q.nbtag],
    scopes=["mt", "et", "tt", "em", "mm", "ee"],
)
bpt_1 = Producer(
    name="bpt_1",
    call="quantities::pt({df}, {output}, {input})",
    input=[q.bjet_p4_1],
    output=[q.bpt_1],
    scopes=["mt", "et", "tt", "em", "mm", "ee"],
)
bpt_2 = Producer(
    name="bpt_2",
    call="quantities::pt({df}, {output}, {input})",
    input=[q.bjet_p4_2],
    output=[q.bpt_2],
    scopes=["mt", "et", "tt", "em", "mm", "ee"],
)
beta_1 = Producer(
    name="beta_1",
    call="quantities::eta({df}, {output}, {input})",
    input=[q.bjet_p4_1],
    output=[q.beta_1],
    scopes=["mt", "et", "tt", "em", "mm", "ee"],
)
beta_2 = Producer(
    name="beta_2",
    call="quantities::eta({df}, {output}, {input})",
    input=[q.bjet_p4_2],
    output=[q.beta_2],
    scopes=["mt", "et", "tt", "em", "mm", "ee"],
)
bphi_1 = Producer(
    name="bphi_1",
    call="quantities::phi({df}, {output}, {input})",
    input=[q.bjet_p4_1],
    output=[q.bphi_1],
    scopes=["mt", "et", "tt", "em", "mm", "ee"],
)
bphi_2 = Producer(
    name="bphi_2",
    call="quantities::phi({df}, {output}, {input})",
    input=[q.bjet_p4_2],
    output=[q.bphi_2],
    scopes=["mt", "et", "tt", "em", "mm", "ee"],
)
btag_value_1 = Producer(
    name="btag_value_1",
    call="quantities::jet::btagValue({df}, {output}, {input}, 0)",
    input=[nanoAOD.BJet_discriminator, q.good_bjet_collection],
    output=[q.btag_value_1],
    scopes=["mt", "et", "tt", "em", "mm", "ee"],
)
btag_value_2 = Producer(
    name="btag_value_2",
    call="quantities::jet::btagValue({df}, {output}, {input}, 1)",
    input=[nanoAOD.BJet_discriminator, q.good_bjet_collection],
    output=[q.btag_value_2],
    scopes=["mt", "et", "tt", "em", "mm", "ee"],
)
BasicBJetQuantities = ProducerGroup(
    name="BasicBJetQuantities",
    call=None,
    input=None,
    output=None,
    scopes=["mt", "et", "tt", "em", "mm", "ee"],
    subproducers=[
        LVBJet1,
        LVBJet2,
        NumberOfBJets,
        bpt_1,
        beta_1,
        bphi_1,
        btag_value_1,
        bpt_2,
        beta_2,
        bphi_2,
        btag_value_2,
    ],
)
'''
