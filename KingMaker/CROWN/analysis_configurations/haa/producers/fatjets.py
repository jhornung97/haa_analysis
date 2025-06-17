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
FatJetPtCut = Producer(
    name="FatJetPtCut",
    call="physicsobject::CutPt({df}, {input}, {output}, {min_fatjet_pt})",
    input=[nanoAOD.FatJet_pt],
    output=[],
    scopes=["global"],
)
'''
BJetPtCut = Producer(
    name="BJetPtCut",
    call="physicsobject::CutPt({df}, {input}, {output}, {min_bjet_pt})",
    input=[q.Jet_pt_corrected],
    output=[],
    scopes=["global"],
)
JetEtaCut = Producer(
    name="JetEtaCut",
    call="physicsobject::CutEta({df}, {input}, {output}, {max_jet_eta})",
    input=[nanoAOD.FatJet_eta],
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
BTagCut = Producer(
    name="BTagCut",
    call="physicsobject::jet::CutRawID({df}, {input}, {output}, {btag_cut})",
    input=[nanoAOD.BJet_discriminator],
    output=[],
    scopes=["global"],
)
'''
GoodFatJets = ProducerGroup(
    name="GoodJets",
    call="physicsobject::CombineMasks({df}, {output}, {input})",
    input=[],
    output=[q.good_fatjets_mask],
    scopes=["global"],
    subproducers=[FatJetPtCut],
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

####################
# Set of producers to apply a veto of jets overlapping with ditaupair candidates and ordering jets by their pt
# 1. check all jets vs the two lepton candidates, if they are not within deltaR = 0.5, keep them --> mask
# 2. Combine mask with good_jets_mask
# 3. Generate JetCollection, an RVec containing all indices of good Jets in pt order
# 4. generate jet quantity outputs
####################
VetoOverlappingJets = Producer(
    name="VetoOverlappingJets",
    call="jet::VetoOverlappingJets({df}, {output}, {input}, {deltaR_jet_veto})",
    input=[nanoAOD.Jet_eta, nanoAOD.Jet_phi, q.p4_1, q.p4_2],
    output=[q.jet_overlap_veto_mask],
    scopes=["mt", "et", "tt", "em", "mm", "ee"],
)

GoodJetsWithVeto = ProducerGroup(
    name="GoodJetsWithVeto",
    call="physicsobject::CombineMasks({df}, {output}, {input})",
    input=[q.good_jets_mask],
    output=[],
    scopes=["mt", "et", "tt", "em", "mm", "ee"],
    subproducers=[VetoOverlappingJets],
)

GoodBJetsWithVeto = Producer(
    name="GoodBJetsWithVeto",
    call="physicsobject::CombineMasks({df}, {output}, {input})",
    input=[q.good_bjets_mask, q.jet_overlap_veto_mask],
    output=[],
    scopes=["mt", "et", "tt", "em", "mm", "ee"],
)
'''
FatJetCollection = ProducerGroup(
    name="FatJetCollection",
    call="jet::OrderJetsByPt({df}, {output}, {input})",
    input=[nanoAOD.FatJet_pt],
    output=[q.good_fatjet_collection],
    scopes=["global"],
    subproducers=[GoodFatJets],
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

##########################
# Basic Jet Quantities
# njets, pt, eta, phi, b-tag value
##########################
'''
LVFatJet1 = Producer(
    name="LVFatFatjet1",
    call="lorentzvectors::build({df}, {input_vec}, 0, {output})",
    input=[
        q.good_fatjet_collection,
        nanoAOD.FatJet_pt,
        nanoAOD.FatJet_eta,
        nanoAOD.FatJet_phi,
        nanoAOD.FatJet_mass,
    ],
    output=[q.fatjet_p4_1],
    scopes=["global"],
)
LVFatJet2 = Producer(
    name="LVFatJet2",
    call="lorentzvectors::build({df}, {input_vec}, 1, {output})",
    input=[
        q.good_fatjet_collection,
        nanoAOD.FatJet_pt,
        nanoAOD.FatJet_eta,
        nanoAOD.FatJet_phi,
        nanoAOD.FatJet_mass,
    ],
    output=[q.fatjet_p4_2],
    scopes=["global"],
)
LVFatJet3 = Producer(
    name="LVFatJet3",
    call="lorentzvectors::build({df}, {input_vec}, 2, {output})",
    input=[
        q.good_fatjet_collection,
        nanoAOD.FatJet_pt,
        nanoAOD.FatJet_eta,
        nanoAOD.FatJet_phi,
        nanoAOD.FatJet_mass,
    ],
    output=[q.fatjet_p4_3],
    scopes=["global"],
)
NumberOfFatJets = Producer(
    name="NumberOfFatJets",
    call="quantities::jet::NumberOfJets({df}, {output}, {input})",
    input=[q.good_fatjet_collection],
    output=[q.nfatjets],
    scopes=["global"],
)
fjpt_1 = Producer(
    name="fjpt_1",
    call="quantities::pt({df}, {output}, {input})",
    input=[q.fatjet_p4_1],
    output=[q.fjpt_1],
    scopes=["global"],
)
fjpt_2 = Producer(
    name="fjpt_2",
    call="quantities::pt({df}, {output}, {input})",
    input=[q.fatjet_p4_2],
    output=[q.fjpt_2],
    scopes=["global"],
)
fjpt_3 = Producer(
    name="fjpt_3",
    call="quantities::pt({df}, {output}, {input})",
    input=[q.fatjet_p4_3],
    output=[q.fjpt_3],
    scopes=["global"],
)
fjeta_1 = Producer(
    name="fjeta_1",
    call="quantities::eta({df}, {output}, {input})",
    input=[q.fatjet_p4_1],
    output=[q.fjeta_1],
    scopes=["global"],
)
fjeta_2 = Producer(
    name="fjeta_2",
    call="quantities::eta({df}, {output}, {input})",
    input=[q.fatjet_p4_2],
    output=[q.fjeta_2],
    scopes=["global"],
)
fjeta_3 = Producer(
    name="fjeta_3",
    call="quantities::eta({df}, {output}, {input})",
    input=[q.fatjet_p4_3],
    output=[q.fjeta_3],
    scopes=["global"],
)
fjphi_1 = Producer(
    name="fjphi_1",
    call="quantities::phi({df}, {output}, {input})",
    input=[q.fatjet_p4_1],
    output=[q.fjphi_1],
    scopes=["global"],
)
fjphi_2 = Producer(
    name="fjphi_2",
    call="quantities::phi({df}, {output}, {input})",
    input=[q.fatjet_p4_2],
    output=[q.fjphi_2],
    scopes=["global"],
)
fjphi_3 = Producer(
    name="fjphi_3",
    call="quantities::phi({df}, {output}, {input})",
    input=[q.fatjet_p4_3],
    output=[q.fjphi_3],
    scopes=["global"],
)
fjmass_1 = Producer(
    name="fjmass_1",
    call="quantities::mass({df}, {output}, {input})",
    input=[q.fatjet_p4_1],
    output=[q.fjmass_1],
    scopes=["global"],
)
fjmass_2 = Producer(
    name="fjmass_2",
    call="quantities::mass({df}, {output}, {input})",
    input=[q.fatjet_p4_2],
    output=[q.fjmass_2],
    scopes=["global"],
)
fjmass_3 = Producer(
    name="fjmass_3",
    call="quantities::mass({df}, {output}, {input})",
    input=[q.fatjet_p4_3],
    output=[q.fjmass_3],
    scopes=["global"],
)
fjmsoftdrop_1 = Producer(
    name="fjmsoftdrop_1",
    call="quantities::jet::btagValue({df}, {output}, {input}, 0)",
    input=[nanoAOD.FatJet_msoftdrop, q.good_fatjet_collection],
    output=[q.fjmsoftdrop_1],
    scopes=["global"],
)
fjmsoftdrop_2 = Producer(
    name="fjmsoftdrop_2",
    call="quantities::jet::btagValue({df}, {output}, {input}, 1)",
    input=[nanoAOD.FatJet_msoftdrop, q.good_fatjet_collection],
    output=[q.fjmsoftdrop_2],
    scopes=["global"],
)
fjmsoftdrop_3 = Producer(
    name="fjmsoftdrop_3",
    call="quantities::jet::btagValue({df}, {output}, {input}, 2)",
    input=[nanoAOD.FatJet_msoftdrop, q.good_fatjet_collection],
    output=[q.fjmsoftdrop_3],
    scopes=["global"],
)
fjn2b1_1 = Producer(
    name="fjn2b1_1",
    call="quantities::jet::btagValue({df}, {output}, {input}, 0)",
    input=[nanoAOD.FatJet_n2b1, q.good_fatjet_collection],
    output=[q.fjn2b1_1],
    scopes=["global"],
)
fjn2b1_2 = Producer(
    name="fjn2b1_2",
    call="quantities::jet::btagValue({df}, {output}, {input}, 1)",
    input=[nanoAOD.FatJet_n2b1, q.good_fatjet_collection],
    output=[q.fjn2b1_2],
    scopes=["global"],
)
fjn2b1_3 = Producer(
    name="fjn2b1_3",
    call="quantities::jet::btagValue({df}, {output}, {input}, 2)",
    input=[nanoAOD.FatJet_n2b1, q.good_fatjet_collection],
    output=[q.fjn2b1_3],
    scopes=["global"],
)
fjn3b1_1 = Producer(
    name="fjn3b1_1",
    call="quantities::jet::btagValue({df}, {output}, {input}, 0)",
    input=[nanoAOD.FatJet_n3b1, q.good_fatjet_collection],
    output=[q.fjn3b1_1],
    scopes=["global"],
)
fjn3b1_2 = Producer(
    name="fjn3b1_2",
    call="quantities::jet::btagValue({df}, {output}, {input}, 1)",
    input=[nanoAOD.FatJet_n3b1, q.good_fatjet_collection],
    output=[q.fjn3b1_2],
    scopes=["global"],
)
fjn3b1_3 = Producer(
    name="fjn3b1_3",
    call="quantities::jet::btagValue({df}, {output}, {input}, 2)",
    input=[nanoAOD.FatJet_n3b1, q.good_fatjet_collection],
    output=[q.fjn3b1_3],
    scopes=["global"],
)
fjparticleNet_H4qvsQCD_1 = Producer(
    name="fjparticleNet_H4qvsQCD_1",
    call="quantities::jet::btagValue({df}, {output}, {input}, 0)",
    input=[nanoAOD.FatJet_particleNet_H4qvsQCD, q.good_fatjet_collection],
    output=[q.fjparticleNet_H4qvsQCD_1],
    scopes=["global"],
)
fjparticleNet_H4qvsQCD_2 = Producer(
    name="fjparticleNet_H4qvsQCD_2",
    call="quantities::jet::btagValue({df}, {output}, {input}, 1)",
    input=[nanoAOD.FatJet_particleNet_H4qvsQCD, q.good_fatjet_collection],
    output=[q.fjparticleNet_H4qvsQCD_2],
    scopes=["global"],
)
fjparticleNet_H4qvsQCD_3 = Producer(
    name="fjparticleNet_H4qvsQCD_3",
    call="quantities::jet::btagValue({df}, {output}, {input}, 2)",
    input=[nanoAOD.FatJet_particleNet_H4qvsQCD, q.good_fatjet_collection],
    output=[q.fjparticleNet_H4qvsQCD_3],
    scopes=["global"],
)
fjparticleNet_QCD_1 = Producer(
    name="fjparticleNet_QCD_1",
    call="quantities::jet::btagValue({df}, {output}, {input}, 0)",
    input=[nanoAOD.FatJet_particleNet_QCD, q.good_fatjet_collection],
    output=[q.fjparticleNet_QCD_1],
    scopes=["global"],
)
fjparticleNet_QCD_2 = Producer(
    name="fjparticleNet_QCD_2",
    call="quantities::jet::btagValue({df}, {output}, {input}, 1)",
    input=[nanoAOD.FatJet_particleNet_QCD, q.good_fatjet_collection],
    output=[q.fjparticleNet_QCD_2],
    scopes=["global"],
)
fjparticleNet_QCD_3 = Producer(
    name="fjparticleNet_QCD_3",
    call="quantities::jet::btagValue({df}, {output}, {input}, 2)",
    input=[nanoAOD.FatJet_particleNet_QCD, q.good_fatjet_collection],
    output=[q.fjparticleNet_QCD_3],
    scopes=["global"],
)
fjparticleNet_mass_1 = Producer(
    name="fjparticleNet_mass_1",
    call="quantities::jet::btagValue({df}, {output}, {input}, 0)",
    input=[nanoAOD.FatJet_particleNet_mass, q.good_fatjet_collection],
    output=[q.fjparticleNet_mass_1],
    scopes=["global"],
)
fjparticleNet_mass_2 = Producer(
    name="fjparticleNet_mass_2",
    call="quantities::jet::btagValue({df}, {output}, {input}, 1)",
    input=[nanoAOD.FatJet_particleNet_mass, q.good_fatjet_collection],
    output=[q.fjparticleNet_mass_2],
    scopes=["global"],
)
fjparticleNet_mass_3 = Producer(
    name="fjparticleNet_mass_3",
    call="quantities::jet::btagValue({df}, {output}, {input}, 2)",
    input=[nanoAOD.FatJet_particleNet_mass, q.good_fatjet_collection],
    output=[q.fjparticleNet_mass_3],
    scopes=["global"],
)
mjj = Producer(
    name="jphi_2",
    call="quantities::m_vis({df}, {output}, {input_vec})",
    input=[q.jet_p4_1, q.jet_p4_2],
    output=[q.mjj],
    scopes=["mt", "et", "tt", "em", "mm", "ee"],
)
BasicFatJetQuantities = ProducerGroup(
    name="BasicFatJetQuantities",
    call=None,
    input=None,
    output=None,
    scopes=["global"],
    subproducers=[
        LVFatJet1,
        LVFatJet2,
        LVFatJet3,
        NumberOfFatJets,
        fjpt_1,
        fjeta_1,
        fjphi_1,
        fjmass_1,
        fjmsoftdrop_1,
        fjn2b1_1,
        fjn3b1_1,
        fjparticleNet_H4qvsQCD_1,
        fjparticleNet_QCD_1,
        fjparticleNet_mass_1,
        fjpt_2,
        fjeta_2,
        fjphi_2,
        fjmass_2,
        fjmsoftdrop_2,
        fjn2b1_2,
        fjn3b1_2,
        fjparticleNet_H4qvsQCD_2,
        fjparticleNet_QCD_2,
        fjparticleNet_mass_2,
        fjpt_3,
        fjeta_3,
        fjphi_3,
        fjmass_3,
        fjmsoftdrop_3,
        fjn2b1_3,
        fjn3b1_3,
        fjparticleNet_H4qvsQCD_3,
        fjparticleNet_QCD_3,
        fjparticleNet_mass_3,
    ],
)
'''
##########################
# Basic b-Jet Quantities
# nbtag, pt, eta, phi, b-tag value
##########################

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
    call="quantities:::jet:::NumberOfJets({df}, {output}, {input})",
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
    call="quantities:::jet:::btagValue({df}, {output}, {input}, 0)",
    input=[nanoAOD.BJet_discriminator, q.good_bjet_collection],
    output=[q.btag_value_1],
    scopes=["mt", "et", "tt", "em", "mm", "ee"],
)
btag_value_2 = Producer(
    name="btag_value_2",
    call="quantities:::jet:::btagValue({df}, {output}, {input}, 1)",
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
