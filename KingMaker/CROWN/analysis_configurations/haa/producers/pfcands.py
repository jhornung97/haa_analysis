from ..quantities import output as q
from ..quantities import nanoAOD as nanoAOD
from code_generation.producer import Producer, ProducerGroup, Filter

pfCandpdgId = Producer(
    name="PFCandsPdgId",
    call="basefunctions::rename<ROOT::RVec<Int_t>>({df}, {input}, {output})",
    input=[nanoAOD.PFCands_pdgId],
    output=[q.pfcands_pdgid],
    scopes=["global"],
)

pfCandpt = Producer(
    name="PFCandspt",
    call="basefunctions::rename<ROOT::RVec<Float_t>>({df}, {input}, {output})",
    input=[nanoAOD.PFCands_pt],
    output=[q.pfcands_pt],
    scopes=["global"],
)

pfCandeta = Producer(
    name="PFCandseta",
    call="basefunctions::rename<ROOT::RVec<Float_t>>({df}, {input}, {output})",
    input=[nanoAOD.PFCands_eta],
    output=[q.pfcands_eta],
    scopes=["global"],
)

pfCandphi = Producer(
    name="PFCandsphi",
    call="basefunctions::rename<ROOT::RVec<Float_t>>({df}, {input}, {output})",
    input=[nanoAOD.PFCands_phi],
    output=[q.pfcands_phi],
    scopes=["global"],
)

pfCandmass = Producer(
    name="PFCandsmass",
    call="basefunctions::rename<ROOT::RVec<Float_t>>({df}, {input}, {output})",
    input=[nanoAOD.PFCands_mass],
    output=[q.pfcands_mass],
    scopes=["global"],
)

####################                                                                                                                                                                                               
# Set of producers used for loosest selection of pfcandss                                                                                                                                                             
####################                                                       
ChargedPFCandsPdgId = Producer(
    name="PFCandsPdgId",
    call="physicsobject::tau::CutDecayModes({df}, {output}, {input}, {vec_open}{pfcands_pdgid}{vec_close})",
    input=[nanoAOD.PFCands_pdgId],
    output=[],
    scopes=["global"],
)                                                                                                                                
PFCandsPtCut = Producer(
    name="PFCandsPtCut",
    call="physicsobject::CutPt({df}, {input}, {output}, {min_pfcands_pt})",
    input=[nanoAOD.PFCands_pt],
    output=[],
    scopes=["global"],
)
PFCandsEtaCut = Producer(
    name="PFCandsEtaCut",
    call="physicsobject::CutEta({df}, {input}, {output}, {max_pfcands_eta})",
    input=[nanoAOD.PFCands_eta],
    output=[],
    scopes=["global"],
)

PFCandsMinMassCut = Producer(
    name="PFCandsMassCut",
    call="physicsobject::CutPt({df}, {input}, {output}, {min_pfcands_mass})",
    input=[nanoAOD.PFCands_mass],
    output=[],
    scopes=["global"],
)

PFCandsMaxMassCut = Producer(
    name="PFCandsMassCut",
    call="physicsobject::CutEta({df}, {input}, {output}, {max_pfcands_mass})",
    input=[nanoAOD.PFCands_mass],
    output=[],
    scopes=["global"],
)

PFCandsFromPVSelection = Producer(
    name="PFCandsFromPVSelection",
    call="physicsobject::tau::CutDecayModes({df}, {output}, {input}, {vec_open}{pfcands_fromPV_flag}{vec_close})",
    input=[nanoAOD.PFCands_fromPV],
    output=[],
    scopes=["global"],
)

VetoOverlappingPFCandsLooseElectrons = Producer(
    name="VetoOverlappingPFCandsLooseElectrons",
    call="jet::VetoOverlappingJetsLooseLeptons({df}, {output}, {input}, {deltaR_jet_veto})",
    input=[nanoAOD.PFCands_eta, nanoAOD.PFCands_phi, q.base_electrons_mask, nanoAOD.Electron_eta, nanoAOD.Electron_phi],
    output=[q.pfcand_electron_overlap_veto_mask],
    scopes=["global"],
)
VetoOverlappingPFCandsLooseMuons = Producer(
    name="VetoOverlappingPFCandsLooseMuons",
    call="jet::VetoOverlappingJetsLooseLeptons({df}, {output}, {input}, {deltaR_jet_veto})",
    input=[nanoAOD.PFCands_eta, nanoAOD.PFCands_phi, q.base_muons_mask, nanoAOD.Muon_eta, nanoAOD.Muon_phi],
    output=[q.pfcand_muon_overlap_veto_mask],
    scopes=["global"],
)
VetoOverlappingPFCandsLoosePhotons = Producer(
    name="VetoOverlappingPFCandsLoosePhotons",
    call="jet::VetoOverlappingJetsLooseLeptons({df}, {output}, {input}, {deltaR_jet_veto})",
    input=[nanoAOD.PFCands_eta, nanoAOD.PFCands_phi, q.base_photons_mask, nanoAOD.Photon_eta, nanoAOD.Photon_phi],
    output=[q.pfcand_photon_overlap_veto_mask],
    scopes=["global"],
)
BasePFCands = ProducerGroup(
    name="BasePFCands",
    call="physicsobject::CombineMasks({df}, {output}, {input})",
    input=[],
    output=[q.base_pfcands_mask],
    scopes=["global"],
    subproducers=[
        ChargedPFCandsPdgId,
        VetoOverlappingPFCandsLooseElectrons,
        VetoOverlappingPFCandsLooseMuons,
        VetoOverlappingPFCandsLoosePhotons,
        PFCandsPtCut,
        PFCandsEtaCut,
        PFCandsMinMassCut,
        PFCandsMaxMassCut,
        #PFCandsFromPVSelection,
    ],
)
ChargedPFCands = Producer(
    name="ChargedPFCands",
    call="physicsobject::tau::CutDecayModes({df}, {output}, {input}, {vec_open}{charged_pfcands_pdgid}{vec_close})",
    input=[nanoAOD.PFCands_pdgId],
    output=[q.ChargedPFCands],
    scopes=["global"],
)
NeutralPFCands = Producer(
    name="NeutralPFCands",
    call="physicsobject::tau::CutDecayModes({df}, {output}, {input}, {vec_open}{neutral_pfcands_pdgid}{vec_close})",
    input=[nanoAOD.PFCands_pdgId],
    output=[q.NeutralPFCands],
    scopes=["global"],
)
PhotonPFCands = Producer(
    name="PhotonPFCands",
    call="physicsobject::tau::CutDecayModes({df}, {output}, {input}, {vec_open}{photon_pfcands_pdgid}{vec_close})",
    input=[nanoAOD.PFCands_pdgId],
    output=[q.PhotonPFCands],
    scopes=["global"],
)
fromPV = Producer(
    name="PFCands_fromPV",
    call="physicsobject::tau::CutDecayModes({df}, {output}, {input}, {vec_open}{fromPV}{vec_close})",
    input=[nanoAOD.PFCands_fromPV],
    output=[q.fromPV],
    scopes=["global"],
)
GenPartPdgId = Producer(
    name="GenPartPdgId",
    call="physicsobject::tau::CutDecayModes({df}, {output}, {input}, {vec_open}{genpart_pdgid}{vec_close})",
    input=[nanoAOD.GenParticle_pdgId],
    output=[],
    scopes=["global"],
)
GenPartStatus = Producer(
    name="GenPartStatus",
    call="physicsobject::tau::CutDecayModes({df}, {output}, {input}, {vec_open}{genpart_status}{vec_close})",
    input=[nanoAOD.GenParticle_status],
    output=[],
    scopes=["global"],
)
BaseZs = ProducerGroup(
    name="BaseZs",
    call="physicsobject::CombineMasks({df}, {output}, {input})",
    input=[],
    output=[q.base_zs_mask],
    scopes=["global"],
    subproducers=[
        GenPartPdgId,
        GenPartStatus
    ],
)

z_pt = Producer(
    name="z_pt",
    call="haa::getGenPt({df}, {input}, {output})",
    input=[q.base_zs_mask, nanoAOD.GenParticle_pt],
    output=[q.z_pt],
    scopes=["mm"],
)

ZFlag = Producer(
    name="ZFlag",
    call="physicsobject::CutNFlag({df}, {output}, {input}, 1)",
    input=[q.base_zs_mask],
    output=[],
    scopes=["global"],
)

ZFilter = Filter(
    name="ZFilter",
    call='basefunctions::FilterFlagsAny({df}, "Zs", {input})',
    input=[],
    scopes=["global"],
    subproducers=[ZFlag],
)

GetGenPartPdgId = Producer(
    name="GetGenPartPdgId",
    call="basefunctions::rename<ROOT::RVec<Int_t>>({df}, {input}, {output})",
    input=[nanoAOD.GenParticle_pdgId],
    output=[q.genpart_pdgid],
    scopes=["global"],
)

GetGenPartStatus = Producer(
    name="GetGenPartStatus",
    call="basefunctions::rename<ROOT::RVec<Int_t>>({df}, {input}, {output})",
    input=[nanoAOD.GenParticle_status],
    output=[q.genpart_status],
    scopes=["global"],
)

GetGenPartEta = Producer(
    name="GetGenPartEta",
    call="basefunctions::rename<ROOT::RVec<Float_t>>({df}, {input}, {output})",
    input=[nanoAOD.GenParticle_eta],
    output=[q.genpart_eta],
    scopes=["global"],
)

GetGenPartPhi = Producer(
    name="GetGenPartPhi",
    call="basefunctions::rename<ROOT::RVec<Float_t>>({df}, {input}, {output})",
    input=[nanoAOD.GenParticle_phi],
    output=[q.genpart_phi],
    scopes=["global"],
)