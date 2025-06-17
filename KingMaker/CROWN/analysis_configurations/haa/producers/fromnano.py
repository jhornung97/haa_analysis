from ..quantities import output as q
from ..quantities import nanoAOD as nanoAOD
from code_generation.producer import Producer, ProducerGroup

nano_jet_pt = Producer(
    name="nano_jet_pt",
    call="basefunctions::rename<ROOT::RVec<Float_t>>({df}, {input}, {output})",
    input=[nanoAOD.Jet_pt],
    output=[q.nano_jet_pt],
    scopes=['global'],
)

nano_jet_eta = Producer(
    name="nano_jet_eta",
    call="basefunctions::rename<ROOT::RVec<Float_t>>({df}, {input}, {output})",
    input=[nanoAOD.Jet_eta],
    output=[q.nano_jet_eta],
    scopes=['global'],
)

nano_jet_phi = Producer(
    name="nano_jet_phi",
    call="basefunctions::rename<ROOT::RVec<Float_t>>({df}, {input}, {output})",
    input=[nanoAOD.Jet_phi],
    output=[q.nano_jet_phi],
    scopes=['global'],
)

nano_jet_mass = Producer(
    name="nano_jet_mass",
    call="basefunctions::rename<ROOT::RVec<Float_t>>({df}, {input}, {output})",
    input=[nanoAOD.Jet_mass],
    output=[q.nano_jet_mass],
    scopes=['global'],
)

nano_jet_btagDeepB = Producer(
    name="nano_jet_btagDeepB",
    call="basefunctions::rename<ROOT::RVec<Float_t>>({df}, {input}, {output})",
    input=[nanoAOD.Jet_btagDeepB],
    output=[q.nano_jet_btagDeepB],
    scopes=['global'],
)

nano_jet_btagDeepFlavQG = Producer(
    name="nano_jet_btagDeepFlavQG",
    call="basefunctions::rename<ROOT::RVec<Float_t>>({df}, {input}, {output})",
    input=[nanoAOD.Jet_btagDeepFlavQG],
    output=[q.nano_jet_btagDeepFlavQG],
    scopes=['global'],
)

nano_jet_chEmEF = Producer(
    name="nano_jet_chEmEf",
    call="basefunctions::rename<ROOT::RVec<Float_t>>({df}, {input}, {output})",
    input=[nanoAOD.Jet_chEmEF],
    output=[q.nano_jet_chEmEF],
    scopes=['global'],
)

nano_jet_jetId = Producer(
    name="nano_jet_electronIdx1",
    call="basefunctions::rename<ROOT::RVec<Int_t>>({df}, {input}, {output})",
    input=[nanoAOD.Jet_jetId],
    output=[q.nano_jet_jetId],
    scopes=['global'],
)

nano_jet_chHEF = Producer(
    name="nano_jet_chHEF",
    call="basefunctions::rename<ROOT::RVec<Int_t>>({df}, {input}, {output})",
    input=[nanoAOD.Jet_chHEF],
    output=[q.nano_jet_chHEF],
    scopes=['global'],
)

nano_jet_neEmEF = Producer(
    name="nano_jet_nEmEF",
    call="basefunctions::rename<ROOT::RVec<Int_t>>({df}, {input}, {output})",
    input=[nanoAOD.Jet_neEmEF],
    output=[q.nano_jet_neEmEF],
    scopes=['global'],
)

nano_jet_neHEF = Producer(
    name="nano_jet_nHEF",
    call="basefunctions::rename<ROOT::RVec<Float_t>>({df}, {input}, {output})",
    input=[nanoAOD.Jet_neHEF],
    output=[q.nano_jet_neHEF],
    scopes=['global'],
)

nano_fatjet_pt = Producer(
    name="nano_fatjet_pt",
    call="basefunctions::rename<ROOT::RVec<Float_t>>({df}, {input}, {output})",
    input=[nanoAOD.FatJet_pt],
    output=[q.nano_fatjet_pt],
    scopes=['global'],
)

nano_fatjet_eta = Producer(
    name="nano_fatjet_eta",
    call="basefunctions::rename<ROOT::RVec<Float_t>>({df}, {input}, {output})",
    input=[nanoAOD.FatJet_eta],
    output=[q.nano_fatjet_eta],
    scopes=['global'],
)

nano_fatjet_phi = Producer(
    name="nano_fatjet_phi",
    call="basefunctions::rename<ROOT::RVec<Float_t>>({df}, {input}, {output})",
    input=[nanoAOD.FatJet_phi],
    output=[q.nano_fatjet_phi],
    scopes=['global'],
)

nano_fatjet_mass = Producer(
    name="nano_fatjet_mass",
    call="basefunctions::rename<ROOT::RVec<Float_t>>({df}, {input}, {output})",
    input=[nanoAOD.FatJet_mass],
    output=[q.nano_fatjet_mass],
    scopes=['global'],
)

nano_fatjet_msoftdrop = Producer(
    name="nano_fatjet_msoftdrop",
    call="basefunctions::rename<ROOT::RVec<Float_t>>({df}, {input}, {output})",
    input=[nanoAOD.FatJet_msoftdrop],
    output=[q.nano_fatjet_msoftdrop],
    scopes=['global'],
)

nano_fatjet_n2b1 = Producer(
    name="nano_fatjet_n2b1",
    call="basefunctions::rename<ROOT::RVec<Float_t>>({df}, {input}, {output})",
    input=[nanoAOD.FatJet_n2b1],
    output=[q.nano_fatjet_n2b1],
    scopes=['global'],
)

nano_fatjet_n3b1 = Producer(
    name="nano_fatjet_n3b1",
    call="basefunctions::rename<ROOT::RVec<Float_t>>({df}, {input}, {output})",
    input=[nanoAOD.FatJet_n3b1],
    output=[q.nano_fatjet_n3b1],
    scopes=['global'],
)

nano_fatjet_particleNet_H4qvsQCD = Producer(
    name="nano_fatjet_particleNet_H4qvsQCD",
    call="basefunctions::rename<ROOT::RVec<Float_t>>({df}, {input}, {output})",
    input=[nanoAOD.FatJet_particleNet_H4qvsQCD],
    output=[q.nano_fatjet_particleNet_H4qvsQCD],
    scopes=['global'],
)

nano_fatjet_particleNet_QCD = Producer(
    name="nano_fatjet_particleNet_QCD",
    call="basefunctions::rename<ROOT::RVec<Float_t>>({df}, {input}, {output})",
    input=[nanoAOD.FatJet_particleNet_QCD],
    output=[q.nano_fatjet_particleNet_QCD],
    scopes=['global'],
)

nano_fatjet_particleNet_mass = Producer(
    name="nano_fatjet_particleNet_mass",
    call="basefunctions::rename<ROOT::RVec<Float_t>>({df}, {input}, {output})",
    input=[nanoAOD.FatJet_particleNet_mass],
    output=[q.nano_fatjet_particleNet_mass],
    scopes=['global'],
)

nano_pfcands_pt = Producer(
    name="nano_pfcands_pt",
    call="basefunctions::rename<ROOT::RVec<Float_t>>({df}, {input}, {output})",
    input=[nanoAOD.PFCands_pt],
    output=[q.nano_pfcands_pt],
    scopes=['global'],
)

nano_pfcands_eta = Producer(
    name="nano_pfcands_eta",
    call="basefunctions::rename<ROOT::RVec<Float_t>>({df}, {input}, {output})",
    input=[nanoAOD.PFCands_eta],
    output=[q.nano_pfcands_eta],
    scopes=['global'],
)

nano_pfcands_phi = Producer(
    name="nano_pfcands_phi",
    call="basefunctions::rename<ROOT::RVec<Float_t>>({df}, {input}, {output})",
    input=[nanoAOD.PFCands_phi],
    output=[q.nano_pfcands_phi],
    scopes=['global'],
)

nano_pfcands_mass = Producer(
    name="nano_pfcands_mass",
    call="basefunctions::rename<ROOT::RVec<Float_t>>({df}, {input}, {output})",
    input=[nanoAOD.PFCands_mass],
    output=[q.nano_pfcands_mass],
    scopes=['global'],
)

nano_pfcands_charge = Producer(
    name="nano_pfcands_charge",
    call="basefunctions::rename<ROOT::RVec<Int_t>>({df}, {input}, {output})",
    input=[nanoAOD.PFCands_charge],
    output=[q.nano_pfcands_charge],
    scopes=['global'],
)

nano_pfcands_pdgId = Producer(
    name="nano_pfcands_pdgId",
    call="basefunctions::rename<ROOT::RVec<Int_t>>({df}, {input}, {output})",
    input=[nanoAOD.PFCands_pdgId],
    output=[q.nano_pfcands_pdgId],
    scopes=['global'],
)

nano_pfcands_d0 = Producer(
    name="nano_pfcands_d0",
    call="basefunctions::rename<ROOT::RVec<Float_t>>({df}, {input}, {output})",
    input=[nanoAOD.PFCands_d0],
    output=[q.nano_pfcands_d0],
    scopes=['global'],
)

nano_pfcands_d0Err = Producer(
    name="nano_pfcands_d0Err",
    call="basefunctions::rename<ROOT::RVec<Float_t>>({df}, {input}, {output})",
    input=[nanoAOD.PFCands_d0Err],
    output=[q.nano_pfcands_d0Err],
    scopes=['global'],
)

nano_pfcands_dz = Producer(
    name="nano_pfcands_dz",
    call="basefunctions::rename<ROOT::RVec<Float_t>>({df}, {input}, {output})",
    input=[nanoAOD.PFCands_dz],
    output=[q.nano_pfcands_dz],
    scopes=['global'],
)

nano_pfcands_dzErr = Producer(
    name="nano_pfcands_dzErr",
    call="basefunctions::rename<ROOT::RVec<Float_t>>({df}, {input}, {output})",
    input=[nanoAOD.PFCands_dzErr],
    output=[q.nano_pfcands_dzErr],
    scopes=['global'],
)

nano_pfcands_trkQuality = Producer(
    name="nano_pfcands_trkQuality",
    call="basefunctions::rename<ROOT::RVec<Int_t>>({df}, {input}, {output})",
    input=[nanoAOD.PFCands_trkQuality],
    output=[q.nano_pfcands_trkQuality],
    scopes=['global'],
)

copy = ProducerGroup(
	name = "copy",
	call = None,
	input = None,
	output = None,
	scopes = ["global"],
subproducers = [nano_jet_pt, nano_jet_eta, nano_jet_phi, nano_jet_mass, nano_jet_btagDeepB, nano_jet_btagDeepFlavQG, nano_jet_chEmEF, nano_jet_jetId, nano_jet_chHEF, nano_jet_neEmEF, nano_jet_neHEF, nano_fatjet_pt, nano_fatjet_eta, nano_fatjet_phi, nano_fatjet_mass, nano_fatjet_msoftdrop, nano_fatjet_n2b1, nano_fatjet_n3b1, nano_fatjet_particleNet_H4qvsQCD, nano_fatjet_particleNet_QCD, nano_fatjet_particleNet_mass, nano_pfcands_pt, nano_pfcands_eta, nano_pfcands_phi, nano_pfcands_mass, nano_pfcands_charge, nano_pfcands_pdgId, nano_pfcands_d0, nano_pfcands_d0Err, nano_pfcands_dz, nano_pfcands_dzErr, nano_pfcands_trkQuality],
)
