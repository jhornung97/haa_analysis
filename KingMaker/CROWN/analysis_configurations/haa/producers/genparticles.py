from ..quantities import output as q
from ..quantities import nanoAOD as nanoAOD
from code_generation.producer import Producer, ProducerGroup

####################
# Set of producers to get the genParticles from the ditaupair
####################
MMGenPair = Producer(
    name="MMGenPair",
    call="ditau_pairselection::buildgenpair({df}, {input}, {output})",
    input=[q.dileptonpair, nanoAOD.Muon_indexToGen, nanoAOD.Muon_indexToGen],
    output=[q.gen_dileptonpair],
    scopes=["mm"],
)
EEGenPair = Producer(
    name="EEGenPair",
    call="ditau_pairselection::buildgenpair({df}, {input}, {output})",
    input=[q.dileptonpair, nanoAOD.Electron_indexToGen, nanoAOD.Electron_indexToGen],
    output=[q.gen_dileptonpair],
    scopes=["ee"],
)
EMGenPair = Producer(
    name="EMGenPair",
    call="ditau_pairselection::buildgenpair({df}, {input}, {output})",
    input=[q.dileptonpair, nanoAOD.Electron_indexToGen, nanoAOD.Muon_indexToGen],
    output=[q.gen_dileptonpair],
    scopes=["em"],
)
MMTrueGenPair = Producer(
    name="GenPair",
    call="ditau_pairselection::buildtruegenpair({df}, {input}, {output}, {truegen_mother_pdgid}, {truegen_daughter_1_pdgid}, {truegen_daughter_2_pdgid})",
    input=[
        nanoAOD.GenParticle_statusFlags,
        nanoAOD.GenParticle_status,
        nanoAOD.GenParticle_pdgId,
        nanoAOD.GenParticle_motherid,
        nanoAOD.GenParticle_pt,
    ],
    output=[q.truegenpair],
    scopes=["mm","ee","em"],
)
####################
# Set of general producers for Gen DiTauPair Quantities
####################

LVGenParticle1 = Producer(
    name="LVGenParticle1",
    call="lorentzvectors::build({df}, {input_vec}, 0, {output})",
    input=[
        q.gen_dileptonpair,
        nanoAOD.GenParticle_pt,
        nanoAOD.GenParticle_eta,
        nanoAOD.GenParticle_phi,
        nanoAOD.GenParticle_mass,
    ],
    output=[q.gen_p4_1],
    scopes=["mm","ee","em"],
)
LVGenParticle2 = Producer(
    name="LVGenParticle2",
    call="lorentzvectors::build({df}, {input_vec}, 1, {output})",
    input=[
        q.gen_dileptonpair,
        nanoAOD.GenParticle_pt,
        nanoAOD.GenParticle_eta,
        nanoAOD.GenParticle_phi,
        nanoAOD.GenParticle_mass,
    ],
    output=[q.gen_p4_2],
    scopes=["mm","ee","em"],
)
LVTrueGenParticle1 = Producer(
    name="LVTrueGenParticle1",
    call="lorentzvectors::build({df}, {input_vec}, 0, {output})",
    input=[
        q.truegenpair,
        nanoAOD.GenParticle_pt,
        nanoAOD.GenParticle_eta,
        nanoAOD.GenParticle_phi,
        nanoAOD.GenParticle_mass,
    ],
    output=[q.gen_p4_1],
    scopes=["mm","ee","em"],
)
LVTrueGenParticle2 = Producer(
    name="LVTrueGenParticle2",
    call="lorentzvectors::build({df}, {input_vec}, 1, {output})",
    input=[
        q.truegenpair,
        nanoAOD.GenParticle_pt,
        nanoAOD.GenParticle_eta,
        nanoAOD.GenParticle_phi,
        nanoAOD.GenParticle_mass,
    ],
    output=[q.gen_p4_2],
    scopes=["mm","ee","em"],
)
gen_pt_1 = Producer(
    name="gen_pt_1",
    call="quantities::pt({df}, {output}, {input})",
    input=[q.gen_p4_1],
    output=[q.gen_pt_1],
    scopes=["mm","ee","em"],
)
gen_pt_2 = Producer(
    name="gen_pt_2",
    call="quantities::pt({df}, {output}, {input})",
    input=[q.gen_p4_2],
    output=[q.gen_pt_2],
    scopes=["mm","ee","em"],
)
gen_eta_1 = Producer(
    name="gen_eta_1",
    call="quantities::eta({df}, {output}, {input})",
    input=[q.gen_p4_1],
    output=[q.gen_eta_1],
    scopes=["mm","ee","em"],
)
gen_eta_2 = Producer(
    name="gen_eta_2",
    call="quantities::eta({df}, {output}, {input})",
    input=[q.gen_p4_2],
    output=[q.gen_eta_2],
    scopes=["mm","ee","em"],
)
gen_phi_1 = Producer(
    name="gen_phi_1",
    call="quantities::phi({df}, {output}, {input})",
    input=[q.gen_p4_1],
    output=[q.gen_phi_1],
    scopes=["mm","ee","em"],
)
gen_phi_2 = Producer(
    name="gen_phi_2",
    call="quantities::phi({df}, {output}, {input})",
    input=[q.gen_p4_2],
    output=[q.gen_phi_2],
    scopes=["mm","ee","em"],
)
gen_mass_1 = Producer(
    name="gen_mass_1",
    call="quantities::mass({df}, {output}, {input})",
    input=[q.gen_p4_1],
    output=[q.gen_mass_1],
    scopes=["mm","ee","em"],
)
gen_mass_2 = Producer(
    name="gen_mass_2",
    call="quantities::mass({df}, {output}, {input})",
    input=[q.gen_p4_2],
    output=[q.gen_mass_2],
    scopes=["mm","ee","em"],
)
gen_pdgid_1 = Producer(
    name="gen_pdgid_1",
    call="quantities::pdgid({df}, {output}, 0, {input})",
    input=[q.gen_dileptonpair, nanoAOD.GenParticle_pdgId],
    output=[q.gen_pdgid_1],
    scopes=["mm","ee","em"],
)
gen_pdgid_2 = Producer(
    name="gen_pdgid_2",
    call="quantities::pdgid({df}, {output}, 1, {input})",
    input=[q.gen_dileptonpair, nanoAOD.GenParticle_pdgId],
    output=[q.gen_pdgid_2],
    scopes=["mm","ee","em"],
)
gen_m_vis = Producer(
    name="gen_m_vis",
    call="quantities::m_vis({df}, {output}, {input_vec})",
    input=[q.gen_p4_1, q.gen_p4_2],
    output=[q.gen_m_vis],
    scopes=["mm","ee","em"],
)
gen_pt_vis = Producer(
    name="gen_pt_vis",
    call="quantities::pt_vis({df}, {output}, {input_vec})",
    input=[q.gen_p4_1, q.gen_p4_2],
    output=[q.gen_pt_vis],
    scopes=["mm","ee","em"],
)
gen_eta_vis = Producer(
    name="gen_eta_vis",
    call="quantities::eta_vis({df}, {output}, {input_vec})",
    input=[q.gen_p4_1, q.gen_p4_2],
    output=[q.gen_eta_vis],
    scopes=["mm","ee","em"],
)
gen_match_2 = Producer(
    name="gen_match_2",
    call="quantities::tau::genmatch({df}, {output}, 1, {input})",
    input=[q.dileptonpair, nanoAOD.Tau_genMatch],
    output=[q.gen_match_2],
    scopes=["mm","ee","em"],
)
gen_taujet_pt_1 = Producer(
    name="gen_taujet_pt_1",
    call="quantities::tau::matching_genjet_pt({df}, {output}, 0, {input})",
    input=[
        q.dileptonpair,
        nanoAOD.Tau_associatedJet,
        nanoAOD.Jet_associatedGenJet,
        nanoAOD.GenJet_pt,
    ],
    output=[q.gen_taujet_pt_1],
    scopes=["mm","ee","em"],
)
gen_taujet_pt_2 = Producer(
    name="gen_taujet_pt_2",
    call="quantities::tau::matching_genjet_pt({df}, {output}, 1, {input})",
    input=[
        q.dileptonpair,
        nanoAOD.Tau_associatedJet,
        nanoAOD.Jet_associatedGenJet,
        nanoAOD.GenJet_pt,
    ],
    output=[q.gen_taujet_pt_2],
    scopes=["mm","ee","em"],
)
UnrollGenMuLV1 = ProducerGroup(
    name="UnrollGenMuLV1",
    call=None,
    input=None,
    output=None,
    scopes=["mm","ee","em"],
    subproducers=[gen_pt_1, gen_eta_1, gen_phi_1, gen_mass_1,] #gen_pdgid_1],
)
UnrollGenMuLV2 = ProducerGroup(
    name="UnrollGenMuLV2",
    call=None,
    input=None,
    output=None,
    scopes=["mm","ee","em"],
    subproducers=[gen_pt_2, gen_eta_2, gen_phi_2, gen_mass_2,] #gen_pdgid_2],
)


MMGenDiTauPairQuantities = ProducerGroup(
    name="MMGenDiTauPairQuantities",
    call=None,
    input=None,
    output=None,
    scopes=["mm"],
    subproducers=[
        MMGenPair,
        LVGenParticle1,
        LVGenParticle2,
        UnrollGenMuLV1,
        UnrollGenMuLV2,
        gen_m_vis,
        gen_pt_vis,
        gen_eta_vis,
    ],
)
EEGenDiTauPairQuantities = ProducerGroup(
    name="EEGenDiTauPairQuantities",
    call=None,
    input=None,
    output=None,
    scopes=["ee"],
    subproducers=[
        EEGenPair,
        LVGenParticle1,
        LVGenParticle2,
        UnrollGenMuLV1,
        UnrollGenMuLV2,
        gen_m_vis,
        gen_pt_vis,
        gen_eta_vis,
    ],
)
EMGenDiTauPairQuantities = ProducerGroup(
    name="EMGenDiTauPairQuantities",
    call=None,
    input=None,
    output=None,
    scopes=["em"],
    subproducers=[
        EMGenPair,
        LVGenParticle1,
        LVGenParticle2,
        UnrollGenMuLV1,
        UnrollGenMuLV2,
        gen_m_vis,
        gen_pt_vis,
        gen_eta_vis,
    ],
)
MMTrueGenDiTauPairQuantities = ProducerGroup(
    name="GenDiTauPairQuantities",
    call=None,
    input=None,
    output=None,
    scopes=["mm","ee","em"],
    subproducers=[
        MMTrueGenPair,
        LVTrueGenParticle1,
        LVTrueGenParticle2,
        UnrollGenMuLV1,
        UnrollGenMuLV2,
        gen_m_vis,
        gen_pt_vis,
    ],
)
