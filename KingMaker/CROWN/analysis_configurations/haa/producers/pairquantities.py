from ..quantities import output as q
from ..quantities import nanoAOD as nanoAOD
from code_generation.producer import Producer, ProducerGroup, ExtendedVectorProducer

####################
# Set of general producers for DiTauPair Quantities
####################

pt_1 = Producer(
    name="pt_1",
    call="quantities::pt({df}, {output}, {input})",
    input=[q.p4_1],
    output=[q.pt_1],
    scopes=["mm","ee","em"],
)
pt_2 = Producer(
    name="pt_2",
    call="quantities::pt({df}, {output}, {input})",
    input=[q.p4_2],
    output=[q.pt_2],
    scopes=["mm","ee","em"],
)
eta_1 = Producer(
    name="eta_1",
    call="quantities::eta({df}, {output}, {input})",
    input=[q.p4_1],
    output=[q.eta_1],
    scopes=["mm","ee","em"],
)
eta_2 = Producer(
    name="eta_2",
    call="quantities::eta({df}, {output}, {input})",
    input=[q.p4_2],
    output=[q.eta_2],
    scopes=["mm","ee","em"],
)
phi_1 = Producer(
    name="phi_1",
    call="quantities::phi({df}, {output}, {input})",
    input=[q.p4_1],
    output=[q.phi_1],
    scopes=["mm","ee","em"],
)
phi_2 = Producer(
    name="phi_2",
    call="quantities::phi({df}, {output}, {input})",
    input=[q.p4_2],
    output=[q.phi_2],
    scopes=["mm","ee","em"],
)
mass_1 = Producer(
    name="mass_1",
    call="quantities::mass({df}, {output}, {input})",
    input=[q.p4_1],
    output=[q.mass_1],
    scopes=["mm","ee","em"],
)
mass_2 = Producer(
    name="mass_2",
    call="quantities::mass({df}, {output}, {input})",
    input=[q.p4_2],
    output=[q.mass_2],
    scopes=["mm","ee","em"],
)
m_vis = Producer(
    name="m_vis",
    call="quantities::m_vis({df}, {output}, {input_vec})",
    input=[q.p4_1, q.p4_2],
    output=[q.m_vis],
    scopes=["mm","ee","em"],
)
pt_vis = Producer(
    name="pt_vis",
    call="quantities::pt_vis({df}, {output}, {input_vec})",
    input=[q.p4_1, q.p4_2],
    output=[q.pt_vis],
    scopes=["mm","ee","em"],
)

eta_vis = Producer(
    name="eta_vis",
    call="quantities::eta_vis({df}, {output}, {input_vec})",
    input=[q.p4_1, q.p4_2],
    output=[q.eta_vis],
    scopes=["mm","ee","em"],
)

phi_vis = Producer(
    name="phi_vis",
    call="quantities::phi_vis({df}, {output}, {input_vec})",
    input=[q.p4_1, q.p4_2],
    output=[q.phi_vis],
    scopes=["mm","ee","em"],
)
####################
# Set of channel specific producers
####################
muon_dxy_1 = Producer(
    name="muon_dxy_1",
    call="quantities::dxy({df}, {output}, 0, {input})",
    input=[q.dileptonpair, nanoAOD.Muon_dxy],
    output=[q.dxy_1],
    scopes=["mm"],
)
muon_dxy_2 = Producer(
    name="muon_dxy_2",
    call="quantities::dxy({df}, {output}, 1, {input})",
    input=[q.dileptonpair, nanoAOD.Muon_dxy],
    output=[q.dxy_2],
    scopes=["mm","em"],
)
muon_dz_1 = Producer(
    name="muon_dz_1",
    call="quantities::dz({df}, {output}, 0, {input})",
    input=[q.dileptonpair, nanoAOD.Muon_dz],
    output=[q.dz_1],
    scopes=["mm"],
)
muon_dz_2 = Producer(
    name="muon_dz_2",
    call="quantities::dz({df}, {output}, 1, {input})",
    input=[q.dileptonpair, nanoAOD.Muon_dz],
    output=[q.dz_2],
    scopes=["mm","em"],
)
muon_q_1 = Producer(
    name="muon_q_1",
    call="quantities::charge({df}, {output}, 0, {input})",
    input=[q.dileptonpair, nanoAOD.Muon_charge],
    output=[q.q_1],
    scopes=["mm"],
)
muon_q_2 = Producer(
    name="muon_q_2",
    call="quantities::charge({df}, {output}, 1, {input})",
    input=[q.dileptonpair, nanoAOD.Muon_charge],
    output=[q.q_2],
    scopes=["mm","em"],
)
muon_iso_1 = Producer(
    name="muon_iso_1",
    call="quantities::isolation({df}, {output}, 0, {input})",
    input=[q.dileptonpair, nanoAOD.Muon_iso],
    output=[q.iso_1],
    scopes=["mm"],
)
muon_iso_2 = Producer(
    name="muon_iso_2",
    call="quantities::isolation({df}, {output}, 1, {input})",
    input=[q.dileptonpair, nanoAOD.Muon_iso],
    output=[q.iso_2],
    scopes=["mm","em"],
)

UnrollMuLV1 = ProducerGroup(
    name="UnrollMuLV1",
    call=None,
    input=None,
    output=None,
    scopes=["mm"],
    subproducers=[
        pt_1,
        eta_1,
        phi_1,
        mass_1,
        muon_dxy_1,
        muon_dz_1,
        muon_q_1,
        muon_iso_1,
    ],
)
UnrollMuLV2 = ProducerGroup(
    name="UnrollMuLV2",
    call=None,
    input=None,
    output=None,
    scopes=["mm","em"],
    subproducers=[
        pt_2,
        eta_2,
        phi_2,
        mass_2,
        muon_dxy_2,
        muon_dz_2,
        muon_q_2,
        muon_iso_2,
    ],
)

electron_dxy_1 = Producer(
    name="electron_dxy_1",
    call="quantities::dxy({df}, {output}, 0, {input})",
    input=[q.dileptonpair, nanoAOD.Electron_dxy],
    output=[q.dxy_1],
    scopes=["ee","em"],
)
electron_dxy_2 = Producer(
    name="electron_dxy_2",
    call="quantities::dxy({df}, {output}, 1, {input})",
    input=[q.dileptonpair, nanoAOD.Electron_dxy],
    output=[q.dxy_2],
    scopes=["ee"],
)
electron_dz_1 = Producer(
    name="electron_dz_1",
    call="quantities::dz({df}, {output}, 0, {input})",
    input=[q.dileptonpair, nanoAOD.Electron_dz],
    output=[q.dz_1],
    scopes=["ee","em"],
)
electron_dz_2 = Producer(
    name="electron_dz_2",
    call="quantities::dz({df}, {output}, 1, {input})",
    input=[q.dileptonpair, nanoAOD.Electron_dz],
    output=[q.dz_2],
    scopes=["ee"],
)
electron_q_1 = Producer(
    name="electron_q_1",
    call="quantities::charge({df}, {output}, 0, {input})",
    input=[q.dileptonpair, nanoAOD.Electron_charge],
    output=[q.q_1],
    scopes=["ee","em"],
)
electron_q_2 = Producer(
    name="electron_q_2",
    call="quantities::charge({df}, {output}, 1, {input})",
    input=[q.dileptonpair, nanoAOD.Electron_charge],
    output=[q.q_2],
    scopes=["ee"],
)
electron_iso_1 = Producer(
    name="electron_iso_1",
    call="quantities::isolation({df}, {output}, 0, {input})",
    input=[q.dileptonpair, nanoAOD.Electron_iso],
    output=[q.iso_1],
    scopes=["ee","em"],
)
electron_iso_2 = Producer(
    name="electron_iso_2",
    call="quantities::isolation({df}, {output}, 1, {input})",
    input=[q.dileptonpair, nanoAOD.Electron_iso],
    output=[q.iso_2],
    scopes=["ee"],
)

UnrollElLV1 = ProducerGroup(
    name="UnrollElLV1",
    call=None,
    input=None,
    output=None,
    scopes=["ee","em"],
    subproducers=[
        pt_1,
        eta_1,
        phi_1,
        mass_1,
        electron_dxy_1,
        electron_dz_1,
        electron_q_1,
        electron_iso_1,
    ],
)
UnrollElLV2 = ProducerGroup(
    name="UnrollElLV2",
    call=None,
    input=None,
    output=None,
    scopes=["ee"],
    subproducers=[
        pt_2,
        eta_2,
        phi_2,
        mass_2,
        electron_dxy_2,
        electron_dz_2,
        electron_q_2,
        electron_iso_2,
    ],
)
#####################
# Producer Groups
#####################

MMDiTauPairQuantities = ProducerGroup(
    name="MMDiTauPairQuantities",
    call=None,
    input=None,
    output=None,
    scopes=["mm"],
    subproducers=[UnrollMuLV1, UnrollMuLV2, m_vis, pt_vis, eta_vis, phi_vis,],
)

EEDiTauPairQuantities = ProducerGroup(
    name="EEDiTauPairQuantities",
    call=None,
    input=None,
    output=None,
    scopes=["ee"],
    subproducers=[UnrollElLV1, UnrollElLV2, m_vis, pt_vis, eta_vis, phi_vis,],
)

EMDiTauPairQuantities = ProducerGroup(
    name="EMDiTauPairQuantities",
    call=None,
    input=None,
    output=None,
    scopes=["em"],
    subproducers=[UnrollElLV1, UnrollMuLV2, m_vis, pt_vis, eta_vis, phi_vis,],
)
