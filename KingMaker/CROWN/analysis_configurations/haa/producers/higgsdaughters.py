from ..quantities import output as q
from ..quantities import nanoAOD as nanoAOD
from code_generation.producer import Producer, ProducerGroup

Daughter_Iso = Producer(
    name = "Daughter_Iso",
    call = "haa::pfCandIso({df}, {output}, {input})",
    input = [
        nanoAOD.PFCands_pt, 
        nanoAOD.PFCands_eta, 
        nanoAOD.PFCands_phi, 
        nanoAOD.PFCands_mass,
        q.higgsdaughters,
        q.ChargedPFCands,
        q.NeutralPFCands,
        q.PhotonPFCands,
        q.fromPV,
        ],
    output = [q.daughter_isos],
    scopes = ["mm", "ee", "em"],
)

# Get kinematics of hardest hadron

LVPFCand1 = Producer(
    name="LVPFCand1",
    call="lorentzvectors::buildFromPFCand({df}, {input_vec}, 0, {output})",
    input=[
        q.higgsdaughters,
        nanoAOD.PFCands_pt,
        nanoAOD.PFCands_eta,
        nanoAOD.PFCands_phi,
        nanoAOD.PFCands_mass,
    ],
    output=[q.d1_p4],
    scopes=["mm", "ee", "em"],
)

d1_pt = Producer(
    name="d1_pt",
    call="quantities::pt({df}, {output}, {input})",
    input=[q.d1_p4],
    output=[q.d1_pt],
    scopes=["mm", "ee", "em"],
)

d1_eta = Producer(
    name="d1_eta",
    call="quantities::eta({df}, {output}, {input})",
    input=[q.d1_p4],
    output=[q.d1_eta],
    scopes=["mm", "ee", "em"],
)

d1_phi = Producer(
    name="d1_phi",
    call="quantities::phi({df}, {output}, {input})",
    input=[q.d1_p4],
    output=[q.d1_phi],
    scopes=["mm", "ee", "em"],
)
d1_mass = Producer(
    name="d1_mass",
    call="quantities::mass({df}, {output}, {input})",
    input=[q.d1_p4],
    output=[q.d1_mass],
    scopes=["mm", "ee", "em"],
)
d1_promptFlag = Producer(
    name="d1_promptFlag",
    call="quantities::charge({df}, {output}, 0, {input})",
    input=[q.higgsdaughters, nanoAOD.PFCands_fromPV],
    output=[q.d1_prompt],
    scopes=["mm", "ee", "em"],
)
d1_iso = Producer(
    name="d1_iso",
    call = "quantities::iso({df}, {output}, 0, {input})",
    input = [q.daughter_isos],
    output = [q.d1_iso],
    scopes = ["mm", "ee", "em"],
)
d1Quantities = ProducerGroup(
	name = "d1Quantities",
	call = None,
	input = None,
	output = None,
	scopes = ["mm", "ee", "em"],
	subproducers = [LVPFCand1, d1_pt, d1_eta, d1_phi, d1_mass, d1_promptFlag, d1_iso],
)

# Get kinematics of 2nd hardest hadron

LVPFCand2 = Producer(
    name="LVPFCand2",
    call="lorentzvectors::buildFromPFCand({df}, {input_vec}, 1, {output})",
    input=[
        q.higgsdaughters,
        nanoAOD.PFCands_pt,
        nanoAOD.PFCands_eta,
        nanoAOD.PFCands_phi,
        nanoAOD.PFCands_mass,
    ],
    output=[q.d2_p4],
    scopes=["mm", "ee", "em"],
)


d2_pt = Producer(
    name="d2_pt",
    call="quantities::pt({df}, {output}, {input})",
    input=[q.d2_p4],
    output=[q.d2_pt],
    scopes=["mm", "ee", "em"],
)

d2_eta = Producer(
    name="d2_eta",
    call="quantities::eta({df}, {output}, {input})",
    input=[q.d2_p4],
    output=[q.d2_eta],
    scopes=["mm", "ee", "em"],
)

d2_phi = Producer(
    name="d2_phi",
    call="quantities::phi({df}, {output}, {input})",
    input=[q.d2_p4],
    output=[q.d2_phi],
    scopes=["mm", "ee", "em"],
)
d2_mass = Producer(
    name="d2_mass",
    call="quantities::mass({df}, {output}, {input})",
    input=[q.d2_p4],
    output=[q.d2_mass],
    scopes=["mm", "ee", "em"],
)
d2_promptFlag = Producer(
    name="d1_promptFlag",
    call="quantities::charge({df}, {output}, 1, {input})",
    input=[q.higgsdaughters, nanoAOD.PFCands_fromPV],
    output=[q.d2_prompt],
    scopes=["mm", "ee", "em"],
)
d2_iso = Producer(
    name="d2_iso",
    call = "quantities::iso({df}, {output}, 1, {input})",
    input = [q.daughter_isos],
    output = [q.d2_iso],
    scopes = ["mm", "ee", "em"],
)
d2Quantities = ProducerGroup(
	name = "d2Quantities",
	call = None,
	input = None,
	output = None,
	scopes = ["mm", "ee", "em"],
	subproducers = [LVPFCand2, d2_pt, d2_eta, d2_phi, d2_mass, d2_promptFlag, d2_iso],
)

# Get kinematics of 3rd hardest hadron

LVPFCand3 = Producer(
    name="LVPFCand3",
    call="lorentzvectors::buildFromPFCand({df}, {input_vec}, 2, {output})",
    input=[
        q.higgsdaughters,
        nanoAOD.PFCands_pt,
        nanoAOD.PFCands_eta,
        nanoAOD.PFCands_phi,
        nanoAOD.PFCands_mass,
    ],
    output=[q.d3_p4],
    scopes=["mm", "ee", "em"],
)


d3_pt = Producer(
    name="d3_pt",
    call="quantities::pt({df}, {output}, {input})",
    input=[q.d3_p4],
    output=[q.d3_pt],
    scopes=["mm", "ee", "em"],
)

d3_eta = Producer(
    name="d3_eta",
    call="quantities::eta({df}, {output}, {input})",
    input=[q.d3_p4],
    output=[q.d3_eta],
    scopes=["mm", "ee", "em"],
)

d3_phi = Producer(
    name="d3_phi",
    call="quantities::phi({df}, {output}, {input})",
    input=[q.d3_p4],
    output=[q.d3_phi],
    scopes=["mm", "ee", "em"],
)
d3_mass = Producer(
    name="d3_mass",
    call="quantities::mass({df}, {output}, {input})",
    input=[q.d3_p4],
    output=[q.d3_mass],
    scopes=["mm", "ee", "em"],
)
d3_promptFlag = Producer(
    name="d1_promptFlag",
    call="quantities::charge({df}, {output}, 2, {input})",
    input=[q.higgsdaughters, nanoAOD.PFCands_fromPV],
    output=[q.d3_prompt],
    scopes=["mm", "ee", "em"],
)
d3_iso = Producer(
    name="d3_iso",
    call = "quantities::iso({df}, {output}, 2, {input})",
    input = [q.daughter_isos],
    output = [q.d3_iso],
    scopes = ["mm", "ee", "em"],
)
d3Quantities = ProducerGroup(
	name = "d3Quantities",
	call = None,
	input = None,
	output = None,
	scopes = ["mm", "ee", "em"],
	subproducers = [LVPFCand3, d3_pt, d3_eta, d3_phi, d3_mass, d3_promptFlag, d3_iso],
)

# Get kinematics of 4th hardest hadron

LVPFCand4 = Producer(
    name="LVPFCand4",
    call="lorentzvectors::buildFromPFCand({df}, {input_vec}, 3, {output})",
    input=[
        q.higgsdaughters,
        nanoAOD.PFCands_pt,
        nanoAOD.PFCands_eta,
        nanoAOD.PFCands_phi,
        nanoAOD.PFCands_mass,
    ],
    output=[q.d4_p4],
    scopes=["mm", "ee", "em"],
)


d4_pt = Producer(
    name="d4_pt",
    call="quantities::pt({df}, {output}, {input})",
    input=[q.d4_p4],
    output=[q.d4_pt],
    scopes=["mm", "ee", "em"],
)

d4_eta = Producer(
    name="d4_eta",
    call="quantities::eta({df}, {output}, {input})",
    input=[q.d4_p4],
    output=[q.d4_eta],
    scopes=["mm", "ee", "em"],
)

d4_phi = Producer(
    name="d4_phi",
    call="quantities::phi({df}, {output}, {input})",
    input=[q.d4_p4],
    output=[q.d4_phi],
    scopes=["mm", "ee", "em"],
)
d4_mass = Producer(
    name="d4_mass",
    call="quantities::mass({df}, {output}, {input})",
    input=[q.d4_p4],
    output=[q.d4_mass],
    scopes=["mm", "ee", "em"],
)
d4_promptFlag = Producer(
    name="d4_promptFlag",
    call="quantities::charge({df}, {output}, 3, {input})",
    input=[q.higgsdaughters, nanoAOD.PFCands_fromPV],
    output=[q.d4_prompt],
    scopes=["mm", "ee", "em"],
)
d4_iso = Producer(
    name="d4_iso",
    call = "quantities::iso({df}, {output}, 3, {input})",
    input = [q.daughter_isos],
    output = [q.d4_iso],
    scopes = ["mm", "ee", "em"],
)
d4Quantities = ProducerGroup(
	name = "d4Quantities",
	call = None,
	input = None,
	output = None,
	scopes = ["mm", "ee", "em"],
	subproducers = [LVPFCand4, d4_pt, d4_eta, d4_phi, d4_mass, d4_promptFlag, d4_iso],
)
dQuantities = ProducerGroup(
	name = "dQuantities",
	call = None,
	input = None,
	output = None,
	scopes = ["mm", "ee", "em"],
	subproducers = [Daughter_Iso, d1Quantities, d2Quantities, d3Quantities, d4Quantities],
)

# Sorted by pseudoscalar pt

ps_1_d_1_p4 = Producer(
    name="ps_1_d_1_p4",
    call="lorentzvectors::buildFromPFCand({df}, {input_vec}, 0, {output})",
    input=[
            q.ps1Pair, 
            nanoAOD.PFCands_pt,
            nanoAOD.PFCands_eta,
            nanoAOD.PFCands_phi,
            nanoAOD.PFCands_mass,
        ],
    output=[q.ps_1_d_1_p4],
    scopes=["mm", "ee", "em"],
)

ps_1_d_2_p4 = Producer(
    name="ps_1_d_2_p4",
    call="lorentzvectors::buildFromPFCand({df}, {input_vec}, 1, {output})",
    input=[
            q.ps1Pair, 
            nanoAOD.PFCands_pt,
            nanoAOD.PFCands_eta,
            nanoAOD.PFCands_phi,
            nanoAOD.PFCands_mass,
        ],
    output=[q.ps_1_d_2_p4],
    scopes=["mm", "ee", "em"],
)

ps_2_d_1_p4 = Producer(
    name="ps_2_d_1_p4",
    call="lorentzvectors::buildFromPFCand({df}, {input_vec}, 0, {output})",
    input=[
            q.ps2Pair, 
            nanoAOD.PFCands_pt,
            nanoAOD.PFCands_eta,
            nanoAOD.PFCands_phi,
            nanoAOD.PFCands_mass,
        ],
    output=[q.ps_2_d_1_p4],
    scopes=["mm", "ee", "em"],
)

ps_2_d_2_p4 = Producer(
    name="ps_2_d_2_p4",
    call="lorentzvectors::buildFromPFCand({df}, {input_vec}, 1, {output})",
    input=[
            q.ps2Pair, 
            nanoAOD.PFCands_pt,
            nanoAOD.PFCands_eta,
            nanoAOD.PFCands_phi,
            nanoAOD.PFCands_mass,
        ],
    output=[q.ps_2_d_2_p4],
    scopes=["mm", "ee", "em"],
)

ps_daughters = ProducerGroup(
    name = "ps_daughters",
    call = None,
    input = None,
    output = None,
    scopes = ["mm", "ee", "em"],
    subproducers = [ps_1_d_1_p4, ps_1_d_2_p4, ps_2_d_1_p4, ps_2_d_2_p4],
)

# Get p4s of true higgs daughters

GetTrueDaughterP4s = Producer(
	name = "GetTrueDaughterP4s",
	call = "haa::GetTrueDaughterP4s({df}, {input}, {output})",
	input = [
			    nanoAOD.nGenParticle,
                nanoAOD.GenParticle_pdgId,
                nanoAOD.GenParticle_motherid,
			    nanoAOD.GenParticle_pt, 
			    nanoAOD.GenParticle_eta, 
			    nanoAOD.GenParticle_phi, 
			    nanoAOD.GenParticle_mass,
                         ],
	output = [q.truth_d1_p4, q.truth_d2_p4, q.truth_d3_p4, q.truth_d4_p4, q.truth_H_p4],
	scopes = ["mm", "ee", "em"],
)

# Get kinematics of hardest hadron

truth_d1_pt = Producer(
    name="truth_d1_pt",
    call="quantities::pt({df}, {output}, {input})",
    input=[q.truth_d1_p4],
    output=[q.truth_d1_pt],
    scopes=["mm", "ee", "em"],
)

truth_d1_eta = Producer(
    name="truth_d1_eta",
    call="quantities::eta({df}, {output}, {input})",
    input=[q.truth_d1_p4],
    output=[q.truth_d1_eta],
    scopes=["mm", "ee", "em"],
)

truth_d1_phi = Producer(
    name="truth_d1_phi",
    call="quantities::phi({df}, {output}, {input})",
    input=[q.truth_d1_p4],
    output=[q.truth_d1_phi],
    scopes=["mm", "ee", "em"],
)

truth_d1Quantities = ProducerGroup(
	name = "truth_d1Quantities",
	call = None,
	input = None,
	output = None,
	scopes = ["mm", "ee", "em"],
	subproducers = [truth_d1_pt, truth_d1_eta, truth_d1_phi],
)

# Get kinematics of 2nd hardest hadron

truth_d2_pt = Producer(
    name="truth_d2_pt",
    call="quantities::pt({df}, {output}, {input})",
    input=[q.truth_d2_p4],
    output=[q.truth_d2_pt],
    scopes=["mm", "ee", "em"],
)

truth_d2_eta = Producer(
    name="truth_d2_eta",
    call="quantities::eta({df}, {output}, {input})",
    input=[q.truth_d2_p4],
    output=[q.truth_d2_eta],
    scopes=["mm", "ee", "em"],
)

truth_d2_phi = Producer(
    name="truth_d2_phi",
    call="quantities::phi({df}, {output}, {input})",
    input=[q.truth_d2_p4],
    output=[q.truth_d2_phi],
    scopes=["mm", "ee", "em"],
)

truth_d2Quantities = ProducerGroup(
	name = "truth_d2Quantities",
	call = None,
	input = None,
	output = None,
	scopes = ["mm", "ee", "em"],
	subproducers = [truth_d2_pt, truth_d2_eta, truth_d2_phi],
)

# Get kinematics of 3rd hardest hadron

truth_d3_pt = Producer(
    name="truth_d3_pt",
    call="quantities::pt({df}, {output}, {input})",
    input=[q.truth_d3_p4],
    output=[q.truth_d3_pt],
    scopes=["mm", "ee", "em"],
)

truth_d3_eta = Producer(
    name="truth_d3_eta",
    call="quantities::eta({df}, {output}, {input})",
    input=[q.truth_d3_p4],
    output=[q.truth_d3_eta],
    scopes=["mm", "ee", "em"],
)

truth_d3_phi = Producer(
    name="truth_d3_phi",
    call="quantities::phi({df}, {output}, {input})",
    input=[q.truth_d3_p4],
    output=[q.truth_d3_phi],
    scopes=["mm", "ee", "em"],
)

truth_d3Quantities = ProducerGroup(
	name = "truth_d3Quantities",
	call = None,
	input = None,
	output = None,
	scopes = ["mm", "ee", "em"],
	subproducers = [truth_d3_pt, truth_d3_eta, truth_d3_phi],
)

# Get kinematics of 4th hardest hadron

truth_d4_pt = Producer(
    name="truth_d4_pt",
    call="quantities::pt({df}, {output}, {input})",
    input=[q.truth_d4_p4],
    output=[q.truth_d4_pt],
    scopes=["mm", "ee", "em"],
)

truth_d4_eta = Producer(
    name="truth_d4_eta",
    call="quantities::eta({df}, {output}, {input})",
    input=[q.truth_d4_p4],
    output=[q.truth_d4_eta],
    scopes=["mm", "ee", "em"],
)

truth_d4_phi = Producer(
    name="truth_d4_phi",
    call="quantities::phi({df}, {output}, {input})",
    input=[q.truth_d4_p4],
    output=[q.truth_d4_phi],
    scopes=["mm", "ee", "em"],
)

truth_d4Quantities = ProducerGroup(
	name = "truth_d4Quantities",
	call = None,
	input = None,
	output = None,
	scopes = ["mm", "ee", "em"],
	subproducers = [truth_d4_pt, truth_d4_eta, truth_d4_phi],
)

truth_H_pt = Producer(
    name="truth_H_pt",
    call="quantities::pt({df}, {output}, {input})",
    input=[q.truth_H_p4],
    output=[q.truth_H_pt],
    scopes=["mm", "ee", "em"],
)

truth_H_eta = Producer(
    name="truth_H_eta",
    call="quantities::eta({df}, {output}, {input})",
    input=[q.truth_H_p4],
    output=[q.truth_H_eta],
    scopes=["mm", "ee", "em"],
)

truth_H_phi = Producer(
    name="truth_H_phi",
    call="quantities::phi({df}, {output}, {input})",
    input=[q.truth_H_p4],
    output=[q.truth_H_phi],
    scopes=["mm", "ee", "em"],
)
truth_H_mass = Producer(
    name="truth_H_mass",
    call="quantities::mass({df}, {output}, {input})",
    input=[q.truth_H_p4],
    output=[q.truth_H_mass],
    scopes=["mm", "ee", "em"],
)

truth_hQuantities = ProducerGroup(
	name = "truth_hQuantities",
	call = None,
	input = None,
	output = None,
	scopes = ["mm", "ee", "em"],
	subproducers = [truth_H_pt, truth_H_eta, truth_H_phi, truth_H_mass],
)


truth_dQuantities = ProducerGroup(
	name = "truth_dQuantities",
	call = None,
	input = None,
	output = None,
	scopes = ["mm", "ee", "em"],
	subproducers = [truth_d1Quantities, truth_d2Quantities, truth_d3Quantities, truth_d4Quantities, truth_hQuantities],
)
