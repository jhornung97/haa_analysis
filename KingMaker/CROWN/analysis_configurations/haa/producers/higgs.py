from ..quantities import output as q
from ..quantities import nanoAOD as nanoAOD
from code_generation.producer import Producer, Filter, ProducerGroup

HiggsSelection = Producer(
	name = "HiggsSelection",
	call = "haa::ClosestToHiggsMassAlgo({df}, {input}, {output})",
	input = [
            nanoAOD.PFCands_pt,
            nanoAOD.PFCands_eta,
            nanoAOD.PFCands_phi,
            nanoAOD.PFCands_mass,
            nanoAOD.PFCands_charge,
            q.base_pfcands_mask,
        ],
	output = [q.higgsdaughters],
	scopes = ["mm", "ee", "em"],
)

FourHardestPFCands = Producer(
        name = "FourHardestPFCands",
        call = "haa::FourHardestPFCandsAlgo({df}, {input}, {output})",
        input = [
            nanoAOD.PFCands_pt,
            q.base_pfcands_mask,
        ],
        output = [q.higgsdaughters],
        scopes = ["mm", "ee", "em"],
)
ChargePairs = Producer(
        name = "ChargePairs",
        call = "haa::ChargePairsAlgo({df}, {input}, {output})",
        input = [
            nanoAOD.PFCands_pt,
            nanoAOD.PFCands_eta,
            nanoAOD.PFCands_phi,
            nanoAOD.PFCands_mass,
            nanoAOD.PFCands_charge,
            q.base_pfcands_mask,
        ],
        output = [q.higgsdaughters],
        scopes = ["mm", "ee", "em"],
)
GoodHiggsDaughtersFlag = Producer(
    name="GoodHiggsDaughtersFlag",
    call="ditau_pairselection::flagGoodPairs({df}, {output}, {input})",
    input=[q.higgsdaughters],
    output=[],
    scopes=["mm", "ee", "em"],
)

GoodHiggsDaughtersFilter = Filter(
    name="GoodHiggsDaughtersFilter",
    call='basefunctions::FilterFlagsAny({df}, "GoodHiggsDaughters", {input})',
    input=[],
    scopes=["mm", "ee", "em"],
    subproducers=[GoodHiggsDaughtersFlag],
)

H_p4 = Producer(
    name="H_p4",
    call="haa::GetHiggsP4({df}, {input}, {output})",
    input=[
        q.d1_p4,
        q.d2_p4,
        q.d3_p4,
        q.d4_p4,
    ],
    output=[q.H_p4],
    scopes=["mm", "ee", "em"],
)

H_pt = Producer(
    name="H_pt",
    call="quantities::pt({df}, {output}, {input})",
    input=[q.H_p4],
    output=[q.H_pt],
    scopes=["mm", "ee", "em"],
)

H_eta = Producer(
    name="H_eta",
    call="quantities::eta({df}, {output}, {input})",
    input=[q.H_p4],
    output=[q.H_eta],
    scopes=["mm", "ee", "em"],
)

H_phi = Producer(
    name="H_phi",
    call="quantities::phi({df}, {output}, {input})",
    input=[q.H_p4],
    output=[q.H_phi],
    scopes=["mm", "ee", "em"],
)

H_mass = Producer(
    name="H_mass",
    call="quantities::mass({df}, {output}, {input})",
    input=[q.H_p4],
    output=[q.H_mass],
    scopes=["mm", "ee", "em"],
)

HiggsQuantities = ProducerGroup(
	name = "HiggsQuantities",
	call = None,
	input = None,
	output = None,
	scopes = ["mm", "ee", "em"],
	subproducers = [H_p4, H_pt, H_eta, H_phi, H_mass],
)

GetPS = Producer(
    name="GetPS",
    call="haa::GetMinMassDiff({df}, {input}, {output})",
    input=[
        q.higgsdaughters,
        nanoAOD.PFCands_pt,
        nanoAOD.PFCands_eta,
        nanoAOD.PFCands_phi,
        nanoAOD.PFCands_mass,
        nanoAOD.PFCands_charge,
        q.base_pfcands_mask,
    ],
    output=[q.ps1Pair, q.ps2Pair],
    scopes=["mm", "ee", "em"],
)

ps_1_mass = Producer(
    name="ps_1_mass",
    call="quantities::m_vis({df}, {output}, {input_vec})",
    input=[q.ps_1_d_1_p4, q.ps_1_d_2_p4],
    output=[q.ps_1_mass],
    scopes=["mm", "ee", "em"],
)

ps_1_deltaR = Producer(
    name="ps_1_deltaR",
    call="haa::getDaughterDeltaR({df}, {input_vec}, {output})",
    input=[q.ps_1_d_1_p4, q.ps_1_d_2_p4],
    output=[q.ps_1_deltaR],
    scopes=["mm", "ee", "em"],
)

ps_2_mass = Producer(
    name="ps_2_mass",
    call="quantities::m_vis({df}, {output}, {input_vec})",
    input=[q.ps_2_d_1_p4, q.ps_2_d_2_p4],
    output=[q.ps_2_mass],
    scopes=["mm", "ee", "em"],
)

ps_2_deltaR = Producer(
    name="ps_2_deltaR",
    call="haa::getDaughterDeltaR({df}, {input_vec}, {output})",
    input=[q.ps_2_d_1_p4, q.ps_2_d_2_p4],
    output=[q.ps_2_deltaR],
    scopes=["mm", "ee", "em"],
)

ps_quantities = ProducerGroup(
    name="ps_quantities",
    call=None,
    input=None,
    output=None,
    scopes=["mm", "ee", "em"],
    subproducers=[ps_1_mass, ps_2_mass, ps_1_deltaR, ps_2_deltaR]
)

GetTruthPS = Producer(
    name="GetTruthPS",
    call="haa::GetTruthDaughterPairs({df}, {input}, {output})",
    input=[
        q.truth_daughters,
        nanoAOD.GenParticle_pdgId,
        nanoAOD.GenParticle_motherid,
    ],
    output=[q.truth_ps_1, q.truth_ps_2],
    scopes=["mm", "ee", "em"],
)

Truth_PS_1_D_1_P4 = Producer(
    name="Truth_PS_1_D_1_P4",
    call="lorentzvectors::build({df}, {input_vec}, 0, {output})",
    input=[
        q.truth_ps_1,
        nanoAOD.GenParticle_pt,
        nanoAOD.GenParticle_eta,
        nanoAOD.GenParticle_phi,
        nanoAOD.GenParticle_mass,
    ],
    output=[q.truth_ps_1_d_1_p4],
    scopes=["mm", "ee", "em"],
)

Truth_PS_1_D_2_P4 = Producer(
    name="Truth_PS_1_D_2_P4",
    call="lorentzvectors::build({df}, {input_vec}, 1, {output})",
    input=[
        q.truth_ps_1,
        nanoAOD.GenParticle_pt,
        nanoAOD.GenParticle_eta,
        nanoAOD.GenParticle_phi,
        nanoAOD.GenParticle_mass,
    ],
    output=[q.truth_ps_1_d_2_p4],
    scopes=["mm", "ee", "em"],
)

truth_ps_1_mass = Producer(
    name="truth_ps_1_mass",
    call="quantities::m_vis({df}, {output}, {input_vec})",
    input=[q.truth_ps_1_d_1_p4, q.truth_ps_1_d_2_p4],
    output=[q.truth_ps_1_mass],
    scopes=["mm", "ee", "em"],
)

truth_ps_1_deltaR = Producer(
    name="truth_ps_1_deltaR",
    call="haa::getDaughterDeltaR({df}, {input_vec}, {output})",
    input=[q.truth_ps_1_d_1_p4, q.truth_ps_1_d_2_p4],
    output=[q.truth_ps_1_deltaR],
    scopes=["mm", "ee", "em"],
)

Truth_PS_2_D_1_P4 = Producer(
    name="Truth_PS_2_D_1_P4",
    call="lorentzvectors::build({df}, {input_vec}, 0, {output})",
    input=[
        q.truth_ps_2,
        nanoAOD.GenParticle_pt,
        nanoAOD.GenParticle_eta,
        nanoAOD.GenParticle_phi,
        nanoAOD.GenParticle_mass,
    ],
    output=[q.truth_ps_2_d_1_p4],
    scopes=["mm", "ee", "em"],
)
Truth_PS_2_D_2_P4 = Producer(
    name="Truth_PS_2_D_2_P4",
    call="lorentzvectors::build({df}, {input_vec}, 1, {output})",
    input=[
        q.truth_ps_2,
        nanoAOD.GenParticle_pt,
        nanoAOD.GenParticle_eta,
        nanoAOD.GenParticle_phi,
        nanoAOD.GenParticle_mass,
    ],
    output=[q.truth_ps_2_d_2_p4],
    scopes=["mm", "ee", "em"],
)
truth_ps_2_mass = Producer(
    name="truth_ps_2_mass",
    call="quantities::m_vis({df}, {output}, {input_vec})",
    input=[q.truth_ps_2_d_1_p4, q.truth_ps_2_d_2_p4],
    output=[q.truth_ps_2_mass],
    scopes=["mm", "ee", "em"],
)
truth_ps_2_deltaR = Producer(
    name="truth_ps_2_deltaR",
    call="haa::getDaughterDeltaR({df}, {input_vec}, {output})",
    input=[q.truth_ps_2_d_1_p4, q.truth_ps_2_d_2_p4],
    output=[q.truth_ps_2_deltaR],
    scopes=["mm", "ee", "em"],
)
TruthPSQuantities = ProducerGroup(
    name="TruthPSQuantities",
    call=None,
    input=None,
    output=None,
    scopes=["mm", "ee", "em"],
    subproducers=[
        Truth_PS_1_D_1_P4,
        Truth_PS_1_D_2_P4,
        truth_ps_1_mass,
        truth_ps_1_deltaR,
        Truth_PS_2_D_1_P4,
        Truth_PS_2_D_2_P4,
        truth_ps_2_mass,
        truth_ps_2_deltaR,
    ],
)
'''
ps_mass_12 = Producer(
    name="ps_mass_12",
    call="quantities::m_vis({df}, {output}, {input_vec})",
    input=[q.d1_p4, q.d2_p4],
    output=[q.ps_mass_12],
    scopes=["mm", "ee", "em"],
)
ps_mass_14 = Producer(
    name="ps_mass_14",
    call="quantities::m_vis({df}, {output}, {input_vec})",
    input=[q.d1_p4, q.d4_p4],
    output=[q.ps_mass_14],
    scopes=["mm", "ee", "em"],
)
ps_mass_23 = Producer(
    name="ps_mass_23",
    call="quantities::m_vis({df}, {output}, {input_vec})",
    input=[q.d2_p4, q.d3_p4],
    output=[q.ps_mass_23],
    scopes=["mm", "ee", "em"],
)
ps_mass_34 = Producer(
    name="ps_mass_34",
    call="quantities::m_vis({df}, {output}, {input_vec})",
    input=[q.d3_p4, q.d4_p4],
    output=[q.ps_mass_34],
    scopes=["mm", "ee", "em"],
)

ps_quanities = ProducerGroup(
    name="ps_quanities",
    call=None,
    input=None,
    output=None,
    scopes=["mm", "ee", "em"],
    subproducers=[ps_mass_12, ps_mass_14, ps_mass_23, ps_mass_34],
)
'''