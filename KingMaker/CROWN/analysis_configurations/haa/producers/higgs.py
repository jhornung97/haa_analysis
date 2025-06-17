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

ps_2_mass = Producer(
    name="ps_2_mass",
    call="quantities::m_vis({df}, {output}, {input_vec})",
    input=[q.ps_2_d_1_p4, q.ps_2_d_2_p4],
    output=[q.ps_2_mass],
    scopes=["mm", "ee", "em"],
)

ps_masses = ProducerGroup(
    name="ps_masses",
    call=None,
    input=None,
    output=None,
    scopes=["mm", "ee", "em"],
    subproducers=[ps_1_mass, ps_2_mass],
)

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