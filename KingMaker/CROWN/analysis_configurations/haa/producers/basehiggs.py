from ..quantities import output as q
from ..quantities import nanoAOD as nanoAOD
from code_generation.producer import Producer, ProducerGroup

HiggsSelection = Producer(
	name = "HiggsSelection",
	call = "haa::HiggsSelection({df}, {input}, {output})",
	input = [
			 nanoAOD.nPFCands, 
			 nanoAOD.PFCands_pt, 
			 nanoAOD.PFCands_eta, 
			 nanoAOD.PFCands_phi, 
			 nanoAOD.PFCands_mass,
			 nanoAOD.PFCands_charge, 
			 nanoAOD.PFCands_pdgId,
			],
	output = [q.H_p4],
	scopes = ["global"],
)

H_pt = Producer(
    name="H_pt",
    call="quantities::pt({df}, {output}, {input})",
    input=[q.H_p4],
    output=[q.H_pt],
    scopes=["global"],
)

H_eta = Producer(
    name="H_eta",
    call="quantities::eta({df}, {output}, {input})",
    input=[q.H_p4],
    output=[q.H_eta],
    scopes=["global"],
)

H_phi = Producer(
    name="H_phi",
    call="quantities::phi({df}, {output}, {input})",
    input=[q.H_p4],
    output=[q.H_phi],
    scopes=["global"],
)

H_mass = Producer(
    name="H_mass",
    call="quantities::mass({df}, {output}, {input})",
    input=[q.H_p4],
    output=[q.H_mass],
    scopes=["global"],
)

HiggsQuantities = ProducerGroup(
	name = "HiggsQuantities",
	call = None,
	input = None,
	output = None,
	scopes = ["global"],
	subproducers = [H_pt, H_eta, H_phi, H_mass],
)
