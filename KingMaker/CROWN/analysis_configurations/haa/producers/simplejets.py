from ..quantities import output as q
from ..quantities import nanoAOD as nanoAOD
from code_generation.producer import Producer, ProducerGroup

JetPtCut = Producer(
    name="JetPtCut",
    call="physicsobject::CutPt({df}, {input}, {output}, {min_jet_pt})",
    input=[nanoAOD.Jet_pt],
    output=[],
    scopes=["global"],
)

GoodJets = ProducerGroup(
    name="GoodJets",
    call="physicsobject::CombineMasks({df}, {output}, {input})",
    input=[],
    output=[q.good_jets_mask],
    scopes=["global"],
    subproducers=[JetPtCut],
)

NumberOfJets = Producer(
    name="NumberOfJets",
    call="quantities::jet::NumberOfJets({df}, {output}, {input})",
    input=[q.good_jets_mask],
    output=[q.njets],
    scopes=["global"],
)
