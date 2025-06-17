from ..quantities import output as q
from ..quantities import nanoAOD as nanoAOD
from code_generation.producer import Producer, ProducerGroup

PFCandsPtCut = Producer(
    name="PFCandsPtCut",
    call="physicsobject::CutPt({df}, {input}, {output}, {min_jet_pt})",
    input=[nanoAOD.PFCands_pt],
    output=[],
    scopes=["global"],
)

GoodPFCands = ProducerGroup(
    name="GoodPFCands",
    call="physicsobject::CombineMasks({df}, {output}, {input})",
    input=[],
    output=[q.good_pfcands_mask],
    scopes=["global"],
    subproducers=[PFCandsPtCut],
)

NumberOfPFCands = Producer(
    name="NumberOfPFCands",
    call="quantities::jet::NumberOfJets({df}, {output}, {input})",
    input=[q.good_pfcands_mask],
    output=[q.npfcands],
    scopes=["global"],
)
