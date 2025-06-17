from ..quantities import output as q
from ..quantities import nanoAOD as nanoAOD
from code_generation.producer import Producer, ProducerGroup

####################                                                                                     
# Set of producers used for loosest selection of photons                                               
#################### 

PhotonIDCut = Producer(
    name="PhotonIDCut",
    call='physicsobject::electron::CutCBID({df}, {output}, "{photon_id}", {photon_id_wp})',
    input=[],
    output=[],
    scopes=["global"],
)
BasePhotons = ProducerGroup(
    name="BasePhotons",
    call="physicsobject::CombineMasks({df}, {output}, {input})",
    input=[],
    output=[q.base_photons_mask],
    scopes=["global"],
    subproducers=[
        PhotonIDCut,
    ],
)
