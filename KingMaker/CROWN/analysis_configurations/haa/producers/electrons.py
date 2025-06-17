from ..quantities import output as q
from ..quantities import nanoAOD as nanoAOD
from code_generation.producer import Producer, ProducerGroup

####################
# Set of producers used for loosest selection of electrons
####################

ElectronPtCut = Producer(
    name="ElectronPtCut",
    call="physicsobject::CutPt({df}, {input}, {output}, {min_electron_pt})",
    input=[nanoAOD.Electron_pt],
    output=[],
    scopes=["global"],
)
ElectronEtaCut = Producer(
    name="ElectronEtaCut",
    call="physicsobject::CutEta({df}, {input}, {output}, {max_electron_eta})",
    input=[nanoAOD.Electron_eta],
    output=[],
    scopes=["global"],
)
ElectronDxyCut = Producer(
    name="ElectronDxyCut",
    call="physicsobject::CutDxy({df}, {input}, {output}, {max_electron_dxy})",
    input=[nanoAOD.Electron_dxy],
    output=[],
    scopes=["global"],
)
ElectronDzCut = Producer(
    name="ElectronDzCut",
    call="physicsobject::CutDz({df}, {input}, {output}, {max_electron_dz})",
    input=[nanoAOD.Electron_dz],
    output=[],
    scopes=["global"],
)

ElectronIDCut = Producer(
    name="ElectronIDCut",
    call='physicsobject::electron::CutCBID({df}, {output}, "{electron_id}", {electron_id_wp})',
    input=[],
    output=[],
    scopes=["global"],
)

ElectronIsoCut = Producer(
    name="ElectronIsoCut",
    call="physicsobject::electron::CutIsolation({df}, {output}, {input}, {electron_iso_cut})",
    input=[nanoAOD.Electron_iso],
    output=[],
    scopes=["global"],
)

BaseElectrons = ProducerGroup(
    name="BaseElectrons",
    call="physicsobject::CombineMasks({df}, {output}, {input})",
    input=[],
    output=[q.base_electrons_mask],
    scopes=["global"],
    subproducers=[
        #ElectronPtCut,
        #ElectronEtaCut,
        #ElectronDxyCut,
        #ElectronDzCut,
        ElectronIDCut,
        #ElectronIsoCut,
    ],
)

####################
# Set of producers used for more specific selection of electrons in channels
####################

GoodElectronPtCut = Producer(
    name="GoodElectronPtCut",
    call="physicsobject::CutPt({df}, {input}, {output}, {min_electron_pt})",
    input=[nanoAOD.Electron_pt],
    output=[],
    scopes=["ee","em"],
)
GoodElectronEtaCut = Producer(
    name="GoodElectronEtaCut",
    call="physicsobject::CutEta({df}, {input}, {output}, {max_electron_eta})",
    input=[nanoAOD.Electron_eta],
    output=[],
    scopes=["ee","em"],
)
GoodElectronDxyCut = Producer(
    name="GoodElectronDxyCut",
    call="physicsobject::CutDxy({df}, {input}, {output}, {max_electron_dxy})",
    input=[nanoAOD.Electron_dxy],
    output=[],
    scopes=["ee","em"],
)
GoodElectronDzCut = Producer(
    name="GoodElectronDzCut",
    call="physicsobject::CutDz({df}, {input}, {output}, {max_electron_dz})",
    input=[nanoAOD.Electron_dz],
    output=[],
    scopes=["ee","em"],
)
GoodElectronIDCut = Producer(
    name="GoodElectronIDCut",
    call='physicsobject::electron::CutCBID({df}, {output}, "{electron_id}", {electron_id_wp})',
    input=[],
    output=[],
    scopes=["ee","em"],
)
GoodElectronIsoCut = Producer(
    name="GoodElectronIsoCut",
    call="physicsobject::electron::CutIsolation({df}, {output}, {input}, {electron_iso_cut})",
    input=[nanoAOD.Electron_iso],
    output=[],
    scopes=["ee","em"],
)
GoodElectrons = ProducerGroup(
    name="GoodElectrons",
    call="physicsobject::CombineMasks({df}, {output}, {input})",
    input=[q.base_electrons_mask],
    output=[q.good_electrons_mask],
    scopes=["ee","em"],
    subproducers=[
        GoodElectronPtCut,
        GoodElectronEtaCut,
        GoodElectronDxyCut,
        GoodElectronDzCut,
        GoodElectronIDCut,
        GoodElectronIsoCut,
    ],
)
NumberOfGoodElectrons = Producer(
    name="NumberOfGoodElectrons",
    call="quantities::NumberOfGoodLeptons({df}, {output}, {input})",
    input=[q.good_electrons_mask],
    output=[q.nelectrons],
    scopes=["ee","em"],
)
VetoElectrons = Producer(
    name="VetoElectrons",
    call="physicsobject::VetoCandInMask({df}, {output}, {input}, {electron_index_in_pair})",
    input=[q.base_electrons_mask, q.dileptonpair],
    output=[q.veto_electrons_mask],
    scopes=["ee"],
)
VetoSecondElectron = Producer(
    name="VetoSecondElectron",
    call="physicsobject::VetoCandInMask({df}, {output}, {input}, {second_electron_index_in_pair})",
    input=[q.veto_electrons_mask, q.dileptonpair],
    output=[q.veto_electrons_mask_2],
    scopes=["ee"],
)

ExtraElectronsVeto = Producer(
    name="ExtraElectronsVeto",
    call="physicsobject::LeptonVetoFlag({df}, {output}, {input})",
    input={
        "ee": [q.veto_electrons_mask_2],
    },
    output=[q.electron_veto_flag],
    scopes=["ee"],
)

####################
# Set of producers used for di-electron veto
####################
'''
DiElectronVetoPtCut = Producer(
    name="DiElectronVetoPtCut",
    call="physicsobject::CutPt({df}, {input}, {output}, {min_dielectronveto_pt})",
    input=[nanoAOD.Electron_pt],
    output=[],
    scopes=["global"],
)
DiElectronVetoIDCut = Producer(
    name="DiElectronVetoIDCut",
    call='physicsobject::electron::CutID({df}, {output}, "{dielectronveto_id}")',
    input=[],
    output=[],
    scopes=["global"],
)
DiElectronVetoElectrons = ProducerGroup(
    name="DiElectronVetoElectrons",
    call="physicsobject::CombineMasks({df}, {output}, {input})",
    input=ElectronEtaCut.output + ElectronDxyCut.output + ElectronDzCut.output + ElectronIsoCut.output,
    output=[],
    scopes=["global"],
    subproducers=[
        DiElectronVetoPtCut,
        DiElectronVetoIDCut,
    ],
)
DiElectronVeto = ProducerGroup(
    name="DiElectronVeto",
    call="physicsobject::CheckForDiLeptonPairs({df}, {output}, {input}, {dileptonveto_dR})",
    input=[
        nanoAOD.Electron_pt,
        nanoAOD.Electron_eta,
        nanoAOD.Electron_phi,
        nanoAOD.Electron_mass,
        nanoAOD.Electron_charge,
    ],
    output=[q.dielectron_veto],
    scopes=["global"],
    subproducers=[DiElectronVetoElectrons],
)
'''
