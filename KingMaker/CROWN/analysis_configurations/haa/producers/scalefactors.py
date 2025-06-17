from ..quantities import output as q
from code_generation.producer import Producer, ProducerGroup
from code_generation.producer import ExtendedVectorProducer


############################
# Muon ID, ISO SF
# The readout is done via correctionlib
############################

Muon_1_ID_SF_RooWorkspace = Producer(
    name="MuonID_SF_RooWorkspace",
    call='scalefactor::muon::id_rooworkspace({df}, {input}, {output}, "{muon_sf_workspace}", "{muon_sf_id_name}", "{muon_sf_id_args}")',
    input=[q.pt_1, q.eta_1],
    output=[q.id_wgt_mu_1],
    scopes=["mt", "mm"],
)
Muon_1_Iso_SF_RooWorkspace = Producer(
    name="MuonIso_SF_RooWorkspace",
    call='scalefactor::muon::iso_rooworkspace({df}, {input}, {output}, "{muon_sf_workspace}", "{muon_sf_iso_name}", "{muon_sf_iso_args}")',
    input=[q.pt_1, q.eta_1, q.iso_1],
    output=[q.iso_wgt_mu_1],
    scopes=["mt", "mm"],
)
Muon_1_ID_SF = Producer(
    name="MuonID_SF",
    call='scalefactor::muon::id({df}, {input}, "{muon_sf_year_id}", "{muon_sf_varation}", {output}, "{muon_sf_file}", "{muon_id_sf_name}")',
    input=[q.pt_1, q.eta_1],
    output=[q.id_wgt_mu_1],
    scopes=["mm"],
)
Muon_1_Iso_SF = Producer(
    name="MuonIso_SF",
    call='scalefactor::muon::iso({df}, {input}, "{muon_sf_year_id}", "{muon_sf_varation}", {output}, "{muon_sf_file}", "{muon_iso_sf_name}")',
    input=[q.pt_1, q.eta_1],
    output=[q.iso_wgt_mu_1],
    scopes=["mm"],
)
Muon_1_Trigger_SF = Producer(
    name="MuonTrigger_SF",
    call='scalefactor::muon::trigger({df}, {input}, "{muon_sf_year_id}", "{muon_sf_varation}", {output}, "{muon_sf_file}", "{muon_trigger_sf_name}")',
    input=[q.pt_1, q.eta_1],
    output=[q.trigger_wgt_mu_1],
    scopes=["mm"],
)
Muon_2_Trigger_SF = Producer(
    name="MuonTrigger_SF",
    call='scalefactor::muon::trigger({df}, {input}, "{muon_sf_year_id}", "{muon_sf_varation}", {output}, "{muon_sf_file}", "{muon_trigger_sf_name}")',
    input=[q.pt_2, q.eta_2],
    output=[q.trigger_wgt_mu_2],
    scopes=["mm","em"],
)
Muon_2_ID_SF_RooWorkspace = Producer(
    name="MuonID_SF_RooWorkspace",
    call='scalefactor::muon::id_rooworkspace({df}, {input}, {output}, "{muon_sf_workspace}", "{muon_sf_id_name}", "{muon_sf_id_args}")',
    input=[q.pt_2, q.eta_2],
    output=[q.id_wgt_mu_2],
    scopes=["em", "mm"],
)
Muon_2_Iso_SF_RooWorkspace = Producer(
    name="MuonIso_SF_RooWorkspace",
    call='scalefactor::muon::iso_rooworkspace({df}, {input}, {output}, "{muon_sf_workspace}", "{muon_sf_iso_name}", "{muon_sf_iso_args}")',
    input=[q.pt_2, q.eta_2, q.iso_2],
    output=[q.iso_wgt_mu_2],
    scopes=["em", "mm"],
)
Muon_2_ID_SF = Producer(
    name="MuonID_SF",
    call='scalefactor::muon::id({df}, {input}, "{muon_sf_year_id}", "{muon_sf_varation}", {output}, "{muon_sf_file}", "{muon_id_sf_name}")',
    input=[q.pt_2, q.eta_2],
    output=[q.id_wgt_mu_2],
    scopes=["em", "mm"],
)
Muon_2_Iso_SF = Producer(
    name="MuonIso_SF",
    call='scalefactor::muon::iso({df}, {input}, "{muon_sf_year_id}", "{muon_sf_varation}", {output}, "{muon_sf_file}", "{muon_iso_sf_name}")',
    input=[q.pt_2, q.eta_2],
    output=[q.iso_wgt_mu_2],
    scopes=["em", "mm"],
)

Ele_1_ID_SF = Producer(
    name="Ele_1_ID_SF",
    call='scalefactor::electron::id({df}, {input}, "{ele_sf_year_id}", "Medium", "{ele_sf_varation}", {output}, "{ele_sf_file}", "{ele_id_sf_name}")',
    input=[q.pt_1, q.eta_1],
    output=[q.id_wgt_ele_1],
    scopes=["em", "ee"],
)
Ele_2_ID_SF = Producer(
    name="Ele_2_ID_SF",
    call='scalefactor::electron::id({df}, {input}, "{ele_sf_year_id}", "Medium", "{ele_sf_varation}", {output}, "{ele_sf_file}", "{ele_id_sf_name}")',
    input=[q.pt_2, q.eta_2],
    output=[q.id_wgt_ele_2],
    scopes=["ee"],
)

Ele_1_Trigger_SF = Producer(
    name="Ele_1_Trigger_SF",
    call='scalefactor::electron::trigger({df}, {input}, {output}, "{nom_ele_trigger_sf_file}", "{nom_ele_trigger_sf_name}", "{syst_ele_trigger_sf_file}", "{syst_ele_trigger_sf_name}")',
    input=[q.pt_1, q.eta_1],
    output=[q.trigger_wgt_ele_1, q.trigger_wgt_ele_1_up, q.trigger_wgt_ele_1_down],
    scopes=["ee","em"],
)

Ele_2_Trigger_SF = Producer(
    name="Ele_2_Trigger_SF",
    call='scalefactor::electron::trigger({df}, {input}, {output}, "{nom_ele_trigger_sf_file}", "{nom_ele_trigger_sf_name}", "{syst_ele_trigger_sf_file}", "{syst_ele_trigger_sf_name}")',
    input=[q.pt_2, q.eta_2],
    output=[q.trigger_wgt_ele_2, q.trigger_wgt_ele_2_up, q.trigger_wgt_ele_2_down],
    scopes=["ee","em"], 
)

MuonIDIsoTrigger_SF = ProducerGroup(
    name="MuonIDIso_SF",
    call=None,
    input=None,
    output=None,
    scopes=["em", "mm"],
    subproducers={
        "em": [
            Muon_2_ID_SF,
            Muon_2_Iso_SF,
            Muon_2_Trigger_SF,
        ],
        "mm": [
            Muon_1_ID_SF,
            Muon_1_Iso_SF,
            Muon_1_Trigger_SF,
            Muon_2_ID_SF,
            Muon_2_Iso_SF,
            Muon_2_Trigger_SF,
        ],
    },
)

ElectronIDTrigger_SF = ProducerGroup(
    name="ElectronIDTrigger_SF",
    call=None,
    input=None,
    output=None,
    scopes=["ee", "em"],
    subproducers={
        "ee": [
            Ele_1_ID_SF,
            Ele_2_ID_SF,
            Ele_1_Trigger_SF,
            Ele_2_Trigger_SF
        ],
        "em": [
            Ele_1_ID_SF,
            Ele_1_Trigger_SF,
        ]
    }
)


MuonIDIso_SF_RooWorkspace = ProducerGroup(
    name="MuonIDIso_SF_RooWorkspace",
    call=None,
    input=None,
    output=None,
    scopes=["mt", "em", "mm"],
    subproducers={
        "mt": [
            Muon_1_ID_SF_RooWorkspace,
            Muon_1_Iso_SF_RooWorkspace,
        ],
        "em": [
            Muon_2_ID_SF_RooWorkspace,
            Muon_2_Iso_SF_RooWorkspace,
        ],
        "mm": [
            Muon_1_ID_SF_RooWorkspace,
            Muon_1_Iso_SF_RooWorkspace,
            Muon_2_ID_SF_RooWorkspace,
            Muon_2_Iso_SF_RooWorkspace,
        ],
    },
)


