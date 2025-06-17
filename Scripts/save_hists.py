import ROOT
import sys

year = sys.argv[1]
scope = sys.argv[2]

luminosities = {
    "2016preVFP": 19.5e3,
    "2016postVFP": 16.8e3,
    "2017": 41.48e3,
    "2018": 59.83e3
}

lumi = luminosities.get(year, 0)

def define_weights(df, scope, lumi):
    nominal_weights = {
        "mm": [
            "id_wgt_mu_1", "id_wgt_mu_2", "iso_wgt_mu_1", "iso_wgt_mu_2",
            "trigger_wgt_mu_1", "trigger_wgt_mu_2", "puweight"
        ],
        "ee": [
            "id_wgt_ele_1", "id_wgt_ele_2", "trigger_wgt_ele_1", "trigger_wgt_ele_2", "puweight"
        ]
    }

    systematics = {
        "pu": ["puweight_up", "puweight_down"],
        "id_1": ["id_wgt_mu_1__MuonIDUp", "id_wgt_mu_1__MuonIDDown"] if scope == "mm" else ["id_wgt_ele_1__ElectronIDUp", "id_wgt_ele_1__ElectronIDDown"],
        "id_2": ["id_wgt_mu_2__MuonIDUp", "id_wgt_mu_2__MuonIDDown"] if scope == "mm" else ["id_wgt_ele_2__ElectronIDUp", "id_wgt_ele_2__ElectronIDDown"],
        "iso_1": ["iso_wgt_mu_1__MuonIsoUp", "iso_wgt_mu_1__MuonIsoDown"] if scope == "mm" else ["1", "1"],
        "iso_2": ["iso_wgt_mu_2__MuonIsoUp", "iso_wgt_mu_2__MuonIsoDown"] if scope == "mm" else ["1", "1"],
        "trigger_1": ["trigger_wgt_mu_1__MuonTriggerUp", "trigger_wgt_mu_1__MuonTriggerDown"] if scope == "mm" else ["trigger_wgt_ele_1_up", "trigger_wgt_ele_1_down"],
        "trigger_2": ["trigger_wgt_mu_2__MuonTriggerUp", "trigger_wgt_mu_2__MuonTriggerDown"] if scope == "mm" else ["trigger_wgt_ele_2_up", "trigger_wgt_ele_2_down"]
    }

    weight_expr = f"{lumi}*signs*evtweight*" + "*".join(nominal_weights[scope])
    print(weight_expr)
    df = df.Define("nominal", weight_expr)

    for var, vals in systematics.items():
        if scope == "ee" and "iso" in var:
            continue

        if scope == "ee" and "trigger_1" in var:
            df = df.Define(f"weight_{var}_up", weight_expr.replace("trigger_wgt_ele_1", "trigger_wgt_ele_1_up"))
            df = df.Define(f"weight_{var}_down", weight_expr.replace("trigger_wgt_ele_1", "trigger_wgt_ele_1_down"))

            print(weight_expr.replace("trigger_wgt_ele_1", "trigger_wgt_ele_1_up"))
            print(weight_expr.replace("trigger_wgt_ele_1", "trigger_wgt_ele_1_down"))
        elif scope == "ee" and "trigger_2" in var:
            df = df.Define(f"weight_{var}_up", weight_expr.replace("trigger_wgt_ele_2", "trigger_wgt_ele_2_up"))
            df = df.Define(f"weight_{var}_down", weight_expr.replace("trigger_wgt_ele_2", "trigger_wgt_ele_2_down"))
            
            print(weight_expr.replace("trigger_wgt_ele_2", "trigger_wgt_ele_2_up"))
            print(weight_expr.replace("trigger_wgt_ele_2", "trigger_wgt_ele_2_down"))
        else:
            nominal_weight = vals[0].split("__")[0] if "__" in vals[0] else vals[0].split("_")[0]
        
            df = df.Define(f"weight_{var}_up", weight_expr.replace(nominal_weight, vals[0]))
            df = df.Define(f"weight_{var}_down", weight_expr.replace(nominal_weight, vals[1]))
        
            print(weight_expr.replace(nominal_weight, vals[0]))
            print(weight_expr.replace(nominal_weight, vals[1]))

    return df

if scope in ["mm", "ee"]:
    file_path = f"/ceph/jhornung/MC_2018/2018/HaaKKKK_M125_13TeV_powheg_pythia8-RunIISummer20UL18MiniAODv2-106X/{scope}/HaaKKKK_M125_13TeV_powheg_pythia8-RunIISummer20UL18MiniAODv2-106X.root"
    df = ROOT.RDataFrame("ntuple", file_path)
    df = df.Define(f"H_mass_{year}_{scope}", "H_mass")
    df = df.Define("signs", "genWeight > 0 ? 1 : -1")
    df = define_weights(df, scope, lumi)
    df = df.Define("ps_mass_mask", "ps_1_mass < 3 && ps_2_mass < 3 && abs(ps_1_mass-ps_2_mass) < 1")
    df = df.Define("iso_mask", "d1_iso < 3 && d2_iso < 3 && d3_iso < 3 && d4_iso < 3")
    df = df.Define("mask", "pt_vis > 30 && (m_vis > 75 && m_vis < 105) && abs(H_eta) < 2.4 && d1_pt > 10 && d2_pt > 10 && ps_mass_mask && iso_mask")

df = df.Filter("mask")

with ROOT.TFile(f"/work/jhornung/Haa/limits/{year}/{scope}/signal_hists.root", "RECREATE") as f:
    hists = {
        "nominal": df.Histo1D(("Haa", "Haa", 65, 70, 200), f"H_mass_{year}_{scope}", "nominal"),
        "pu_up": df.Histo1D(("Haa_puUp", "Haa_puUp", 65, 70, 200), f"H_mass_{year}_{scope}", "weight_pu_up"),
        "pu_down": df.Histo1D(("Haa_puDown", "Haa_puDown", 65, 70, 200), f"H_mass_{year}_{scope}", "weight_pu_down"),
        "id_1_up": df.Histo1D(("Haa_id_1Up", "Haa_id_1Up", 65, 70, 200), f"H_mass_{year}_{scope}", "weight_id_1_up"),
        "id_1_down": df.Histo1D(("Haa_id_1Down", "Haa_id_1Down", 65, 70, 200), f"H_mass_{year}_{scope}", "weight_id_1_down"),
        "id_2_up": df.Histo1D(("Haa_id_2Up", "Haa_id_2Up", 65, 70, 200), f"H_mass_{year}_{scope}", "weight_id_2_up"),
        "id_2_down": df.Histo1D(("Haa_id_2Down", "Haa_id_2Down", 65, 70, 200), f"H_mass_{year}_{scope}", "weight_id_2_down"),
        "trigger_1_up": df.Histo1D(("Haa_trigger_1Up", "Haa_trigger_1Up", 65, 70, 200), f"H_mass_{year}_{scope}", "weight_trigger_1_up"),
        "trigger_1_down": df.Histo1D(("Haa_trigger_1Down", "Haa_trigger_1Down", 65, 70, 200), f"H_mass_{year}_{scope}", "weight_trigger_1_down"),
        "trigger_2_up": df.Histo1D(("Haa_trigger_2Up", "Haa_trigger_2Up", 65, 70, 200), f"H_mass_{year}_{scope}", "weight_trigger_2_up"),
        "trigger_2_down": df.Histo1D(("Haa_trigger_2Down", "Haa_trigger_2Down", 65, 70, 200), f"H_mass_{year}_{scope}", "weight_trigger_2_down"),
    }

    if scope == "mm":
        hists.update({
            "iso_1_up": df.Histo1D(("Haa_iso_1Up", "Haa_iso_1Up", 65, 70, 200), f"H_mass_{year}_{scope}", "weight_iso_1_up"),
            "iso_1_down": df.Histo1D(("Haa_iso_1Down", "Haa_iso_1Down", 65, 70, 200), f"H_mass_{year}_{scope}", "weight_iso_1_down"),
            "iso_2_up": df.Histo1D(("Haa_iso_2Up", "Haa_iso_2Up", 65, 70, 200), f"H_mass_{year}_{scope}", "weight_iso_2_up"),
            "iso_2_down": df.Histo1D(("Haa_iso_2Down", "Haa_iso_2Down", 65, 70, 200), f"H_mass_{year}_{scope}", "weight_iso_2_down"),
        })

    for name, hist in hists.items():
        hist.Write()

    nominal_hist = hists["nominal"].GetValue()
    print(f"Nominal histogram integral: {nominal_hist.Integral()}")