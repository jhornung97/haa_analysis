import sys
import ROOT
from ROOT import RooFit
import numpy as np
import os

scope = sys.argv[1]
year = sys.argv[2]

path = None

if scope == "mm" and int(year[:4]) <= 2018:
    path = f"/ceph/jhornung/Data_2018/{year}/single_muon_data/{scope}/"
    exclude = "single_muon_data.root"
    files = [path + f for f in os.listdir(path) if f.endswith(".root") and f != exclude]

    ROOT.ROOT.EnableImplicitMT(len(files))

    df = ROOT.RDataFrame("ntuple", files)

elif scope == "mm" and int(year[:4]) > 2018:
    path = f"/ceph/jhornung/Data_2018/{year}/muon_data/{scope}/"
    exclude = "muon_data.root"
    files = [path + f for f in os.listdir(path) if f.endswith(".root") and f != exclude]

    ROOT.ROOT.EnableImplicitMT(len(files))

    df = ROOT.RDataFrame("ntuple", files)

elif scope == "ee" and year == "2018":
    path = f"/ceph/jhornung/Data_2018/{year}/egamma_data/{scope}/"
    exclude = "egamma_data.root"
    files = [path + f for f in os.listdir(path) if f.endswith(".root") and f != exclude]

    ROOT.ROOT.EnableImplicitMT(len(files))

    df = ROOT.RDataFrame("ntuple", files)

elif scope == "ee" and (year == "2017" or year == "2016preVFP" or year == "2016postVFP"):
    path = f"/ceph/jhornung/Data_2018/{year}/single_electron_data/{scope}/"
    exclude = "single_electron_data.root"
    files = [path + f for f in os.listdir(path) if f.endswith(".root") and f != exclude]

    ROOT.ROOT.EnableImplicitMT(len(files))

    df = ROOT.RDataFrame("ntuple", files)

df = df.Define(f"H_mass_{year}_{scope}", "H_mass")
df = df.Define("nominal", "1")
df = df.Define("ps_mass_mask", "ps_1_mass < 3 && ps_2_mass < 3 && abs(ps_1_mass-ps_2_mass) < 0.06")
df = df.Define("H_mass_mask", "(H_mass > 70 && H_mass < 110) || (H_mass > 140 && H_mass < 200)")
df = df.Define("base_mask", "pt_vis > 30 && (m_vis > 75 && m_vis < 105) && abs(H_eta) < 2.4 && d1_pt > 10 && d2_pt > 10 && ps_mass_mask && H_mass_mask")

ref_count = df.Filter("base_mask").Count().GetValue()

threshold = 500 / 3000
cut = 10
step = 0.5

csv_path = f"/work/jhornung/Haa/alphas_{year}_{scope}.csv"

with open(csv_path, "w") as csv_file:
    continue_cutting = True
    while continue_cutting:
        iso_df = df.Define("iso_mask", f"d1_iso < {cut} && d2_iso < {cut}")
        iso_df = iso_df.Define("mask", "base_mask && iso_mask")

        n_pass = iso_df.Filter("mask").Count().GetValue()
        efficiency = n_pass / ref_count

        if efficiency < threshold:
            continue_cutting = False
            continue

        iso_df.Snapshot("ntuple", f"bkg_{scope}.root")

        f = ROOT.TFile.Open(f"bkg_{scope}.root")
        t = f.Get("ntuple")

        mass = ROOT.RooRealVar(f"H_mass_{year}_{scope}", "H_mass", 125, 70, 200)
        mask = ROOT.RooRealVar("mask", "mask", 0, 1)
        weight = ROOT.RooRealVar("nominal", "weight", -10, 10)

        data = ROOT.RooDataSet(f"data_H_mass_{year}_{scope}", "data", ROOT.RooArgSet(mass, mask), RooFit.Import(t), RooFit.Cut("mask==1"))

        print(f"iso cut {cut}: {data.sumEntries()} events")

        nbins = 65
        binning = ROOT.RooFit.Binning(nbins, 70, 200)
        mass.setRange("full", 70, 200)
        mass.setRange("loSB", 70, 110)
        mass.setRange("hiSB", 140, 200)
        fitRange = "loSB,hiSB"

        decorrelated_alpha = ROOT.RooRealVar(f"alpha_{year}_{scope}", "coeff of the exp func", -.02, -.025, -.015)
        decorrelated_model_bkg = ROOT.RooExponential(f"model_bkg_mass_H_mass_{year}_{scope}", "exp func", mass, decorrelated_alpha)
        decorrelated_result = decorrelated_model_bkg.fitTo(data, RooFit.Range(fitRange), RooFit.Save())

        alpha_val = decorrelated_result.floatParsFinal().find(f"alpha_{year}_{scope}").getVal()
        alpha_err = decorrelated_result.floatParsFinal().find(f"alpha_{year}_{scope}").getError()

        print(f"Fitted alpha: {alpha_val} +/- {alpha_err}")

        csv_file.write(f"{cut} {alpha_val} {alpha_err}\n")
        csv_file.flush()

        f.Close()

        cut -= step
