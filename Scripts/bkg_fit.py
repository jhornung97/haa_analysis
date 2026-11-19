import sys
import ROOT
from ROOT import RooFit
import numpy as np
import utility
import os
import json

ROOT.ROOT.EnableImplicitMT()

scope = sys.argv[1]
year = sys.argv[2]
#iso_mask = sys.argv[3] 

path = None 

if scope == "mm" and int(year) <= 2018:
    path = f"/ceph/jhornung/Data_2018/{year}/single_muon_data/{scope}/"
    exclude = "single_muon_data.root"
    files = [path + f for f in os.listdir(path) if f.endswith(".root") and f != exclude]

    df = ROOT.RDataFrame("ntuple", files)

if scope == "mm" and int(year) > 2018:
    path = f"/ceph/jhornung/Data_2018/{year}/muon_data/{scope}/"
    exclude = "muon_data.root"
    files = [path + f for f in os.listdir(path) if f.endswith(".root") and f != exclude]

    df = ROOT.RDataFrame("ntuple", files)

elif scope == "ee" and year == "2018":
    path = f"/ceph/jhornung/Data_2018/{year}/egamma_data/{scope}/"
    exclude = "egamma_data.root"
    files = [path + f for f in os.listdir(path) if f.endswith(".root") and f != exclude]
    
    df = ROOT.RDataFrame("ntuple", files)

elif scope == "ee" and (year == "2017" or year == "2016preVFP" or year == "2016postVFP"):
    path = f"/ceph/jhornung/Data_2018/{year}/single_electron_data/{scope}/"
    exclude = "single_electron_data.root"
    files = [path + f for f in os.listdir(path) if f.endswith(".root") and f != exclude]

    df = ROOT.RDataFrame("ntuple", files)

elif scope == "ee" and year == "2024":
    path = f"/ceph/jhornung/Data_2018/{year}/egamma_data/{scope}/"
    exclude = "egamma_data.root"
    files = [path + f for f in os.listdir(path) if f.endswith(".root") and f != exclude]

    df = ROOT.RDataFrame("ntuple", files)

df = df.Define(f"H_mass_{year}_{scope}", "H_mass")
df = df.Define("nominal", "1")
df = df.Define("ps_mass_mask", "ps_1_mass < 3 && ps_2_mass < 3 && abs(ps_1_mass-ps_2_mass) < 0.06")
#df = df.Define("iso_mask", f"d1_iso < {iso_mask} && d2_iso < {iso_mask}") # && d3_iso < {iso_mask} && d4_iso < {iso_mask}")
trigger_mask = "trg_single_mu" if scope == "mm" else "trg_single_el"
df = df.Define("mask", "pt_vis > 30 && (m_vis > 75 && m_vis < 105) && abs(H_eta) < 2.4 && d1_pt > 10 && d2_pt > 10 && ps_mass_mask && " + trigger_mask) #&& iso_mask")
snap = df.Snapshot("ntuple", f"/work/jhornung/Haa/limits/{year}/{scope}/bkg_{scope}.root")

f = ROOT.TFile.Open(f"/work/jhornung/Haa/limits/{year}/{scope}/bkg_{scope}.root")
t = f.Get("ntuple")

mass = ROOT.RooRealVar(f"H_mass_{year}_{scope}", "H_mass", 125, 70, 200)
mask = ROOT.RooRealVar("mask", "mask", 0, 1)
weight = ROOT.RooRealVar("nominal", "weight", -10, 10) 

data = ROOT.RooDataSet(f"data_H_mass_{year}_{scope}", "data", ROOT.RooArgSet(mass, mask), RooFit.Import(t),  RooFit.Cut("mask==1"))#, RooFit.WeightVar(weight))

print(data.sumEntries())

nbins = 65
binning = ROOT.RooFit.Binning(nbins, 70, 200)
mass.setRange("full", 70, 200)
mass.setRange("loSB", 70, 110)
mass.setRange("hiSB", 140, 200)
fitRange = "loSB,hiSB"

alpha = ROOT.RooRealVar("alpha", "coeff of the exp func", -.05, -2, 0)
decorrelated_alpha = ROOT.RooRealVar(f"alpha_{year}_{scope}", "coeff of the exp func", -.019, -0.03, -0.01)

model_bkg = ROOT.RooExponential(f"model_bkg_mass_H_mass_{year}_{scope}", "exp func", mass, alpha)
decorrelated_model_bkg = ROOT.RooExponential(f"model_bkg_mass_H_mass_{year}_{scope}", "exp func", mass, decorrelated_alpha)

# Load normalization values from JSON file
json_path = "/work/jhornung/Haa/iso_extrapolation/bkg_normalization.json"
with open(json_path, 'r') as f:
    norm_data = json.load(f)

# Get the normalization value for the current year and scope
if year in norm_data and scope in norm_data[year]:
    norm_value = norm_data[year][scope]["normalization"]
    # stat (+) extrapolating-function systematic; older JSON entries only have the stat error
    norm_error = norm_data[year][scope].get("normalization_total_error", norm_data[year][scope]["normalization_error"])
    norm_from_json = True
    print(f"Using sideband normalization from JSON: {norm_value} +/- {norm_error}")
else:
    # Fallback to data.sumEntries() if not found in JSON
    norm_value = data.sumEntries()
    norm_error = 0
    norm_from_json = False
    print(f"Warning: No normalization found in JSON for {year}/{scope}, using data.sumEntries(): {norm_value}")

result = model_bkg.fitTo(data, RooFit.Range(fitRange), RooFit.Save())
decorrelated_result = decorrelated_model_bkg.fitTo(data, RooFit.Range(fitRange), RooFit.Save())

# the iso extrapolation counts events in the sidebands only, while _norm is the yield over the full
# range the pdf is defined on -> scale up by the fitted shape: N_full = N_SB / (fraction of pdf in SB)
def sideband_fraction(pdf):
    norm_set = ROOT.RooArgSet(mass)
    return sum(pdf.createIntegral(norm_set, RooFit.NormSet(norm_set), RooFit.Range(r)).getVal() for r in fitRange.split(","))

def make_norm(pdf):
    scale = 1 / sideband_fraction(pdf) if norm_from_json else 1
    value, error = norm_value * scale, norm_error * scale
    print(f"{pdf.GetName()}: sideband -> full range scale {scale:.3f}, norm = {value:.1f} +/- {error:.1f}")
    norm_var = ROOT.RooRealVar(f"{pdf.GetName()}_norm", "number of background events", value, value - 3*error, value + 3*error)
    norm_var.setError(error)
    return norm_var

norm = make_norm(model_bkg)
decorrelated_norm = make_norm(decorrelated_model_bkg)

f_out = ROOT.TFile(f"/work/jhornung/Haa/limits/{year}/{scope}/workspace_bkg_{scope}.root", "RECREATE")
w_bkg = ROOT.RooWorkspace(f"workspace_bkg_{scope}",f"workspace_bkg_{scope}]")
getattr(w_bkg, "import")(data)
getattr(w_bkg, "import")(model_bkg)
getattr(w_bkg, "import")(norm)
w_bkg.Print()
w_bkg.Write()

decorrelated_w_bkg = ROOT.RooWorkspace(f"decorrelated_workspace_bkg_{scope}",f"decorrelated_workspace_bkg_{scope}]")
getattr(decorrelated_w_bkg, "import")(data)
getattr(decorrelated_w_bkg, "import")(decorrelated_model_bkg)
getattr(decorrelated_w_bkg, "import")(decorrelated_norm)
decorrelated_w_bkg.Print()
decorrelated_w_bkg.Write()

f_out.Close()

# full-range normalisations for the datacards (create_mass_datacards.py adds them as Gaussian constraints)
norm_summary = {
    f"{ws.GetName()}:{n.GetName()}": {"value": n.getVal(), "error": n.getError()}
    for ws, n in [(w_bkg, norm), (decorrelated_w_bkg, decorrelated_norm)]
}
with open(f"/work/jhornung/Haa/limits/{year}/{scope}/bkg_norm_{scope}.json", "w") as f:
    json.dump(norm_summary, f, indent=2)

# Apply global ROOT style once — must be done before any canvas is created
ROOT.gROOT.SetStyle("Plain")
ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetTextFont(42)        # Helvetica regular (Plain sets this to bold 62)
ROOT.gStyle.SetLabelFont(42, "xyz")
ROOT.gStyle.SetTitleFont(42, "xyz")
ROOT.gStyle.SetLegendFont(42)

legend = ROOT.TLegend(0.60, 0.44, 0.93, 0.68)

xframe = mass.frame()
xframe.SetTitle("")

data.plotOn(xframe, binning, ROOT.RooFit.MarkerColor(0), ROOT.RooFit.LineColor(0) )

model_bkg.plotOn(xframe, ROOT.RooFit.NormRange(fitRange), ROOT.RooFit.Range("full"), ROOT.RooFit.LineColor(ROOT.kRed))
sig_fit = xframe.getObject(int(xframe.numItems()-1))

data.plotOn(xframe, ROOT.RooFit.CutRange(fitRange), binning )
hist = xframe.getObject(int(xframe.numItems()-1))

legend.AddEntry(hist, "Data", "ep")
legend.AddEntry(sig_fit, "Bkg. Fit", "l")

legend.SetBorderSize(0)
legend.SetFillStyle(0)
legend.SetTextSize(0.042)

#xframe.SetMaximum(250)

model_bkg.Print("t")

num_params = result.floatParsFinal().getSize()
chi2 = xframe.chiSquare(num_params)

hpull = xframe.pullHist()
xframe2 = mass.frame()
xframe2.SetTitle("")
xframe2.addPlotable(hpull, "P")

xframe2.SetMinimum(-3)
xframe2.SetMaximum(3)

c = ROOT.TCanvas("canvas", "canvas", 1200, 1200)
c.Divide(1,2,0.1,0)

c.cd(1)
c.cd(1).SetTopMargin(0.08)
c.cd(1).SetRightMargin(0.05)
c.cd(1).SetBottomMargin(0.001)
c.cd(1).SetLeftMargin(0.12)
ROOT.gPad.SetBottomMargin(0.001)
c.cd(1).SetPad(0,0.3,1,1)

xframe.GetXaxis().SetLabelSize(0.05)
xframe.GetYaxis().SetTitleSize(0.06)
xframe.GetYaxis().SetLabelSize(0.06)
xframe.GetYaxis().SetTitleOffset(1.1)
xframe.GetYaxis().SetTitle("Events / 2 GeV")
xframe.SetMinimum(0)

xframe.Draw()
legend.Draw()

# CMS-style header labels
pw_label = ROOT.TLatex()
pw_label.SetNDC()
pw_label.SetTextFont(42)
pw_label.SetTextSize(0.055)
pw_label.DrawLatex(0.12, 0.93, "#it{Private Work} #bf{(CMS data)}")

lumi_map = {
    "2018": "59.8 fb^{-1} (13 TeV)",
    "2017": "41.5 fb^{-1} (13 TeV)",
    "2016preVFP": "19.5 fb^{-1} (13 TeV)",
    "2016postVFP": "16.8 fb^{-1} (13 TeV)",
    "2024": "110.0 fb^{-1} (13.6 TeV)",
}
lumi_label = ROOT.TLatex()
lumi_label.SetNDC()
lumi_label.SetTextFont(42)
lumi_label.SetTextSize(0.055)
lumi_label.SetTextAlign(31)
lumi_label.DrawLatex(0.95, 0.93, lumi_map.get(year, ""))

scope_str = "ee SR" if scope == "ee" else "#mu#mu SR"
ann_label = ROOT.TLatex()
ann_label.SetNDC()
ann_label.SetTextFont(42)
ann_label.SetTextSize(0.048)
for y_pos, text in [
    (0.86, "H #rightarrow aa #rightarrow KKKK"),
    (0.80, "Bkg Est."),
    (0.74, scope_str),
    (0.68, "#chi^{2}/dof = " + f"{chi2:.2f}"),
]:
    ann_label.DrawLatex(0.60, y_pos, text)

c.cd(2)
c.cd(2).SetTopMargin(0.001)
c.cd(2).SetRightMargin(0.05)
c.cd(2).SetBottomMargin(0.3)
c.cd(2).SetLeftMargin(0.12)
c.cd(2).SetPad(0,0,1,0.3)

xframe2.GetXaxis().SetTitle("m_{4K} in GeV")
xframe2.GetYaxis().SetTitle("#frac{Data - Fit}{Error}")
xframe2.GetXaxis().SetTitleSize(0.14)
xframe2.GetXaxis().SetTitleOffset(1)
xframe2.GetXaxis().SetLabelSize(0.11)
xframe2.GetYaxis().SetTitleSize(0.11)
xframe2.GetYaxis().SetTitleOffset(0.35)
xframe2.GetYaxis().SetLabelSize(0.11)

xframe2.Draw()

zero = ROOT.TLine(70, 0, 200, 0)
zero.SetLineColor(ROOT.kRed)
zero.SetLineWidth(2)
zero.SetLineStyle(2)
zero.Draw()

c.SaveAs(f"/web/jhornung/public_html/analysis_plots/{year}/bkg_fit_{year}_{scope}.png")
c.SaveAs(f"/web/jhornung/public_html/analysis_plots/{year}/bkg_fit_{year}_{scope}.pdf")

print(data.sumEntries())
print(data.numEntries())