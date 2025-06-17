import sys
import ROOT
from ROOT import RooFit
import numpy as np
import utility
import os

ROOT.ROOT.EnableImplicitMT()

scope = sys.argv[1]
year = sys.argv[2]
iso_mask = sys.argv[3] 

path = None 

if scope == "mm":
    path = f"/ceph/jhornung/Data_2018/{year}/single_muon_data/{scope}/"
    exclude = "single_muon_data.root"
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
    
df = df.Define(f"H_mass_{year}_{scope}", "H_mass")
df = df.Define("nominal", "1")
df = df.Define("ps_mass_mask", "ps_1_mass < 3 && ps_2_mass < 3 && abs(ps_1_mass-ps_2_mass) < 0.06")
df = df.Define("iso_mask", f"d1_iso < {iso_mask} && d2_iso < {iso_mask}") # && d3_iso < {iso_mask} && d4_iso < {iso_mask}")
df = df.Define("mask", "pt_vis > 30 && (m_vis > 75 && m_vis < 105) && abs(H_eta) < 2.4 && d1_pt > 10 && d2_pt > 10 && ps_mass_mask && iso_mask")
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
mass.setRange("hiSB", 130, 200)
fitRange = "loSB,hiSB"

alpha = ROOT.RooRealVar("alpha", "coeff of the exp func", -.05, -2, 0)
decorrelated_alpha = ROOT.RooRealVar(f"alpha_{year}_{scope}", "coeff of the exp func", -.019, -.0205, -.0175)

model_bkg = ROOT.RooExponential(f"model_bkg_mass_H_mass_{year}_{scope}", "exp func", mass, alpha)
decorrelated_model_bkg = ROOT.RooExponential(f"model_bkg_mass_H_mass_{year}_{scope}", "exp func", mass, decorrelated_alpha)

norm = ROOT.RooRealVar(f"model_bkg_mass_H_mass_{year}_{scope}_norm", "number of background events", data.sumEntries(), 0, 3*data.sumEntries())
decorrelated_norm = ROOT.RooRealVar(f"model_bkg_mass_H_mass_{year}_{scope}_norm", "number of background events", data.sumEntries(), 0, 3*data.sumEntries())

result = model_bkg.fitTo(data, RooFit.Range(fitRange), RooFit.Save())
decorrelated_result = decorrelated_model_bkg.fitTo(data, RooFit.Range(fitRange), RooFit.Save())

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

legend = ROOT.TLegend(0.6, 0.35, 0.875, 0.6)

xframe = mass.frame()
xframe.SetTitle("")

data.plotOn(xframe, binning, ROOT.RooFit.MarkerColor(0), ROOT.RooFit.LineColor(0) )

model_bkg.plotOn(xframe, ROOT.RooFit.NormRange(fitRange), ROOT.RooFit.Range("full"), ROOT.RooFit.LineColor(2))
sig_fit = xframe.getObject(int(xframe.numItems()-1))
legend.AddEntry(sig_fit, "Bkg. Fit", "l")

data.plotOn(xframe, ROOT.RooFit.CutRange(fitRange), binning )
hist = xframe.getObject(int(xframe.numItems()-1))
legend.AddEntry(hist, "Data", "ep")

legend.SetBorderSize(0)
legend.SetFillStyle(0)
legend.SetTextSize(0.04)


model_bkg.Print("t")

num_params = result.floatParsFinal().getSize()
chi2 = xframe.chiSquare(num_params)

hpull = xframe.pullHist()
xframe2 = mass.frame()
xframe2.SetTitle("")
xframe2.addPlotable(hpull, "P")

xframe2.SetMinimum(-3)
xframe2.SetMaximum(3)
 
ROOT.gROOT.SetStyle("Plain")
ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

c = ROOT.TCanvas("canvas", "canvas", 1200, 1200)

c.Divide(1,2,0.1,0)

if year == "2018":
    lumilabel = "59.83"
elif year == "2017":
    lumilabel = "41.53"

tl = ROOT.TLatex()
tl.SetTextAlign(12)
tl.SetTextSize(0.04)
tl.DrawLatex(0.1,0.95, "#bf{#it{Private Work}} (CMS data)")
if year == "2018":
    tl.DrawLatex(0.57,0.95, "#bf{58.8 fb^{-1}, 2018 (13 TeV)}")
elif year == "2017":
    tl.DrawLatex(0.55,0.95, "41.48 fb^{-1}, 2017 (13 TeV)")
tl.DrawLatex(0.6,0.9, "H #rightarrow aa #rightarrow KKKK")
tl.DrawLatex(0.6,0.85, "Bkg Est.")
if scope == "ee":
    tl.DrawLatex(0.6,0.8, "ee SR")
elif scope == "mm":
    tl.DrawLatex(0.6,0.8, "#mu#mu SR")
tl.DrawLatex(0.6,0.75, "#chi^{2}/dof = " + str('%.2f'%((chi2))))

c.cd(1)
c.cd(1).SetTopMargin(0.1)
c.cd(1).SetRightMargin(0.1)
c.cd(1).SetBottomMargin(0.001)
c.cd(1).SetLeftMargin(0.1)
ROOT.gPad.SetBottomMargin(0.001)
c.cd(1).SetPad(0,0.3,1,1)

xframe.GetXaxis().SetLabelSize(0.05)
#xframe.GetXaxis().SetTitleSize(0.04)
xframe.GetYaxis().SetTitleSize(0.06)
xframe.GetYaxis().SetLabelSize(0.06)
xframe.GetYaxis().SetTitleOffset(0.9)
xframe.GetYaxis().SetTitle("Events / 2 GeV")

xframe.Draw()
legend.Draw()

c.cd(2)
c.cd(2).SetTopMargin(0.001)
c.cd(2).SetRightMargin(0.1)
c.cd(2).SetBottomMargin(0.3)
c.cd(2).SetLeftMargin(0.1)
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

c.SaveAs(f"/web/jhornung/public_html/analysis_plots/bkg_fit_{year}_{scope}.png")
c.SaveAs(f"/web/jhornung/public_html/analysis_plots/bkg_fit_{year}_{scope}.pdf")

print(data.sumEntries())
print(data.numEntries())