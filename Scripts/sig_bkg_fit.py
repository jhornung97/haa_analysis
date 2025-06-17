import sys
import ROOT
from ROOT import RooFit
import numpy as np
import utility
'''
bkg_df = ROOT.RDataFrame("ntuple", sys.argv[1])
bkg_df = bkg_df.Define("signs", "genWeight > 0 ? 1 : -1")
bkg_df = bkg_df.Define("Weight", "59e3*signs*evtweight")
bkg_df = bkg_df.Define("mask", "pt_vis > 30 && (m_vis > 75 && m_vis < 105)")
bkg_snapshot = bkg_df.Snapshot("ntuple", "bkg_snapshot.root")
'''
bkg_file = ROOT.TFile.Open("bkg_snapshot.root")
bkg_tree = bkg_file.Get("ntuple")

sig_root_file = ROOT.TFile("/ceph/jhornung/analysis/2018/HaaKKKK_M125_13TeV_powheg_pythia8-RunIISummer20UL18MiniAODv2-106X/mm/HaaKKKK_M125_13TeV_powheg_pythia8-RunIISummer20UL18MiniAODv2-106X.root")
sig_tree = sig_root_file.Get("ntuple")

sig_H_mass = sig_tree.AsMatrix(columns=["H_mass"]).flatten()
sig_weight = np.sign(sig_tree.AsMatrix(columns=["genWeight"]).flatten())

x = ROOT.RooRealVar("H_mass", "H_mass", 110, 140)
w = ROOT.RooRealVar("Weight", "Weight", -26000, 26000)
mask = ROOT.RooRealVar("mask", "mask", 0, 1)

data = ROOT.RooDataSet("data", "data", ROOT.RooArgSet(x, w, mask), RooFit.Import(bkg_tree), RooFit.WeightVar(w), RooFit.Cut("mask==1"))

for i, weight in enumerate(sig_weight[(sig_H_mass > 110) & (sig_H_mass < 140)][:5000]):
    w.setVal(weight)
    x.setVal(sig_H_mass[(sig_H_mass > 110) & (sig_H_mass < 140)][i])
    data.add(ROOT.RooArgSet(x))

#define bkg model
alpha = ROOT.RooRealVar("alpha", "coeff of the exp func", -0.1, -10.0, 0.0)
truebkg = ROOT.RooExponential("truebkg", "exp func", x, alpha)

#define sig model
#true signal
MH = ROOT.RooRealVar("mHiggs", "mHiggs", 125.0, 120.0, 130.0)
dMH = ROOT.RooRealVar("dMHiggs", "dMHiggs", 0.0, -1.0, 1.0)
mean = ROOT.RooFormulaVar("mean", "mean", "@0+@1", ROOT.RooArgList(MH, dMH))
sigma1 = ROOT.RooRealVar("sigma1", "width of gaussians", 1.0, 0.0, 2.0)
sigma2 = ROOT.RooRealVar("sigma2", "width of gaussians", 2.0, 0.0, 3.0)

sig1 = ROOT.RooGaussian("sig1", "Signal component 1", x, mean, sigma1)
sig2 = ROOT.RooGaussian("sig2", "Signal component 2", x, mean, sigma2)

sig1frac = ROOT.RooRealVar("sig1frac", "fraction of component 1 in signal", 0.8, 0.0, 1.0)
sig = ROOT.RooAddPdf("true signal", "true gaussian Signal", ROOT.RooArgList(sig1, sig2), ROOT.RooArgList(sig1frac))

#signal shoulder
a0 = ROOT.RooRealVar("a0", "a0", 0.5, -1, 1)
a1 = ROOT.RooRealVar("a1", "a1", -0.2, -2, 2.0)
a2 = ROOT.RooRealVar("a2", "a2", .5, -1, 1.0)
a3 = ROOT.RooRealVar("a3", "a3", -.5, -1, 1.0)
a4 = ROOT.RooRealVar("a4", "a4", -.5, -1, 1.0)
a5 = ROOT.RooRealVar("a5", "a5", -.5, -1, 1.0)
a6 = ROOT.RooRealVar("a6", "a6", -.5, -1, 1.0)
shoulder = ROOT.RooChebychev("shoulder", "shoulder from misrecoed signal events", x, ROOT.RooArgList(a0))#, a1, a2, a3, a4))#, a5, a6))

combi_mean = ROOT.RooRealVar("combi_mean", "combi_mean", 118.0, 115.0, 120.0)
combisigma = ROOT.RooRealVar("combisigma", "combisigma", 10.0, 0.0, 20.0)
combi_peak = ROOT.RooGaussian("combi_peak", "combi_peak", x, combi_mean, combisigma)

#full model
ntruesig = ROOT.RooRealVar("nsig", "number of signal events", 250, 0., data.numEntries())
nshoulder = ROOT.RooRealVar("nshoulder", "number of shoulder events", 1000, 0., data.numEntries())
ncombi = ROOT.RooRealVar("ncombi", "number of combi events", 1000, 0., data.numEntries())
nbkg = ROOT.RooRealVar("nbkg", "number of non signal events", 200000, 0., data.numEntries())
model = ROOT.RooAddPdf("model", "bkg+(gauss1+gauss2)", ROOT.RooArgList(truebkg, shoulder, sig, combi_peak), ROOT.RooArgList(nbkg, nshoulder, ntruesig, ncombi))

data.Print()

result = model.fitTo(data, RooFit.Save(), RooFit.SumW2Error(True), RooFit.PrintLevel(-1))
model.Print("t")

ROOT.gROOT.SetStyle("Plain")
ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

c = ROOT.TCanvas("canvas", "canvas", 1200, 1200)
c.Divide(1,2,0.1,0)

legend = ROOT.TLegend(0.6, 0.6, 0.9, 0.8)

xframe = x.frame()
xframe.SetTitle("")

hist = data.plotOn(xframe, RooFit.Binning(60))
legend.AddEntry(hist.getObject(0), "Data", "lep")

#model.plotOn(xframe, ROOT.RooFit.VisualizeError(result, 2), ROOT.RooFit.FillColor(ROOT.kYellow), ROOT.RooFit.FillStyle(3001), ROOT.RooFit.LineColor(ROOT.kYellow))
#model.plotOn(xframe, ROOT.RooFit.VisualizeError(result, 1), ROOT.RooFit.FillColor(ROOT.kGreen), ROOT.RooFit.FillStyle(3001), ROOT.RooFit.LineColor(ROOT.kGreen))

fit = model.plotOn(xframe, RooFit.LineColor(ROOT.kRed), RooFit.MoveToBack())
legend.AddEntry(fit.getObject(0), "S+B fit", "l")

bkg_fit = model.plotOn(xframe, ROOT.RooFit.Components("truebkg"), ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.MoveToBack())
legend.AddEntry(bkg_fit.getObject(0), "B component", "l")

legend.SetBorderSize(0)
legend.SetFillStyle(0)
legend.SetTextSize(0.04)

num_params = result.floatParsFinal().getSize()
chi2 = xframe.chiSquare(num_params)
tl = ROOT.TLatex()
tl.SetTextAlign(12)
tl.SetTextSize(0.04)
tl.DrawLatex(0.1,0.95, "CMS #bf{#it{Simulation Work in progress}}")
tl.DrawLatex(0.725,0.95, "#sqrt{s} = 13 GeV")
tl.DrawLatex(0.6,0.9, "H #rightarrow aa #rightarrow KKKK")
#tl.DrawLatex(0.6,0.85, "Sig MC + Bkg MC")
tl.DrawLatex(0.6,0.7, "#chi^{2}/dof = " + str('%.2f'%((chi2))))

c.cd(1)
c.cd(1).SetTopMargin(0.1)
c.cd(1).SetRightMargin(0.1)
c.cd(1).SetBottomMargin(0.001)
c.cd(1).SetLeftMargin(0.1)
ROOT.gPad.SetBottomMargin(0.001)
c.cd(1).SetPad(0,0.3,1,1)

xframe.GetXaxis().SetLabelSize(0.05)
xframe.GetXaxis().SetTitleSize(0.04)
xframe.GetYaxis().SetTitleSize(0.04)
xframe.GetYaxis().SetLabelSize(0.035)
xframe.GetYaxis().SetTitleOffset(0.9)

xframe.Draw()
legend.Draw()

c.cd(2)
c.cd(2).SetTopMargin(0.001)
c.cd(2).SetRightMargin(0.1)
c.cd(2).SetBottomMargin(0.3)
c.cd(2).SetLeftMargin(0.1)
c.cd(2).SetPad(0,0,1,0.3)

hpull = xframe.pullHist()
xframe2 = x.frame()
xframe2.SetTitle("")
xframe2.addPlotable(hpull, "P")

xframe2.GetXaxis().SetTitle("M_{KKKK} in GeV")
xframe2.GetYaxis().SetTitle("MC - Fit")
xframe2.GetXaxis().SetTitleSize(0.075)
xframe2.GetXaxis().SetTitleOffset(1)
xframe2.GetXaxis().SetLabelSize(0.08)
xframe2.GetYaxis().SetTitleSize(0.075)
xframe2.GetYaxis().SetTitleOffset(0.35)
xframe2.GetYaxis().SetLabelSize(0.08)

xframe2.Draw()

zero = ROOT.TLine(110, 0, 140, 0)
zero.SetLineColor(ROOT.kRed)
zero.SetLineWidth(2)
zero.SetLineStyle(2)
zero.Draw()

c.SaveAs(sys.argv[2] + ".png")
c.SaveAs(sys.argv[2] + ".pdf")
