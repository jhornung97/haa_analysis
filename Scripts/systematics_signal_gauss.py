import sys
import subprocess
import ROOT
from ROOT import RooFit
import numpy as np
import matplotlib.pyplot as plt
import utility
from array import array

def prep_data(scope, year):
    lumi = 0
    if year == "2018":
        lumi = 59.83e3
    elif year == "2017":
        lumi = 41.48e3
    elif year == "2016preVFP":
        lumi = 19.5e3
    elif year == "2016postVFP":
        lumi = 16.8e3
    if scope == "mm":
        df = ROOT.RDataFrame("ntuple", "/ceph/jhornung/MC_2018/2018/HaaKKKK_M125_13TeV_powheg_pythia8-RunIISummer20UL18MiniAODv2-106X/mm/HaaKKKK_M125_13TeV_powheg_pythia8-RunIISummer20UL18MiniAODv2-106X.root")
        df = df.Define(f"H_mass_{year}_{scope}", "H_mass")
        df = df.Define("signs", "genWeight > 0 ? 1 : -1")
        df = df.Define("nominal", f'{lumi}*signs*evtweight*id_wgt_mu_1*id_wgt_mu_2*iso_wgt_mu_1*iso_wgt_mu_2*trigger_wgt_mu_1*trigger_wgt_mu_2*puweight')#*evtweight
        df = df.Define("weight_pu_up", f'{lumi}*signs*evtweight*id_wgt_mu_1*id_wgt_mu_2*iso_wgt_mu_1*iso_wgt_mu_2*trigger_wgt_mu_1*trigger_wgt_mu_2*puweight_up')#*evtweight
        df = df.Define("weight_pu_down", f'{lumi}*signs*evtweight*id_wgt_mu_1*id_wgt_mu_2*iso_wgt_mu_1*iso_wgt_mu_2*trigger_wgt_mu_1*trigger_wgt_mu_2*puweight_down')#*evtweight
        df = df.Define("weight_id_1_up", f'{lumi}*signs*evtweight*id_wgt_mu_1__MuonIDUp*id_wgt_mu_2*iso_wgt_mu_1*iso_wgt_mu_2*trigger_wgt_mu_1*trigger_wgt_mu_2*puweight')#*evtweight
        df = df.Define("weight_id_1_down", f'{lumi}*signs*evtweight*id_wgt_mu_1__MuonIDDown*id_wgt_mu_2*iso_wgt_mu_1*iso_wgt_mu_2*trigger_wgt_mu_1*trigger_wgt_mu_2*puweight')#*evtweight
        df = df.Define("weight_id_2_up", f'{lumi}*signs*evtweight*id_wgt_mu_1*id_wgt_mu_2__MuonIDUp*iso_wgt_mu_1*iso_wgt_mu_2*trigger_wgt_mu_1*trigger_wgt_mu_2*puweight')#*evtweight
        df = df.Define("weight_id_2_down", f'{lumi}*signs*evtweight*id_wgt_mu_1*id_wgt_mu_2__MuonIDDown*iso_wgt_mu_1*iso_wgt_mu_2*trigger_wgt_mu_1*trigger_wgt_mu_2*puweight')#*evtweight
        df = df.Define("weight_iso_1_up", f'{lumi}*signs*evtweight*id_wgt_mu_1*id_wgt_mu_2*iso_wgt_mu_1__MuonIsoUp*iso_wgt_mu_2*trigger_wgt_mu_1*trigger_wgt_mu_2*puweight')#*evtweight
        df = df.Define("weight_iso_1_down", f'{lumi}*signs*evtweight*id_wgt_mu_1*id_wgt_mu_2*iso_wgt_mu_1__MuonIsoDown*iso_wgt_mu_2*trigger_wgt_mu_1*trigger_wgt_mu_2*puweight')#*evtweight
        df = df.Define("weight_iso_2_up", f'{lumi}*signs*evtweight*id_wgt_mu_1*id_wgt_mu_2*iso_wgt_mu_1*iso_wgt_mu_2__MuonIsoUp*trigger_wgt_mu_1*trigger_wgt_mu_2*puweight')#*evtweight
        df = df.Define("weight_iso_2_down", f'{lumi}*signs*evtweight*id_wgt_mu_1*id_wgt_mu_2*iso_wgt_mu_1*iso_wgt_mu_2__MuonIsoDown*trigger_wgt_mu_1*trigger_wgt_mu_2*puweight')#*evtweight
        df = df.Define("weight_trigger_1_up", f'{lumi}*signs*evtweight*id_wgt_mu_1*id_wgt_mu_2*iso_wgt_mu_1*iso_wgt_mu_2*trigger_wgt_mu_1__MuonTriggerUp*trigger_wgt_mu_2*puweight')#*evtweight
        df = df.Define("weight_trigger_1_down", f'{lumi}*signs*evtweight*id_wgt_mu_1*id_wgt_mu_2*iso_wgt_mu_1*iso_wgt_mu_2*trigger_wgt_mu_1__MuonTriggerDown*trigger_wgt_mu_2*puweight')#*evtweight
        df = df.Define("weight_trigger_2_up", f'{lumi}*signs*evtweight*id_wgt_mu_1*id_wgt_mu_2*iso_wgt_mu_1*iso_wgt_mu_2*trigger_wgt_mu_1*trigger_wgt_mu_2__MuonTriggerUp*puweight')#*evtweight
        df = df.Define("weight_trigger_2_down", f'{lumi}*signs*evtweight*id_wgt_mu_1*id_wgt_mu_2*iso_wgt_mu_1*iso_wgt_mu_2*trigger_wgt_mu_1*trigger_wgt_mu_2__MuonTriggerDown*puweight')#*evtweight
        df = df.Define("ps_mass_mask", "ps_1_mass < 3 && ps_2_mass < 3 && abs(ps_1_mass-ps_2_mass) < 0.06")
        df = df.Define("mask", "pt_vis > 30 && (m_vis > 75 && m_vis < 105) && abs(H_eta) < 2.4 && d1_pt > 10 && d2_pt > 10 && ps_mass_mask")
        snap = df.Snapshot("ntuple", f"/work/jhornung/Haa/limits/{year}/{scope}/signal_mm.root")
    elif scope == "ee":
        df = ROOT.RDataFrame("ntuple", "/ceph/jhornung/MC_2018/2018/HaaKKKK_M125_13TeV_powheg_pythia8-RunIISummer20UL18MiniAODv2-106X/ee/HaaKKKK_M125_13TeV_powheg_pythia8-RunIISummer20UL18MiniAODv2-106X.root")
        df = df.Define(f"H_mass_{year}_{scope}", "H_mass")
        df = df.Define("signs", "genWeight > 0 ? 1 : -1")
        df = df.Define("sfs", 'id_wgt_ele_1*id_wgt_ele_2*trigger_wgt_ele_1*trigger_wgt_ele_2')
        df = df.Define("nominal", f'{lumi}*signs*evtweight*sfs*puweight')
        df = df.Define("weight_pu_up", f'{lumi}*signs*evtweight*sfs*puweight_up')
        df = df.Define("weight_pu_down", f'{lumi}*signs*evtweight*sfs*puweight_down')
        df = df.Define("weight_id_1_up", f'{lumi}*signs*evtweight*id_wgt_ele_1__ElectronIDUp*id_wgt_ele_2*trigger_wgt_ele_1*trigger_wgt_ele_2*puweight')
        df = df.Define("weight_id_1_down", f'{lumi}*signs*evtweight*id_wgt_ele_1__ElectronIDDown*id_wgt_ele_2*trigger_wgt_ele_1*trigger_wgt_ele_2*puweight')
        df = df.Define("weight_id_2_up", f'{lumi}*signs*evtweight*id_wgt_ele_1*id_wgt_ele_2__ElectronIDUp*trigger_wgt_ele_1*trigger_wgt_ele_2*puweight')
        df = df.Define("weight_id_2_down", f'{lumi}*signs*evtweight*id_wgt_ele_1*id_wgt_ele_2__ElectronIDDown*trigger_wgt_ele_1*trigger_wgt_ele_2*puweight')
        df = df.Define("weight_iso_1_up", f'{lumi}*signs*evtweight*id_wgt_ele_1*id_wgt_ele_2*trigger_wgt_ele_1*trigger_wgt_ele_2*puweight')
        df = df.Define("weight_iso_1_down", f'{lumi}*signs*evtweight*id_wgt_ele_1*id_wgt_ele_2*trigger_wgt_ele_1*trigger_wgt_ele_2*puweight')
        df = df.Define("weight_trigger_1_up", f'{lumi}*signs*evtweight*id_wgt_ele_1*id_wgt_ele_2*trigger_wgt_ele_1_up*trigger_wgt_ele_2*puweight')
        df = df.Define("weight_trigger_1_down", f'{lumi}*signs*evtweight*id_wgt_ele_1*id_wgt_ele_2*trigger_wgt_ele_1_down*trigger_wgt_ele_2*puweight')
        df = df.Define("weight_trigger_2_up", f'{lumi}*signs*evtweight*id_wgt_ele_1*id_wgt_ele_2*trigger_wgt_ele_1*trigger_wgt_ele_2_up*puweight')
        df = df.Define("weight_trigger_2_down", f'{lumi}*signs*evtweight*id_wgt_ele_1*id_wgt_ele_2*trigger_wgt_ele_1*trigger_wgt_ele_2_down*puweight')
        df = df.Define("ps_mass_mask", "ps_1_mass < 3 && ps_2_mass < 3 && abs(ps_1_mass-ps_2_mass) < 1")
        df = df.Define("mask", "pt_vis > 30 && (m_vis > 75 && m_vis < 105) && abs(H_eta) < 2.4 && d1_pt > 10 && d2_pt > 10 && ps_mass_mask")
        snap = df.Snapshot("ntuple", f"/work/jhornung/Haa/limits/{year}/{scope}/signal_ee.root")

scope = sys.argv[1]
year = sys.argv[2]

prep_data(scope, year)

mass = ROOT.RooRealVar(f"H_mass_{year}_{scope}", "H_mass", 125, 115, 135)
'''
bin_edges_1 = np.linspace(70, 120, 25, endpoint=False)
bin_edges_2 = np.linspace(120, 130, 20, endpoint=False)
bin_edges_3 = np.linspace(130, 200, 36)

bin_edges = np.concatenate((bin_edges_1, bin_edges_2, bin_edges_3))

custom_binning = ROOT.RooBinning(len(bin_edges)-1, array('d', bin_edges))
'''
mass.setBins(40)

mass.getBinning().Print()

mask = ROOT.RooRealVar("mask", "mask", 0, 1)
weight_nominal = ROOT.RooRealVar("nominal", "weight", -10, 10) 
weight_pu_up = ROOT.RooRealVar("weight_pu_up", "weight_pu_up", -10, 10)
weight_pu_down = ROOT.RooRealVar("weight_pu_down", "weight_pu_down", -10, 10)
weight_id_1_up = ROOT.RooRealVar("weight_id_1_up", "weight_id_1_up", -10, 10)
weight_id_1_down = ROOT.RooRealVar("weight_id_1_down", "weight_id_1_down", -10, 10)
weight_id_2_up = ROOT.RooRealVar("weight_id_2_up", "weight_id_2_up", -10, 10)
weight_id_2_down = ROOT.RooRealVar("weight_id_2_down", "weight_id_2_down", -10, 10)
weight_iso_1_up = ROOT.RooRealVar("weight_iso_1_up", "weight_iso_1_up", -10, 10)
weight_iso_1_down = ROOT.RooRealVar("weight_iso_1_down", "weight_iso_1_down", -10, 10)
weight_iso_2_up = ROOT.RooRealVar("weight_iso_2_up", "weight_iso_2_up", -10, 10)
weight_iso_2_down = ROOT.RooRealVar("weight_iso_2_down", "weight_iso_2_down", -10, 10)
weight_trigger_1_up = ROOT.RooRealVar("weight_trigger_1_up", "weight_trigger_1_up", -10, 10)
weight_trigger_1_down = ROOT.RooRealVar("weight_trigger_1_down", "weight_trigger_1_down", -10, 10)
weight_trigger_2_up = ROOT.RooRealVar("weight_trigger_2_up", "weight_trigger_2_up", -10, 10)
weight_trigger_2_down = ROOT.RooRealVar("weight_trigger_2_down", "weight_trigger_2_down", -10, 10)

if scope == "mm":

    f = ROOT.TFile.Open(f"/work/jhornung/Haa/limits/{year}/{scope}/signal_{scope}.root")
    t = f.Get("ntuple")
    mc = {}
    yield_vars = {}
    yields = {}

    mc["nominal"] = ROOT.RooDataSet('data_nominal', 'data_nominal', ROOT.RooArgSet(mass, weight_nominal, mask), RooFit.Import(t), RooFit.WeightVar(weight_nominal), RooFit.Cut("mask==1"))
    yields["nominal"] = mc["nominal"].sumEntries()

    systs = [weight_pu_up, weight_pu_down, weight_id_1_up, weight_id_1_down, weight_id_2_up, weight_id_2_down, weight_iso_1_up, weight_iso_1_down, weight_iso_2_up, weight_iso_2_down, weight_trigger_1_up, weight_trigger_1_down, weight_trigger_2_up, weight_trigger_2_down]

    for syst in systs:
        mc[syst.GetName()] = ROOT.RooDataSet(f'data_{syst.GetName()}', f'data_{syst.GetName()}', ROOT.RooArgSet(mass, syst, mask), RooFit.Import(t), RooFit.WeightVar(syst), RooFit.Cut("mask==1"))
        yields[syst.GetName()] = mc[syst.GetName()].sumEntries()

        yield_var = mc[syst.GetName()].sumEntries()/mc["nominal"].sumEntries()
        yield_vars[syst.GetName()] = yield_var

if scope == "ee":
    f = ROOT.TFile.Open(f"/work/jhornung/Haa/limits/{year}/{scope}/signal_{scope}.root")
    t = f.Get("ntuple")
    mc = {}
    yield_vars = {}
    yields = {}
    
    mc["nominal"] = ROOT.RooDataSet('data_nominal', 'data_nominal', ROOT.RooArgSet(mass, weight_nominal, mask), RooFit.Import(t), RooFit.WeightVar(weight_nominal), RooFit.Cut("mask==1"))
    yields["nominal"] = mc["nominal"].sumEntries()
    systs = [weight_pu_up, weight_pu_down, weight_id_1_up, weight_id_1_down, weight_id_2_up, weight_id_2_down, weight_trigger_1_up, weight_trigger_1_down, weight_trigger_2_up, weight_trigger_2_down]

    for syst in systs:
        mc[syst.GetName()] = ROOT.RooDataSet(f'data_{syst.GetName()}', f'data_{syst.GetName()}', ROOT.RooArgSet(mass, syst, mask), RooFit.Import(t), RooFit.WeightVar(syst), RooFit.Cut("mask==1"))
        yields[syst.GetName()] = mc[syst.GetName()].sumEntries()

        yield_var = mc[syst.GetName()].sumEntries()/mc["nominal"].sumEntries()
        yield_vars[syst.GetName()] = yield_var


MH = ROOT.RooRealVar("MH", "MH", 125.0, 120.0, 130.0)
sigma = ROOT.RooRealVar("sigma", "sigma", 1.5, 0.1, 10.0)
xi = ROOT.RooRealVar("xi", "xi", -0.2, -2.0, 0.0)
rhoL = ROOT.RooRealVar("rhoL", "rhoL", 0.1, 0.0, 2.0)
rhoR = ROOT.RooRealVar("rhoR", "rhoR", 0.0, -0.5, 0.5)
 
sig = ROOT.RooBukinPdf("sig", "sig", mass, MH, sigma, xi, rhoL, rhoR)

result = sig.fitTo(mc["nominal"], RooFit.Save(), RooFit.PrintLevel(-1), RooFit.Strategy(2), RooFit.SumW2Error(True), RooFit.PrintEvalErrors(-1))


c = ROOT.TCanvas("canvas", "canvas", 1200, 1200)
c.Divide(1,2,0.1,0)

legend = ROOT.TLegend(0.625, 0.35, 0.875, 0.6)

xframe = mass.frame()

mc["nominal"].plotOn(xframe)

#model_nom.plotOn(xframe, ROOT.RooFit.Components("offset"), ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.LineStyle(ROOT.kDashed))
#bkg_fit = xframe.getObject(int(xframe.numItems()-1))
#legend.AddEntry(bkg_fit, "Combinatorics Bkg", "l")

#model.plotOn(xframe, ROOT.RooFit.Components("combinatorics"), ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.LineStyle(ROOT.kDotted))
#comb_peak = xframe.getObject(int(xframe.numItems()-1))
#legend.AddEntry(comb_peak, "Combinatorics Peak", "l")

#model.plotOn(xframe, ROOT.RooFit.Components("bkg+combinatorics"), ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.LineStyle(ROOT.kDashDotted))
#comb_cont = xframe.getObject(int(xframe.numItems()-1))
#legend.AddEntry(comb_cont, "Combinatorics Cont", "l")

sig.plotOn(xframe, RooFit.LineColor(ROOT.kRed), RooFit.LineStyle(ROOT.kSolid))
sig_fit = xframe.getObject(int(xframe.numItems()-1))
legend.AddEntry(sig_fit, "Signal Fit", "l")

mc["nominal"].plotOn(xframe)
hist = xframe.getObject(int(xframe.numItems()-1))
legend.AddEntry(hist, "Simulated Data", "ep")

legend.SetBorderSize(0)
legend.SetFillStyle(0)
legend.SetTextSize(0.04)

num_params = int(result.constPars().getSize()/2)
chi2 = xframe.chiSquare(num_params)

hpull = xframe.pullHist()
xframe2 = mass.frame()
xframe2.SetTitle("")
xframe2.addPlotable(hpull, "P")

#hpull_error1 = xframe.pullHist("model_Norm[mass]_VisualizeError[result,1]")
#hpull_error2 = xframe.pullHist("model_Norm[mass]_VisualizeError[result,2]")
#
#xframe2.addPlotable(hpull_error1, "E_3")
#xframe2.addPlotable(hpull_error2, "E_3")

ROOT.gROOT.SetStyle("Plain")
ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetTitleSize(1, "x")
ROOT.gStyle.SetTitleSize(1, "y")

tl = ROOT.TLatex()
tl.SetTextAlign(12)
tl.SetTextSize(0.04)
tl.DrawLatex(0.1,0.95, "#bf{#it{Private Work}} (CMS simulation)")
if year == "2018":
    tl.DrawLatex(0.57,0.95, "#bf{58.8 fb^{-1}, 2018 (13 TeV)}")
elif year == "2017":
    tl.DrawLatex(0.55,0.95, "41.48 fb^{-1}, 2017 (13 TeV)")
tl.DrawLatex(0.625,0.9, "H #rightarrow aa #rightarrow KKKK")
tl.DrawLatex(0.625,0.85, "Signal MC")
if scope == "ee":
    tl.DrawLatex(0.625,0.8, "ee SR")
elif scope == "mm":
    tl.DrawLatex(0.625,0.8, "#mu#mu SR")
#tl.DrawLatex(0.15,0.8, "#mu = " + str('%.2f'%((muval))) + " #pm " + str('%.2f'%(delmu)))
#tl.DrawLatex(0.15,0.75, "#sigma = " + str('%.2f'%((sigval))) + " #pm " + str('%.2f'%(delsig)))
tl.DrawLatex(0.625,0.75, "#chi^{2}/dof = " + str('%.2f'%((chi2))))

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
xframe2.GetYaxis().SetTitle("#frac{MC - Fit}{Error}")
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

sig.Print("t")


c.SaveAs(f"/web/jhornung/public_html/analysis_plots/sig_fit_{year}_{scope}_bukin_narrow_range.png")
c.SaveAs(f"/web/jhornung/public_html/analysis_plots/sig_fit_{year}_{scope}_bukin_narrow_range.pdf")

#for syst in systs:
#    print(f"Yield variation for {syst.GetName()}: {yield_vars[syst.GetName()]}")
#
for y in yields.keys():
    print(f"Yield for {y}: {yields[y]/yields['nominal']}")

params = result.floatParsFinal()
param_names = [params[i].GetName() for i in range(params.getSize())]
print("Parameter order in covariance matrix:")
print(param_names)

cov = result.covarianceMatrix()
cov.Print()