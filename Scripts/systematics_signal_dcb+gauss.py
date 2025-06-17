import sys
import subprocess
import ROOT
from ROOT import RooFit
import numpy as np
import matplotlib.pyplot as plt
import utility
from array import array

def prep_data(scope, year, iso_mask):
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
        df = df.Define("iso_mask", f"d1_iso < {iso_mask} && d2_iso < {iso_mask}") # && d3_iso < {iso_mask} && d4_iso < {iso_mask}")
        df = df.Define("mask", "pt_vis > 30 && (m_vis > 75 && m_vis < 105) && abs(H_eta) < 2.4 && d1_pt > 10 && d2_pt > 10 && ps_mass_mask") #&& iso_mask")
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
        df = df.Define("iso_mask", f"d1_iso < {iso_mask} && d2_iso < {iso_mask}") # && d3_iso < {iso_mask} && d4_iso < {iso_mask}")
        df = df.Define("mask", "pt_vis > 30 && (m_vis > 75 && m_vis < 105) && abs(H_eta) < 2.4 && d1_pt > 10 && d2_pt > 10 && ps_mass_mask && iso_mask")
        snap = df.Snapshot("ntuple", f"/work/jhornung/Haa/limits/{year}/{scope}/signal_ee.root")

    mass = ROOT.RooRealVar(f"H_mass_{year}_{scope}", "H_mass", 125, 115, 135)

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
        yields = {}
        yield_vars = {}
        

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
        yields = {}
        yield_vars = {}
        

        mc["nominal"] = ROOT.RooDataSet('data_nominal', 'data_nominal', ROOT.RooArgSet(mass, weight_nominal, mask), RooFit.Import(t), RooFit.WeightVar(weight_nominal), RooFit.Cut("mask==1"))
        yields["nominal"] = mc["nominal"].sumEntries()
        systs = [weight_pu_up, weight_pu_down, weight_id_1_up, weight_id_1_down, weight_id_2_up, weight_id_2_down, weight_trigger_1_up, weight_trigger_1_down, weight_trigger_2_up, weight_trigger_2_down]

        for syst in systs:
            mc[syst.GetName()] = ROOT.RooDataSet(f'data_{syst.GetName()}', f'data_{syst.GetName()}', ROOT.RooArgSet(mass, syst, mask), RooFit.Import(t), RooFit.WeightVar(syst), RooFit.Cut("mask==1"))
            yields[syst.GetName()] = mc[syst.GetName()].sumEntries()

            yield_var = mc[syst.GetName()].sumEntries()/mc["nominal"].sumEntries()
            yield_vars[syst.GetName()] = yield_var

    return mc, yields, yield_vars, systs, mass

scope = sys.argv[1]
year = sys.argv[2]
iso_mask = float(sys.argv[3])

mc, yields, yield_vars, systs, mass = prep_data(scope, year, iso_mask) 


MH = ROOT.RooRealVar(f"{scope}_MH", "MH", 125.0, 120.0, 130.0)
sigma_g = ROOT.RooRealVar(f"{scope}_sigma", "width of gaussian", 1.0, 0.5, 1.5)
 
gauss = ROOT.RooGaussian(f"{scope}_sig", "gauss", mass, MH, sigma_g)

sigma_cb = ROOT.RooRealVar(f"{scope}_sigma_cb", "width of cb", 5, 0, 10)
alpha = ROOT.RooRealVar(f"{scope}_alpha", "alpha", 1.0, 0.0, 5.0)
n = ROOT.RooRealVar(f"{scope}_n", "n", 1.0, 0.0, 10.0)

cb = ROOT.RooCrystalBall(f"{scope}_sig_cb", "Crystal Ball", mass, MH, sigma_cb, alpha, n)

frac = ROOT.RooRealVar(f"{scope}_sigfrac", "sigfrac", 0.5, 0.0, 1.0)

model = ROOT.RooAddPdf(f"{scope}_model", "gauss+cb", ROOT.RooArgList(gauss,cb), ROOT.RooArgList(frac))

mean_vals = {}
sigma_g_vals = {}
sigma_cb_vals = {}
alpha_vals = {}
n_vals = {}

model.fitTo(mc["nominal"], RooFit.Save(), RooFit.SumW2Error(False), RooFit.Minimizer("Minuit2", "migrad"), RooFit.PrintLevel(-1))
model.fitTo(mc["nominal"], RooFit.Save(), RooFit.SumW2Error(False), RooFit.Minimizer("Minuit2", "migrad"), RooFit.PrintLevel(-1))

mean_vals["nominal"] = [MH.getVal(), MH.getError()]
sigma_g_vals["nominal"] = [sigma_g.getVal(), sigma_g.getError()]
sigma_cb_vals["nominal"] = [sigma_cb.getVal(), sigma_cb.getError()]
alpha_vals["nominal"] = [alpha.getVal(), alpha.getError()]
n_vals["nominal"] = [n.getVal(), n.getError()]


for syst in systs:
    model.fitTo(mc[syst.GetName()], RooFit.Save(), RooFit.SumW2Error(False), RooFit.Minimizer("Minuit2", "migrad"), RooFit.PrintLevel(-1))
    model.fitTo(mc[syst.GetName()], RooFit.Save(), RooFit.SumW2Error(False), RooFit.Minimizer("Minuit2", "migrad"), RooFit.PrintLevel(-1))

    mean_vals[syst.GetName()] = [MH.getVal(), MH.getError()]
    sigma_g_vals[syst.GetName()] = [sigma_g.getVal(), sigma_g.getError()]
    sigma_cb_vals[syst.GetName()] = [sigma_cb.getVal(), sigma_cb.getError()]
    alpha_vals[syst.GetName()] = [alpha.getVal(), alpha.getError()]
    n_vals[syst.GetName()] = [n.getVal(), n.getError()]


mean_variations = []
sigma_g_variations = []
sigma_cb_variations = []
alpha_variations = []
n_variations = []

for syst in systs:
    mean_var = np.abs(mean_vals[syst.GetName()][0] - mean_vals["nominal"][0])/mean_vals["nominal"][0]
    sigma_g_var = np.abs(sigma_g_vals[syst.GetName()][0] - sigma_g_vals["nominal"][0])/sigma_g_vals["nominal"][0]
    sigma_cb_var = np.abs(sigma_cb_vals[syst.GetName()][0] - sigma_cb_vals["nominal"][0])/sigma_cb_vals["nominal"][0]
    alpha_var = np.abs(alpha_vals[syst.GetName()][0] - alpha_vals["nominal"][0])/alpha_vals["nominal"][0]
    n_var = np.abs(n_vals[syst.GetName()][0] - n_vals["nominal"][0])/n_vals["nominal"][0]

    mean_variations.append(mean_var)
    sigma_g_variations.append(sigma_g_var)
    sigma_cb_variations.append(sigma_cb_var)
    alpha_variations.append(alpha_var)
    n_variations.append(n_var)

max_mean_var = np.max(mean_variations)
max_sigma_g_var = np.max(sigma_g_variations)
max_sigma_cb_var = np.max(sigma_cb_variations)
max_alpha_var = np.max(alpha_variations)
max_n_var = np.max(n_variations)

mean_nuisance = ROOT.RooRealVar(f"{scope}_mean_nuisance", "mean nuisance", 0, -1, 1)
sigma_g_nuisance = ROOT.RooRealVar(f"{scope}_sigma_g_nuisance", "sigma_g nuisance", 0, -1, 1)
sigma_cb_nuisance = ROOT.RooRealVar(f"{scope}_sigma_cb_nuisance", "sigma_cb nuisance", 0, -1, 1)
alpha_nuisance = ROOT.RooRealVar(f"{scope}_alpha_nuisance", "alpha nuisance", 0, -1, 1)
n_nuisance = ROOT.RooRealVar(f"{scope}_n_nuisance", "n nuisance", 0, -1, 1)

mean_nuisance.setConstant(True)
sigma_g_nuisance.setConstant(True)
sigma_cb_nuisance.setConstant(True)
alpha_nuisance.setConstant(True)
n_nuisance.setConstant(True)  

mean_formula = ROOT.RooFormulaVar(f"{scope}_mean_formula", f"mean formula", f"@0*(1+{max_mean_var}*@1)", ROOT.RooArgList(MH, mean_nuisance))
sigma_g_formula = ROOT.RooFormulaVar(f"{scope}_sigma_g_formula", f"sigma_g formula", f"@0*(1+{max_sigma_g_var}*@1)", ROOT.RooArgList(sigma_g, sigma_g_nuisance))
sigma_cb_formula = ROOT.RooFormulaVar(f"{scope}_sigma_cb_formula", f"sigma_cb formula", f"@0*(1+{max_sigma_cb_var}*@1)", ROOT.RooArgList(sigma_cb, sigma_cb_nuisance))
alpha_formula = ROOT.RooFormulaVar(f"{scope}_alpha_formula", f"alpha formula", f"@0*(1+{max_alpha_var}*@1)", ROOT.RooArgList(alpha, alpha_nuisance))
n_formula = ROOT.RooFormulaVar(f"{scope}_n_formula", f"n formula", f"@0*(1+{max_n_var}*@1)", ROOT.RooArgList(n, n_nuisance))

gauss = ROOT.RooGaussian(f"{scope}_gauss", "gauss", mass, mean_formula, sigma_g_formula)
cb = ROOT.RooCrystalBall(f"{scope}_cb", "Crystal Ball", mass, mean_formula, sigma_cb_formula, alpha_formula, n_formula)

model_nom = ROOT.RooAddPdf(f"model_Haa_H_mass_{year}_{scope}", "gauss+cb", ROOT.RooArgList(gauss, cb), ROOT.RooArgList(frac))
result = model_nom.fitTo(mc["nominal"], RooFit.Save(), RooFit.SumW2Error(False), RooFit.Minimizer("Minuit2", "migrad"), RooFit.PrintLevel(-1))

eff_nom = float("%.2f"%(mc["nominal"].sumEntries()))
eff_nom = ROOT.RooRealVar(f"model_Haa_H_mass_{year}_{scope}_norm", "eff_nom", eff_nom)

MH.setConstant(True)
sigma_g.setConstant(True)
sigma_cb.setConstant(True)
alpha.setConstant(True)
n.setConstant(True)
frac.setConstant(True)

f_out = ROOT.TFile(f"/work/jhornung/Haa/limits/{year}/{scope}/workspace_sig_{scope}.root", "RECREATE")
w_sig_nom = ROOT.RooWorkspace(f"workspace_sig_nom", f"workspace_sig")
getattr(w_sig_nom, 'import')(model_nom)
getattr(w_sig_nom, 'import')(eff_nom)
w_sig_nom.Print()
w_sig_nom.Write()
f_out.Close()

c = ROOT.TCanvas("canvas", "canvas", 1200, 1200)
c.Divide(1,2,0.1,0)

legend = ROOT.TLegend(0.625, 0.35, 0.875, 0.6)

xframe = mass.frame()

mc["nominal"].plotOn(xframe)

model_nom.plotOn(xframe, RooFit.LineColor(ROOT.kRed), RooFit.LineStyle(ROOT.kSolid))
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
tl.DrawLatex(0.625,0.75, "#chi^{2}/dof = " + str('%.2f'%((chi2))))

c.cd(1)
c.cd(1).SetTopMargin(0.1)
c.cd(1).SetRightMargin(0.1)
c.cd(1).SetBottomMargin(0.001)
c.cd(1).SetLeftMargin(0.1)
ROOT.gPad.SetBottomMargin(0.001)
c.cd(1).SetPad(0,0.3,1,1)

xframe.GetXaxis().SetLabelSize(0.05)
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

model_nom.Print("t")


c.SaveAs(f"/web/jhornung/public_html/analysis_plots/sig_fit_{year}_{scope}_test.png")
c.SaveAs(f"/web/jhornung/public_html/analysis_plots/sig_fit_{year}_{scope}_test.pdf")

for y in yields.keys():
    print(f"Yield for {y}: {yields[y]/yields['nominal']}")

print(f"Max mean variation: {max_mean_var}")
print(f"Max sigma_g variation: {max_sigma_g_var}")
print(f"Max sigma_cb variation: {max_sigma_cb_var}")
print(f"Max alpha variation: {max_alpha_var}")
print(f"Max n variation: {max_n_var}")

params = result.floatParsFinal()
param_names = [params[i].GetName() for i in range(params.getSize())]
print("Parameter order in covariance matrix:")
print(param_names)

cov = result.covarianceMatrix()
cov.Print()

print(eff_nom)