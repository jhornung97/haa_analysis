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
        df = df.Define("iso_mask", "d1_iso < 3 && d2_iso < 3 && d3_iso < 3 && d4_iso < 3")
        df = df.Define("mask", "pt_vis > 30 && (m_vis > 75 && m_vis < 105) && abs(H_eta) < 2.4 && d1_pt > 10 && d2_pt > 10 && ps_mass_mask && iso_mask")
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
        df = df.Define("iso_mask", "d1_iso < 3 && d2_iso < 3 && d3_iso < 3 && d4_iso < 3")
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

mc, yields, yield_vars, systs, mass = prep_data(scope, year)


MH = ROOT.RooRealVar(f"{scope}_MH", "MH", 125.0, 120.0, 130.0)
#MH.setConstant(True)
#dMH = ROOT.RooRealVar(f"{year}_{scope}_dMH", "dMH", 0.0, -1.0, 1.0)
#mean = ROOT.RooFormulaVar(f"{year}_{scope}_mean", "mean of gaussians", "@0+@1", ROOT.RooArgList(MH, dMH))
sigmaL = ROOT.RooRealVar(f"{scope}_sigmaL", "width of gaussians", 1.0, 0.5, 1.5)
sigmaR = ROOT.RooRealVar(f"{scope}_sigmaR", "width of gaussians", 1.0, 0.7, 1.3)
alphaL = ROOT.RooRealVar(f"{scope}_alphaL", "alphaL", 2, 1.5, 2.5)
alphaR = ROOT.RooRealVar(f"{scope}_alphaR", "alphaR", 1.5, 0.4, 2.7)
nL = ROOT.RooRealVar(f"{scope}_nL", "nL", 0.9, 0.3, 1.5)
nR = ROOT.RooRealVar(f"{scope}_nR", "nR", 3.7, 3.4, 4)
 
sig = ROOT.RooCrystalBall(f"{scope}_sig", "Signal component 1", mass, MH, sigmaL, sigmaR, alphaL, nL, alphaR, nR)

model = sig

a0 = ROOT.RooRealVar(f"{scope}_a0", "a0", -1, -1.54, -0.46)
a1 = ROOT.RooRealVar(f"{scope}_a1", "a1", 0.2, -.28, .68)
a2 = ROOT.RooRealVar(f"{scope}_a2", "a2", 0.1, -1, 1)

comb_offset = ROOT.RooChebychev(f"{scope}_offset", "offset", mass, ROOT.RooArgList(a0, a1))

frac = ROOT.RooRealVar(f"{scope}_sigfrac", "sigfrac", 0.5, 0.0, 1.0)

model = ROOT.RooAddPdf(f"{scope}_model", "cb+cbcv", ROOT.RooArgList(sig,comb_offset), ROOT.RooArgList(frac))

mean_vals = {}
sigmaL_vals = {}
sigmaR_vals = {}
alphaL_vals = {}
alphaR_vals = {}
nL_vals = {}
nR_vals = {}
a0_vals = {}
a1_vals = {}

model.fitTo(mc["nominal"], RooFit.Save(), RooFit.SumW2Error(False), RooFit.Minimizer("Minuit2", "migrad"), RooFit.PrintLevel(-1))
model.fitTo(mc["nominal"], RooFit.Save(), RooFit.SumW2Error(False), RooFit.Minimizer("Minuit2", "migrad"), RooFit.PrintLevel(-1))

mean_vals["nominal"] = [MH.getVal(), MH.getError()]
sigmaL_vals["nominal"] = [sigmaL.getVal(), sigmaL.getError()]
sigmaR_vals["nominal"] = [sigmaR.getVal(), sigmaR.getError()]
alphaL_vals["nominal"] = [alphaL.getVal(), alphaL.getError()]
alphaR_vals["nominal"] = [alphaR.getVal(), alphaR.getError()]
nL_vals["nominal"] = [nL.getVal(), nL.getError()]
nR_vals["nominal"] = [nR.getVal(), nR.getError()]
a0_vals["nominal"] = [a0.getVal(), a0.getError()]
a1_vals["nominal"] = [a1.getVal(), a1.getError()]


for syst in systs:
    model.fitTo(mc[syst.GetName()], RooFit.Save(), RooFit.SumW2Error(False), RooFit.Minimizer("Minuit2", "migrad"), RooFit.PrintLevel(-1))
    model.fitTo(mc[syst.GetName()], RooFit.Save(), RooFit.SumW2Error(False), RooFit.Minimizer("Minuit2", "migrad"), RooFit.PrintLevel(-1))

    mean_vals[syst.GetName()] = [MH.getVal(), MH.getError()]
    sigmaL_vals[syst.GetName()] = [sigmaL.getVal(), sigmaL.getError()]
    sigmaR_vals[syst.GetName()] = [sigmaR.getVal(), sigmaR.getError()]
    alphaL_vals[syst.GetName()] = [alphaL.getVal(), alphaL.getError()]
    alphaR_vals[syst.GetName()] = [alphaR.getVal(), alphaR.getError()]
    nL_vals[syst.GetName()] = [nL.getVal(), nL.getError()]
    nR_vals[syst.GetName()] = [nR.getVal(), nR.getError()]
    a0_vals[syst.GetName()] = [a0.getVal(), a0.getError()]
    a1_vals[syst.GetName()] = [a1.getVal(), a1.getError()]

mean_variations = []

for syst in systs:
    mean_var = np.abs(mean_vals[syst.GetName()][0] - mean_vals["nominal"][0])/mean_vals["nominal"][0]
    mean_variations.append(mean_var)

max_mean_var = np.max(mean_variations)

MH.setConstant(True)
dMH = ROOT.RooRealVar(f"{scope}_dMH", "dMH", 0.0, -1.0, 1.0)
eta = ROOT.RooRealVar(f"{scope}_nusiance_scale", "nusiances_scale", 0.0, -5.0, 5.0)
eta.setConstant(True)
mean_formula = ROOT.RooFormulaVar(f"{scope}_mean_formula", "mean of gaussians", f"(@0+@1)*(1+{max_mean_var}*@2)", ROOT.RooArgList(MH, dMH, eta))

sigmaL_variations = []

print(f"sigmaL_nominal = {sigmaL_vals['nominal'][0]} +- {sigmaL_vals['nominal'][1]}")

for syst in systs:
    print(f"sigmaL_{syst.GetName()} = {sigmaL_vals[syst.GetName()][0]} +- {sigmaL_vals[syst.GetName()][1]}")
    sigmaL_var = np.abs(sigmaL_vals[syst.GetName()][0] - sigmaL_vals["nominal"][0])/sigmaL_vals["nominal"][0]
    sigmaL_variations.append(sigmaL_var)

max_sigmaL_var = np.max(sigmaL_variations)

sigmaL_nuis = ROOT.RooRealVar(f"{scope}_sigmaL_nuis", "width of gaussians", 0.0, -5.0, 5.0)
sigmaL_nuis.setConstant(True)
sigmaL_formula  = ROOT.RooFormulaVar(f"{scope}_sigmaL_formula", "width of gaussians", f"@0*(1+{max_sigmaL_var}*@1)", ROOT.RooArgList(sigmaL, sigmaL_nuis))

sigmaR_variations = []

for syst in systs:
    sigmaR_var = np.abs(sigmaR_vals[syst.GetName()][0] - sigmaR_vals["nominal"][0])/sigmaR_vals["nominal"][0]
    sigmaR_variations.append(sigmaR_var)

max_sigmaR_var = np.max(sigmaR_variations)

sigmaR_nuis = ROOT.RooRealVar(f"{scope}_sigmaR_nuis", "width of gaussians", 0.0, -5.0, 5.0)
sigmaR_nuis.setConstant(True)
sigmaR_formula  = ROOT.RooFormulaVar(f"{scope}_sigmaR_formula", "width of gaussians", f"@0*(1+{max_sigmaR_var}*@1)", ROOT.RooArgList(sigmaR, sigmaR_nuis))

alphaR_variations = []
alphaL_variations = []
nL_variations = []
nR_variations = []
a0_variations = []
a1_variations = []

for syst in systs:
    alphaL_var = np.abs(alphaL_vals[syst.GetName()][0] - alphaL_vals["nominal"][0])/alphaL_vals["nominal"][0]
    alphaR_var = np.abs(alphaR_vals[syst.GetName()][0] - alphaR_vals["nominal"][0])/alphaR_vals["nominal"][0]
    nL_var = np.abs(nL_vals[syst.GetName()][0] - nL_vals["nominal"][0])/nL_vals["nominal"][0]
    nR_var = np.abs(nR_vals[syst.GetName()][0] - nR_vals["nominal"][0])/nR_vals["nominal"][0]
    a0_var = np.abs((a0_vals[syst.GetName()][0] - a0_vals["nominal"][0])/a0_vals["nominal"][0])
    a1_var = np.abs(a1_vals[syst.GetName()][0] - a1_vals["nominal"][0])/a1_vals["nominal"][0]

    alphaL_variations.append(alphaL_var)
    alphaR_variations.append(alphaR_var)
    nL_variations.append(nL_var)
    nR_variations.append(nR_var)
    a0_variations.append(a0_var)
    a1_variations.append(a1_var)

max_alphaL_var = np.max(alphaL_variations)
max_alphaR_var = np.max(alphaR_variations)
max_nL_var = np.max(nL_variations)
max_nR_var = np.max(nR_variations)
max_a0_var = np.max(a0_variations)
max_a1_var = np.max(a1_variations)

alphaL_nuis = ROOT.RooRealVar(f"{scope}_alphaL_nuis", "alphaL", 0.0, -5.0, 5.0)
alphaL_nuis.setConstant(True)
alphaL_formula  = ROOT.RooFormulaVar(f"{scope}_alphaL_formula", "alphaL", f"@0*(1+{max_alphaL_var}*@1)", ROOT.RooArgList(alphaL, alphaL_nuis))

alphaR_nuis = ROOT.RooRealVar(f"{scope}_alphaR_nuis", "alphaR", 0.0, -5.0, 5.0)
alphaR_nuis.setConstant(True)
alphaR_formula  = ROOT.RooFormulaVar(f"{scope}_alphaR_formula", "alphaR", f"@0*(1+{max_alphaR_var}*@1)", ROOT.RooArgList(alphaR, alphaR_nuis))

nL_nuis = ROOT.RooRealVar(f"{scope}_nL_nuis", "nL", 0.0, -5.0, 5.0)
nL_nuis.setConstant(True)
nL_formula  = ROOT.RooFormulaVar(f"{scope}_nL_formula", "nL", f"@0*(1+{max_nL_var}*@1)", ROOT.RooArgList(nL, nL_nuis))

nR_nuis = ROOT.RooRealVar(f"{scope}_nR_nuis", "nR", 0.0, -5.0, 5.0)
nR_nuis.setConstant(True)
nR_formula  = ROOT.RooFormulaVar(f"{scope}_nR_formula", "nR", f"@0*(1+{max_nR_var}*@1)", ROOT.RooArgList(nR, nR_nuis))

a0_nuis = ROOT.RooRealVar(f"{scope}_a0_nuis", "a0", 0.0, -5.0, 5.0)
a0_nuis.setConstant(True)
a0_formula  = ROOT.RooFormulaVar(f"{scope}_a0_formula", "a0", f"@0*(1+{max_a0_var}*@1)", ROOT.RooArgList(a0, a0_nuis))
#
a1_nuis = ROOT.RooRealVar(f"{scope}_a1_nuis", "a1", 0.0, -5.0, 5.0)
a1_nuis.setConstant(True)
a1_formula  = ROOT.RooFormulaVar(f"{scope}_a1_formula", "a1", f"@0*(1+{max_a1_var}*@1)", ROOT.RooArgList(a1, a1_nuis))

sig = ROOT.RooCrystalBall(f"{scope}_sig1", "Signal component 1", mass, mean_formula, sigmaL_formula, sigmaR_formula, alphaL_formula, nL_formula, alphaR_formula, nR_formula)

comb_offset = ROOT.RooChebychev(f"{scope}_offset", "offset", mass, ROOT.RooArgList(a0_formula, a1_formula))

model_nom = ROOT.RooAddPdf(f"model_Haa_H_mass_{year}_{scope}", "cb+cbcv", ROOT.RooArgList(sig,comb_offset), ROOT.RooArgList(frac))

result = model_nom.fitTo(mc["nominal"], RooFit.Save(), RooFit.SumW2Error(False), RooFit.Minimizer("Minuit2", "migrad"), RooFit.PrintLevel(-1))

#xs_ZH = ROOT.RooRealVar(f"xs_ZH", "xs of ZH [pb]", 1.0) #pb
#br_Haa = ROOT.RooRealVar(f"br_Haa", "br of Haa", 1.0)

eff_nom = float("%.2f"%(mc["nominal"].sumEntries()))
#eff_nom_signal_H_mass = ROOT.RooRealVar(f"{year}_{scope}_eff_signal_H_mass", "eff of signal", eff_nom)
#
#
#xs_ZH.setConstant(True)
#br_Haa.setConstant(True)
#eff_nom_signal_H_mass.setConstant(True)
#
#r = ROOT.RooRealVar(f"r", "r", 1.0, 0.0, 10.0)
#norm_sig_nom = ROOT.RooProduct(f"model_signal_H_mass_{year}_{scope}_norm", "sig normalization", ROOT.RooArgList(xs_ZH, br_Haa, eff_nom_signal_H_mass, r))

dMH.setConstant(True)

sigmaL.setConstant(True)
sigmaR.setConstant(True)
alphaL.setConstant(True)
alphaR.setConstant(True)
nL.setConstant(True)
nR.setConstant(True)
a0.setConstant(True)
a1.setConstant(True)
frac.setConstant(True)

f_out = ROOT.TFile(f"/work/jhornung/Haa/limits/{year}/{scope}/workspace_sig_{scope}.root", "RECREATE")
w_sig_nom = ROOT.RooWorkspace(f"workspace_sig_nom", f"workspace_sig")
getattr(w_sig_nom, 'import')(model_nom)
#getattr(w_sig_nom, 'import')(norm_sig_nom)
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

#for syst in systs:
#    print(f"Yield variation for {syst.GetName()}: {yield_vars[syst.GetName()]}")
#
for y in yields.keys():
    print(f"Yield for {y}: {yields[y]/yields['nominal']}")

print(f"Max mean variation: {max_mean_var}")
print(f"Max sigmaL variation: {max_sigmaL_var}")
print(f"Max sigmaR variation: {max_sigmaR_var}")
print(f"Max alphaL variation: {max_alphaL_var}")
print(f"Max alphaR variation: {max_alphaR_var}")
print(f"Max nL variation: {max_nL_var}")
print(f"Max nR variation: {max_nR_var}")
print(f"Max a0 variation: {max_a0_var}")
print(f"Max a1 variation: {max_a1_var}")

params = result.floatParsFinal()
param_names = [params[i].GetName() for i in range(params.getSize())]
print("Parameter order in covariance matrix:")
print(param_names)

cov = result.covarianceMatrix()
cov.Print()

print(eff_nom)