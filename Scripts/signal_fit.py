import sys
import subprocess
import ROOT
from ROOT import RooFit
import numpy as np
#import matplotlib.pyplot as plt
import utility

def prep_data(scope):
    if scope == "mm":
        df = ROOT.RDataFrame("ntuple", "/ceph/jhornung/MC_2018/2018/HaaKKKK_M125_13TeV_powheg_pythia8-RunIISummer20UL18MiniAODv2-106X/mm/HaaKKKK_M125_13TeV_powheg_pythia8-RunIISummer20UL18MiniAODv2-106X.root")
        df = df.Define("signs", "genWeight > 0 ? 1 : -1")
        df = df.Define("weight", 'signs*id_wgt_mu_1*id_wgt_mu_2*iso_wgt_mu_1*iso_wgt_mu_2*trigger_wgt_mu_1*trigger_wgt_mu_2*puweight')#*evtweight

        
        df = df.Define("ps_mass_mask", "ps_1_mass < 10 && ps_2_mass < 10 && abs(ps_1_mass-ps_2_mass) < 1")
        df = df.Define("mask", "pt_vis > 30 && (m_vis > 75 && m_vis < 105) && abs(H_eta) < 2.4 && d1_pt > 10 && d2_pt > 10 && ps_mass_mask")
        snap = df.Snapshot("ntuple", f"/work/jhornung/Haa/limits/2018/mm/signal_mm.root")
    elif scope == "ee":
        df = ROOT.RDataFrame("ntuple", "/ceph/jhornung/MC_2018/2018/HaaKKKK_M125_13TeV_powheg_pythia8-RunIISummer20UL18MiniAODv2-106X/ee/HaaKKKK_M125_13TeV_powheg_pythia8-RunIISummer20UL18MiniAODv2-106X.root")
        df = df.Define("signs", "genWeight > 0 ? 1 : -1")
        df = df.Define("sfs", 'id_wgt_ele_1*id_wgt_ele_2*trigger_wgt_ele_1*trigger_wgt_ele_2')
        df = df.Define("weight", 'signs*sfs*evtweight*puweight')
        df = df.Define("ps_mass_mask", "ps_1_mass < 10 && ps_2_mass < 10 && abs(ps_1_mass-ps_2_mass) < 1")
        df = df.Define("mask", "pt_vis > 30 && (m_vis > 75 && m_vis < 105) && abs(H_eta) < 2.4 && d1_pt > 10 && d2_pt > 10 && ps_mass_mask")
        snap = df.Snapshot("ntuple", f"/work/jhornung/Haa/limits/2018/ee/signal_ee.root")
    elif scope == "combined":
        df_mm = ROOT.RDataFrame("ntuple", "/ceph/jhornung/MC_2018/2018/HaaKKKK_M125_13TeV_powheg_pythia8-RunIISummer20UL18MiniAODv2-106X/mm/HaaKKKK_M125_13TeV_powheg_pythia8-RunIISummer20UL18MiniAODv2-106X.root")
        df_mm = df_mm.Define("signs", "genWeight > 0 ? 1 : -1")
        df_mm = df_mm.Define("sfs", 'id_wgt_mu_1*id_wgt_mu_2*iso_wgt_mu_1*iso_wgt_mu_2*trigger_wgt_mu_1*trigger_wgt_mu_2')
        df_mm = df_mm.Define("weight", 'signs*sfs*evtweight*puweight')
        df_mm = df_mm.Define("ps_mass_mask", "ps_1_mass < 10 && ps_2_mass < 10 && abs(ps_1_mass-ps_2_mass) < 1")
        df_mm = df_mm.Define("mask", "pt_vis > 30 && (m_vis > 75 && m_vis < 105) && abs(H_eta) < 2.4 && d1_pt > 10 && d2_pt > 10 && ps_mass_mask")
        snap_mm = df_mm.Snapshot("ntuple", f"/work/jhornung/Haa/limits/2018/combined/signal_mm.root")

        df_ee = ROOT.RDataFrame("ntuple", "/ceph/jhornung/MC_2018/2018/HaaKKKK_M125_13TeV_powheg_pythia8-RunIISummer20UL18MiniAODv2-106X/ee/HaaKKKK_M125_13TeV_powheg_pythia8-RunIISummer20UL18MiniAODv2-106X.root")
        df_ee = df_ee.Define("signs", "genWeight > 0 ? 1 : -1")
        df_ee = df_ee.Define("sfs", 'id_wgt_ele_1*id_wgt_ele_2*trigger_wgt_ele_1*trigger_wgt_ele_2')
        df_ee = df_ee.Define("weight", 'signs*sfs*evtweight*puweight')
        df_ee = df_ee.Define("ps_mass_mask", "ps_1_mass < 10 && ps_2_mass < 10 && abs(ps_1_mass-ps_2_mass) < 1")
        df_ee = df_ee.Define("mask", "pt_vis > 30 && (m_vis > 75 && m_vis < 105) && abs(H_eta) < 2.4 && d1_pt > 10 && d2_pt > 10 && ps_mass_mask")
        snap_ee = df_ee.Snapshot("ntuple", f"/work/jhornung/Haa/limits/2018/combined/signal_ee.root")

    return

scope = sys.argv[1]

prep_data(scope)

mass = ROOT.RooRealVar("H_mass", "H_mass", 125, 70, 200)
mask = ROOT.RooRealVar("mask", "mask", 0, 1)
weight = ROOT.RooRealVar("weight", "weight", -10, 10) 

if scope == "mm" or scope == "ee":

    f = ROOT.TFile.Open(f"/work/jhornung/Haa/limits/2018/{scope}/signal_{scope}.root")
    t = f.Get("ntuple")

    data = ROOT.RooDataSet('data', 'data', ROOT.RooArgSet(mass, weight, mask), RooFit.Import(t), RooFit.WeightVar(weight), RooFit.Cut("mask==1"))

elif scope == "combined":

    f_mm = ROOT.TFile.Open(f"/work/jhornung/Haa/limits/2018/combined/signal_mm.root")
    t_mm = f_mm.Get("ntuple")
    data_mm = ROOT.RooDataSet('data_mm', 'data_mm', ROOT.RooArgSet(mass, weight, mask), RooFit.Import(t_mm), RooFit.WeightVar(weight), RooFit.Cut("mask==1"))

    f_ee = ROOT.TFile.Open(f"/work/jhornung/Haa/limits/2018/combined/signal_ee.root")
    t_ee = f_ee.Get("ntuple")

    data_ee = ROOT.RooDataSet('data_ee', 'data_ee', ROOT.RooArgSet(mass, weight, mask), RooFit.Import(t_ee), RooFit.WeightVar(weight), RooFit.Cut("mask==1"))

    data = data_mm
    data.append(data_ee)

data.Print()

MH = ROOT.RooRealVar("MH", "MH", 125.0, 120.0, 130.0)
MH.setConstant(True)
dMH = ROOT.RooRealVar("dMH", "dMH", 0.0, -1.0, 1.0)
mean = ROOT.RooFormulaVar("mean", "mean of gaussians", "@0+@1", ROOT.RooArgList(MH, dMH))
sigmaL = ROOT.RooRealVar("sigma1", "width of gaussians", 1.0, 0.0, 5.0)
sigmaR = ROOT.RooRealVar("sigma2", "width of gaussians", 1.0, 0.0, 5.0)
alphaL = ROOT.RooRealVar("alphaL", "alphaL", 1.0, 0.0, 5.0)
alphaR = ROOT.RooRealVar("alphaR", "alphaR", 1.0, 0.0, 5.0)
nL = ROOT.RooRealVar("nL", "nL", 1.0, 0.0, 5.0)
nR = ROOT.RooRealVar("nR", "nR", 1.0, 0.0, 5.0)
 
sig = ROOT.RooCrystalBall("sig1", "Signal component 1", mass, mean, sigmaL, sigmaR, alphaL, nL, alphaR, nR)

a0 = ROOT.RooRealVar("a0", "a0", 0.5, -1, 1)
a1 = ROOT.RooRealVar("a1", "a1", -0.2, -2, 2.0)
a2 = ROOT.RooRealVar("a2", "a2", .5, -1, 1.0)
a3 = ROOT.RooRealVar("a3", "a3", -.5, -1, 1.0)
a4 = ROOT.RooRealVar("a4", "a4", -.5, -1, 1.0)
comb_offset = ROOT.RooChebychev("offset", "offset", mass, ROOT.RooArgList(a0, a1))

nsig = ROOT.RooRealVar("nsig", "number of signal events", 1.0)
noffset = ROOT.RooRealVar("noffset", "number of background events", 1.5, 0.0, 2.0)
ncomb = ROOT.RooRealVar("ncomb", "number of combinatorial events", 0.6, -0.5, 1)
sigfrac = ROOT.RooFormulaVar("sigfrac", "sigfrac", "@0/(@0+@1)", ROOT.RooArgList(nsig, noffset))
bkgfrac = ROOT.RooFormulaVar("bkgfrac", "bkgfrac", "@1/(@0+@1)", ROOT.RooArgList(nsig, noffset))

model = ROOT.RooAddPdf("model", "(g1)+cbcv+gc", ROOT.RooArgList(sig,comb_offset), ROOT.RooArgList(sigfrac,bkgfrac))

result_one = model.fitTo(data, RooFit.Save(), RooFit.SumW2Error(True), RooFit.Minimizer("Minuit2", "migrad"), RooFit.PrintLevel(-1))

dMH.setConstant(True)
sigmaL.setConstant(True)
sigmaR.setConstant(True)
alphaL.setConstant(True)
alphaR.setConstant(True)
nL.setConstant(True)
nR.setConstant(True)
a0.setConstant(True)
a1.setConstant(True)
nsig.setConstant(True)
noffset.setConstant(True)

n_mc = 355000

n_mc_sel = data.sumEntries()

eff = n_mc_sel/n_mc

xs_Haa = ROOT.RooRealVar("xs_Haa", "xs of Haa [pb]", 1.0) #pb
br_ZH = ROOT.RooRealVar("br_ZH", "br of ZH", 1.0)

eff_signal_H_mass = ROOT.RooRealVar("eff_signal_H_mass", "eff of signal", eff) 

xs_Haa.setConstant(True)
br_ZH.setConstant(True)
eff_signal_H_mass.setConstant(True)

norm_sig = ROOT.RooProduct("model_norm", "sig normalization", ROOT.RooArgList(xs_Haa, br_ZH, eff_signal_H_mass))

f_out = ROOT.TFile(f"/work/jhornung/Haa/limits/2018/{scope}/workspace_sig_{scope}.root", "RECREATE")
w_sig = ROOT.RooWorkspace(f"workspace_sig_weight", f"workspace_sig_weight")
getattr(w_sig, 'import')(model)
getattr(w_sig, 'import')(norm_sig)
w_sig.Print()
w_sig.Write()
f_out.Close()

xframe = mass.frame()
xframe.SetTitle("")
xframe.SetMinimum(0)
#xframe.SetMaximum(7000)

c = ROOT.TCanvas("canvas", "canvas", 1200, 1200)
c.Divide(1,2,0.1,0)

legend = ROOT.TLegend(0.575, 0.4, 0.875, 0.65)

data.plotOn(xframe, RooFit.Binning(65))
hist = xframe.getObject(int(xframe.numItems()-1))
legend.AddEntry(hist, "Data", "ep")

model.plotOn(xframe, ROOT.RooFit.Components("offset"), ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.LineStyle(ROOT.kDashed))
bkg_fit = xframe.getObject(int(xframe.numItems()-1))
legend.AddEntry(bkg_fit, "Combinatorics Bkg", "l")

#model.plotOn(xframe, ROOT.RooFit.Components("combinatorics"), ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.LineStyle(ROOT.kDotted))
#comb_peak = xframe.getObject(int(xframe.numItems()-1))
#legend.AddEntry(comb_peak, "Combinatorics Peak", "l")

#model.plotOn(xframe, ROOT.RooFit.Components("bkg+combinatorics"), ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.LineStyle(ROOT.kDashDotted))
#comb_cont = xframe.getObject(int(xframe.numItems()-1))
#legend.AddEntry(comb_cont, "Combinatorics Cont", "l")

model.plotOn(xframe, RooFit.LineColor(ROOT.kRed), RooFit.LineStyle(ROOT.kSolid))
sig_fit = xframe.getObject(int(xframe.numItems()-1))
legend.AddEntry(sig_fit, "Signal Fit", "l")

data.plotOn(xframe, RooFit.Binning(65))

legend.SetBorderSize(0)
legend.SetFillStyle(0)
legend.SetTextSize(0.04)

num_params = result_one.floatParsFinal().getSize()
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

tl = ROOT.TLatex()
tl.SetTextAlign(12)
tl.SetTextSize(0.04)
tl.DrawLatex(0.1,0.95, "CMS #bf{#it{Simulation Preliminary}}")
tl.DrawLatex(0.7,0.95, "#sqrt{s} = 13 TeV")
tl.DrawLatex(0.575,0.9, "H #rightarrow aa #rightarrow KKKK")
tl.DrawLatex(0.575,0.85, "Signal MC")
#tl.DrawLatex(0.15,0.8, "#mu = " + str('%.2f'%((muval))) + " #pm " + str('%.2f'%(delmu)))
#tl.DrawLatex(0.15,0.75, "#sigma = " + str('%.2f'%((sigval))) + " #pm " + str('%.2f'%(delsig)))
tl.DrawLatex(0.575,0.8, "#chi^{2}/dof = " + str('%.2f'%((chi2))))

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

xframe2.GetXaxis().SetTitle("M_{KKKK} in GeV")
xframe2.GetYaxis().SetTitle("#frac{MC - Fit}{Error}")
xframe2.GetXaxis().SetTitleSize(0.075)
xframe2.GetXaxis().SetTitleOffset(1)
xframe2.GetXaxis().SetLabelSize(0.08)
xframe2.GetYaxis().SetTitleSize(0.075)
xframe2.GetYaxis().SetTitleOffset(0.35)
xframe2.GetYaxis().SetLabelSize(0.08)

xframe2.Draw()

zero = ROOT.TLine(70, 0, 200, 0)
zero.SetLineColor(ROOT.kRed)
zero.SetLineWidth(2)
zero.SetLineStyle(2)
zero.Draw()

model.Print("t")

c.SaveAs(f"/web/jhornung/public_html/signal_fit_{scope}_good_pf_cands.pdf")
c.SaveAs(f"/web/jhornung/public_html/signal_fit_{scope}_good_pf_cands.png")

#print(n_mc)
#print(n_mc_sel)
print(eff)
#print(59e3*eff)

'''
min_mass = mass.getMin()
max_mass = mass.getMax()

masses = np.linspace(min_mass, max_mass, 1000)

fit_vals = np.zeros(masses.shape[0])
comb_offset_fit_vals = np.zeros(masses.shape[0])
combinatorics_fit_vals = np.zeros(masses.shape[0])

for i, m in enumerate(masses):
    mass.setVal(m)
    fit_vals[i] = model.getVal(ROOT.RooArgSet(mass))
    comb_offset_fit_vals[i] = comb_offset.getVal(ROOT.RooArgSet(mass))
    combinatorics_fit_vals[i] = combinatorics.getVal(ROOT.RooArgSet(mass))

hist_file = uproot.open(f"/work/jhornung/Haa/limits/signal_{scope}.root")
ntuple = hist_file["ntuple"]
branches = ntuple.arrays()

fig, ax = plt.subplots()
ax.plot(masses, fit_vals, label="Fit")
ax.plot(masses, comb_offset_fit_vals, label="Combinatorics Offset", linestyle="--")
ax.plot(masses, combinatorics_fit_vals, label="Combinatorics Peak", linestyle="--")
ax.hist(branches["H_mass"][branches["mask"] == 1], bins=np.linspace(min_mass, max_mass, 66), density=True, histtype='step', weights=branches["weight"][branches["mask"] == 1], label="Data")

ax.legend()

ax.set_xlim(min_mass, max_mass)

plt.show()
'''