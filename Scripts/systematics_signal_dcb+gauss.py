import sys
import subprocess
import ROOT
from ROOT import RooFit
import numpy as np
import matplotlib.pyplot as plt
import utility
from array import array
import glob
import os
import json

def prep_data(scope, year, ma_value, signal_path):
    lumi = 0
    if year == "2018":
        lumi = 59.83e3
    elif year == "2017":
        lumi = 41.48e3
    elif year == "2016preVFP":
        lumi = 19.5e3
    elif year == "2016postVFP":
        lumi = 16.8e3
    elif year == "2024":
        lumi = 109.95e3
    if scope == "mm":
        df = ROOT.RDataFrame("ntuple", f"{signal_path}/mm/{os.path.basename(signal_path)}.root")
        df = df.Define(f"H_mass_{year}_{scope}", "H_mass")
        df = df.Define("signs", "genWeight > 0 ? 1 : -1")
        df = df.Define("trigger_sfs", "1-(1 - trigger_wgt_mu_1)*(1 - trigger_wgt_mu_2)")
        df = df.Define("trigger_sfs_up_1", "1-(1 - trigger_wgt_mu_1__MuonTriggerUp)*(1 - trigger_wgt_mu_2)")
        df = df.Define("trigger_sfs_down_1", "1-(1 - trigger_wgt_mu_1__MuonTriggerDown)*(1 - trigger_wgt_mu_2)")
        df = df.Define("trigger_sfs_up_2", "1-(1 - trigger_wgt_mu_1)*(1 - trigger_wgt_mu_2__MuonTriggerUp)")
        df = df.Define("trigger_sfs_down_2", "1-(1 - trigger_wgt_mu_1)*(1 - trigger_wgt_mu_2__MuonTriggerDown)")
        df = df.Define("nominal", f'{lumi}*signs*evtweight*id_wgt_mu_1*id_wgt_mu_2*iso_wgt_mu_1*iso_wgt_mu_2*trigger_sfs*puweight')#*evtweight
        df = df.Define("weight_pu_up", f'{lumi}*signs*evtweight*id_wgt_mu_1*id_wgt_mu_2*iso_wgt_mu_1*iso_wgt_mu_2*trigger_sfs*puweight__PUUp')#*evtweight
        df = df.Define("weight_pu_down", f'{lumi}*signs*evtweight*id_wgt_mu_1*id_wgt_mu_2*iso_wgt_mu_1*iso_wgt_mu_2*trigger_sfs*puweight__PUDown')#*evtweight
        df = df.Define("weight_id_1_up", f'{lumi}*signs*evtweight*id_wgt_mu_1__MuonIDUp*id_wgt_mu_2*iso_wgt_mu_1*iso_wgt_mu_2*trigger_sfs*puweight')#*evtweight
        df = df.Define("weight_id_1_down", f'{lumi}*signs*evtweight*id_wgt_mu_1__MuonIDDown*id_wgt_mu_2*iso_wgt_mu_1*iso_wgt_mu_2*trigger_sfs*puweight')#*evtweight
        df = df.Define("weight_id_2_up", f'{lumi}*signs*evtweight*id_wgt_mu_1*id_wgt_mu_2__MuonIDUp*iso_wgt_mu_1*iso_wgt_mu_2*trigger_sfs*puweight')#*evtweight
        df = df.Define("weight_id_2_down", f'{lumi}*signs*evtweight*id_wgt_mu_1*id_wgt_mu_2__MuonIDDown*iso_wgt_mu_1*iso_wgt_mu_2*trigger_sfs*puweight')#*evtweight
        df = df.Define("weight_iso_1_up", f'{lumi}*signs*evtweight*id_wgt_mu_1*id_wgt_mu_2*iso_wgt_mu_1__MuonIsoUp*iso_wgt_mu_2*trigger_sfs*puweight')#*evtweight
        df = df.Define("weight_iso_1_down", f'{lumi}*signs*evtweight*id_wgt_mu_1*id_wgt_mu_2*iso_wgt_mu_1__MuonIsoDown*iso_wgt_mu_2*trigger_sfs*puweight')#*evtweight
        df = df.Define("weight_iso_2_up", f'{lumi}*signs*evtweight*id_wgt_mu_1*id_wgt_mu_2*iso_wgt_mu_1*iso_wgt_mu_2__MuonIsoUp*trigger_sfs*puweight')#*evtweight
        df = df.Define("weight_iso_2_down", f'{lumi}*signs*evtweight*id_wgt_mu_1*id_wgt_mu_2*iso_wgt_mu_1*iso_wgt_mu_2__MuonIsoDown*trigger_sfs*puweight')#*evtweight
        df = df.Define("weight_trigger_1_up", f'{lumi}*signs*evtweight*id_wgt_mu_1*id_wgt_mu_2*iso_wgt_mu_1*iso_wgt_mu_2*trigger_sfs_up_1*puweight')#*evtweight
        df = df.Define("weight_trigger_1_down", f'{lumi}*signs*evtweight*id_wgt_mu_1*id_wgt_mu_2*iso_wgt_mu_1*iso_wgt_mu_2*trigger_sfs_down_1*puweight')#*evtweight
        df = df.Define("weight_trigger_2_up", f'{lumi}*signs*evtweight*id_wgt_mu_1*id_wgt_mu_2*iso_wgt_mu_1*iso_wgt_mu_2*trigger_sfs_up_2*puweight')#*evtweight
        df = df.Define("weight_trigger_2_down", f'{lumi}*signs*evtweight*id_wgt_mu_1*id_wgt_mu_2*iso_wgt_mu_1*iso_wgt_mu_2*trigger_sfs_down_2*puweight')#*evtweight
        df = df.Define("ps_mass_mask", "ps_1_mass < 3 && ps_2_mass < 3 && abs(ps_1_mass-ps_2_mass) < 0.06")
        df = df.Define("iso_mask", f"d1_iso < 0.5 && d2_iso < 0.5") # && d3_iso < {iso_mask} && d4_iso < {iso_mask}")
        df = df.Define("mask", "pt_vis > 30 && (m_vis > 75 && m_vis < 105) && abs(H_eta) < 2.4 && d1_pt > 10 && d2_pt > 10 && ps_mass_mask && iso_mask")
        snap = df.Snapshot("ntuple", f"/work/jhornung/Haa/limits/{year}/{scope}/signal_mm_ma{ma_value}.root")
    elif scope == "ee":
        df = ROOT.RDataFrame("ntuple", f"{signal_path}/ee/{os.path.basename(signal_path)}.root")
        df = df.Define(f"H_mass_{year}_{scope}", "H_mass")
        df = df.Define("signs", "genWeight > 0 ? 1 : -1")
        #df = df.Define("trigger_sfs", "1-(1 - trigger_wgt_ele_1)*(1 - trigger_wgt_ele_2)")
        #df = df.Define("trigger_sfs_up_1", "1-(1 - trigger_wgt_ele_1_up)*(1 - trigger_wgt_ele_2)")
        #df = df.Define("trigger_sfs_down_1", "1-(1 - trigger_wgt_ele_1_down)*(1 - trigger_wgt_ele_2)")
        #df = df.Define("trigger_sfs_up_2", "1-(1 - trigger_wgt_ele_1)*(1 - trigger_wgt_ele_2_up)")
        #df = df.Define("trigger_sfs_down_2", "1-(1 - trigger_wgt_ele_1)*(1 - trigger_wgt_ele_2_down)")
        df = df.Define("sfs", 'id_wgt_ele_1*id_wgt_ele_2')#*trigger_sfs')
        df = df.Define("nominal", f'{lumi}*signs*evtweight*sfs*puweight')
        df = df.Define("weight_pu_up", f'{lumi}*signs*evtweight*sfs*puweight__PUUp')
        df = df.Define("weight_pu_down", f'{lumi}*signs*evtweight*sfs*puweight__PUDown')
        df = df.Define("weight_id_1_up", f'{lumi}*signs*evtweight*id_wgt_ele_1__ElectronIDUp*id_wgt_ele_2*1*puweight')
        df = df.Define("weight_id_1_down", f'{lumi}*signs*evtweight*id_wgt_ele_1__ElectronIDDown*id_wgt_ele_2*1*puweight')
        df = df.Define("weight_id_2_up", f'{lumi}*signs*evtweight*id_wgt_ele_1*id_wgt_ele_2__ElectronIDUp*1*puweight')
        df = df.Define("weight_id_2_down", f'{lumi}*signs*evtweight*id_wgt_ele_1*id_wgt_ele_2__ElectronIDDown*1*puweight')
        df = df.Define("weight_trigger_1_up", f'{lumi}*signs*evtweight*id_wgt_ele_1*id_wgt_ele_2*1*puweight')
        df = df.Define("weight_trigger_1_down", f'{lumi}*signs*evtweight*id_wgt_ele_1*id_wgt_ele_2*1*puweight')
        df = df.Define("weight_trigger_2_up", f'{lumi}*signs*evtweight*id_wgt_ele_1*id_wgt_ele_2*1*puweight')
        df = df.Define("weight_trigger_2_down", f'{lumi}*signs*evtweight*id_wgt_ele_1*id_wgt_ele_2*1*puweight')
        df = df.Define("ps_mass_mask", "ps_1_mass < 3 && ps_2_mass < 3 && abs(ps_1_mass-ps_2_mass) < 1")
        df = df.Define("iso_mask", f"d1_iso < 0.5 && d2_iso < 0.5") # && d3_iso < {iso_mask} && d4_iso < {iso_mask}")
        df = df.Define("mask", "pt_vis > 30 && (m_vis > 75 && m_vis < 105) && abs(H_eta) < 2.4 && d1_pt > 10 && d2_pt > 10 && ps_mass_mask && iso_mask")
        snap = df.Snapshot("ntuple", f"/work/jhornung/Haa/limits/{year}/{scope}/signal_ee_ma{ma_value}.root")

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

        f = ROOT.TFile.Open(f"/work/jhornung/Haa/limits/{year}/{scope}/signal_{scope}_ma{ma_value}.root")
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
        f = ROOT.TFile.Open(f"/work/jhornung/Haa/limits/{year}/{scope}/signal_{scope}_ma{ma_value}.root")
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

# Apply global ROOT style once — must be done before any canvas is created
# Only show errors from ROOT/RooFit, everything else is printed explicitly
ROOT.gErrorIgnoreLevel = ROOT.kError
ROOT.RooMsgService.instance().setGlobalKillBelow(RooFit.ERROR)
ROOT.RooMsgService.instance().setSilentMode(True)

ROOT.gROOT.SetStyle("Plain")
ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetTextFont(42)        # Helvetica regular (Plain sets this to bold 62)
ROOT.gStyle.SetLabelFont(42, "xyz")
ROOT.gStyle.SetTitleFont(42, "xyz")
ROOT.gStyle.SetLegendFont(42)

# Global dictionary to track boundary hits across all mass points
global_boundary_hits = {}


def adjust_variation_for_boundary(nominal_val, min_val, max_val, max_variation):
    """If nominal is within 20% of a boundary, scale down the variation."""
    range_val = max_val - min_val
    dist_from_lower = (nominal_val - min_val) / range_val
    dist_from_upper = (max_val - nominal_val) / range_val
    scale_factor = 1.0
    if dist_from_lower < 0.2:
        scale_factor = min(scale_factor, dist_from_lower / 0.2)
        print(f"  Scaling down variation due to proximity to lower boundary: {scale_factor:.2f}")
    if dist_from_upper < 0.2:
        scale_factor = min(scale_factor, dist_from_upper / 0.2)
        print(f"  Scaling down variation due to proximity to upper boundary: {scale_factor:.2f}")
    return max_variation * scale_factor


def scan_profile_nll(model, data, par, param_vars, nominal_vals, nll_min, n_points=40, scan_range=None):
    """
    Profile likelihood scan of par on the weighted MC: at each point par is fixed and all other
    parameters are refitted. The weighted NLL is scaled by sum(w)/sum(w^2) so that a step of 0.5
    corresponds to the 1 sigma interval from the MC statistics, independent of the weight normalisation.
    scan_range = (lo, hi) overrides the default range of nominal/10 to nominal*10.
    Returns the scan values, the scaled delta NLL and the 1 sigma interval (lo, hi).
    """
    sumw = data.sumEntries()
    sumw2 = 0.0
    for i in range(data.numEntries()):
        data.get(i)
        sumw2 += data.weight()**2
    scale = sumw / sumw2

    # log-spaced grid from nominal/10 to nominal*10 (clipped to the parameter range), dense near the nominal value
    par_nom = nominal_vals[param_vars.index(par)]
    lo = max(par.getMin(), abs(par_nom) / 10)
    hi = min(par.getMax(), abs(par_nom) * 10)
    if scan_range is not None:
        lo, hi = scan_range
    values = np.geomspace(lo, hi, n_points)
    if lo <= par_nom <= hi:
        values = np.unique(np.append(values, par_nom))

    was_constant = par.isConstant()
    nll_vals = []
    statuses = []
    for val in values:
        for var, nom in zip(param_vars, nominal_vals):
            var.setVal(nom)
        par.setVal(val)
        par.setConstant(True)
        # no Hesse: only the minimum is needed, and a failed Hesse would flag a good minimum as failed
        res = model.fitTo(data, RooFit.Save(), RooFit.Hesse(False), RooFit.SumW2Error(False), RooFit.Minimizer("Minuit2", "migrad"), RooFit.PrintLevel(-1), RooFit.PrintEvalErrors(-1))
        nll_vals.append(res.minNll())
        statuses.append(res.status())
    par.setConstant(was_constant)
    for var, nom in zip(param_vars, nominal_vals):
        var.setVal(nom)

    # drop scan points whose refit did not converge, their minNll is meaningless
    statuses = np.array(statuses)
    failed = values[statuses != 0]
    if len(failed) > 0:
        print(f"  {par.GetName()} scan: {len(failed)}/{len(values)} refits failed, dropped: {np.array2string(failed, precision=3)}, status: {statuses[statuses != 0]}")
    values = values[statuses == 0]
    nll_vals = np.array(nll_vals)[statuses == 0]

    dnll = (nll_vals - min(nll_min, nll_vals.min())) * scale
    inside = values[dnll < 0.5]
    if inside.min() <= par.getMin() * 1.001 or inside.max() >= par.getMax() * 0.999:
        print(f"  {par.GetName()} scan: 1 sigma interval reaches the parameter range limit, the parameter is not constrained there")
    if nll_vals.min() < nll_min - 1e-6:
        print(f"  {par.GetName()} scan: found a lower NLL than the nominal fit at {values[np.argmin(nll_vals)]:.4g}, the nominal fit is not at the minimum")
    return values, dnll, (inside.min(), inside.max())


def plot_nll_scan(values, dnll, interval, syst_vals, par_name, scope, year, ma_value):
    """Draws the scaled profile NLL of one parameter, with the fitted values of each systematic variation marked."""
    c = ROOT.TCanvas(f"c_scan_{par_name}", "scan", 1000, 800)
    c.SetLeftMargin(0.12)
    c.SetRightMargin(0.05)
    c.SetLogx()
    graph = ROOT.TGraph(len(values), array("d", values), array("d", dnll))
    graph.SetLineWidth(2)
    graph.SetLineColor(ROOT.kBlue)
    graph.SetMarkerStyle(20)
    graph.SetMarkerSize(0.8)
    graph.SetMarkerColor(ROOT.kBlue)
    y_max = max(2.0, min(dnll.max(), 10.0))
    graph.GetHistogram().SetMinimum(0)
    graph.GetHistogram().SetMaximum(y_max)
    # x axis exactly covers the scanned range, otherwise the log axis adds empty space
    graph.GetXaxis().SetLimits(values.min(), values.max())
    graph.GetXaxis().SetTitle(par_name)
    graph.GetXaxis().SetMoreLogLabels()
    graph.GetYaxis().SetTitle("#Delta NLL (scaled to MC stat.)")
    graph.Draw("ALP")

    keep = []
    half = ROOT.TLine(values.min(), 0.5, values.max(), 0.5)
    half.SetLineStyle(2)
    half.Draw()
    keep.append(half)
    for val in syst_vals:
        line = ROOT.TLine(val, 0, val, y_max)
        line.SetLineColor(ROOT.kRed)
        line.SetLineStyle(3)
        line.Draw()
        keep.append(line)

    label = ROOT.TLatex()
    label.SetNDC()
    label.SetTextFont(42)
    label.SetTextSize(0.04)
    if interval[0] <= values.min() and interval[1] >= values.max():
        interval_text = f"#Delta NLL < 0.5 over the full scan: not constrained"
    else:
        interval_text = f"1#sigma: [{interval[0]:.3g}, {interval[1]:.3g}]"
    label.DrawLatex(0.15, 0.85, f"{scope}, m_{{a}} = {ma_value} GeV, {interval_text}")
    if syst_vals:
        label.DrawLatex(0.15, 0.80, "#color[2]{dotted: fitted value per systematic variation}")

    for ext in ["png", "pdf"]:
        c.SaveAs(f"/web/jhornung/public_html/analysis_plots/{year}/ma{ma_value}/nll_scan_{par_name}_{year}_{scope}_ma{ma_value}.{ext}")


def plot_systematic_fits(model, mass, mc, systs, param_vars, snapshots, scope, year, ma_value):
    """
    One plot per systematic source: nominal, up and down MC with their own signal fits on top, and
    the shape ratio (fit_up / fit_nominal, fit_down / fit_nominal, each normalised to unit area) below.
    snapshots maps dataset name -> fitted values of param_vars.
    """
    def set_params(name):
        for var, val in zip(param_vars, snapshots[name]):
            var.setVal(val)

    groups = {}
    for syst in systs:
        name = syst.GetName()
        base = name[:-3] if name.endswith("_up") else name[:-5]
        groups.setdefault(base, {})["up" if name.endswith("_up") else "down"] = name

    norm_set = ROOT.RooArgSet(mass)
    x_values = np.linspace(mass.getMin(), mass.getMax(), 201)

    def shape(name):
        set_params(name)
        values = []
        for x in x_values:
            mass.setVal(x)
            values.append(model.getVal(norm_set))
        return np.array(values)

    for base, variations in groups.items():
        label = base.replace("weight_", "")
        styles = [("nominal", ROOT.kBlack, 20, "nominal")]
        if "up" in variations:
            styles.append((variations["up"], ROOT.kBlue, 21, "up"))
        if "down" in variations:
            styles.append((variations["down"], ROOT.kRed, 22, "down"))

        c = ROOT.TCanvas(f"canvas_{label}", f"canvas_{label}", 1200, 1200)
        c.Divide(1, 2, 0.1, 0)
        frame = mass.frame()
        legend = ROOT.TLegend(0.60, 0.62, 0.93, 0.90)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        legend.SetTextSize(0.038)

        for name, color, marker, tag in styles:
            set_params(name)
            mc[name].plotOn(frame, RooFit.MarkerColor(color), RooFit.LineColor(color), RooFit.MarkerStyle(marker), RooFit.Name(f"data_{tag}"))
            model.plotOn(frame, RooFit.LineColor(color), RooFit.Name(f"fit_{tag}"),
                         RooFit.Normalization(mc[name].sumEntries(), ROOT.RooAbsReal.NumEvent))
            yield_ratio = mc[name].sumEntries() / mc["nominal"].sumEntries()
            text = "Nominal" if tag == "nominal" else f"{label} {tag} (yield #times {yield_ratio:.3f})"
            legend.AddEntry(frame.findObject(f"fit_{tag}"), text, "lp")

        c.cd(1)
        c.cd(1).SetTopMargin(0.08)
        c.cd(1).SetRightMargin(0.05)
        c.cd(1).SetBottomMargin(0.001)
        c.cd(1).SetLeftMargin(0.12)
        c.cd(1).SetPad(0, 0.3, 1, 1)
        frame.GetXaxis().SetLabelSize(0.05)
        frame.GetYaxis().SetTitleSize(0.06)
        frame.GetYaxis().SetLabelSize(0.06)
        frame.GetYaxis().SetTitleOffset(1.1)
        frame.GetYaxis().SetTitle("Events / 0.5 GeV")
        frame.SetMinimum(0)
        frame.Draw()
        legend.Draw()

        header = ROOT.TLatex()
        header.SetNDC()
        header.SetTextFont(42)
        header.SetTextSize(0.055)
        header.DrawLatex(0.12, 0.93, "#it{Private Work} #bf{(CMS Simulation)}")
        ann = ROOT.TLatex()
        ann.SetNDC()
        ann.SetTextFont(42)
        ann.SetTextSize(0.048)
        scope_str = "ee SR" if scope == "ee" else "#mu#mu SR"
        for y_pos, text in [(0.86, f"Systematic: {label}"), (0.80, scope_str), (0.74, f"m_{{a}} = {ma_value} GeV")]:
            ann.DrawLatex(0.15, y_pos, text)

        # shape ratio of the fitted curves
        nominal_shape = shape("nominal")
        graphs = []
        max_dev = 0.01
        for name, color, marker, tag in styles[1:]:
            ratio = shape(name) / nominal_shape
            max_dev = max(max_dev, np.max(np.abs(ratio - 1)))
            g = ROOT.TGraph(len(x_values), x_values.astype("d"), ratio.astype("d"))
            g.SetLineColor(color)
            g.SetLineWidth(2)
            graphs.append(g)

        c.cd(2)
        c.cd(2).SetTopMargin(0.001)
        c.cd(2).SetRightMargin(0.05)
        c.cd(2).SetBottomMargin(0.3)
        c.cd(2).SetLeftMargin(0.12)
        c.cd(2).SetPad(0, 0, 1, 0.3)
        frame2 = mass.frame()
        frame2.SetMinimum(1 - 1.3 * max_dev)
        frame2.SetMaximum(1 + 1.3 * max_dev)
        frame2.GetXaxis().SetTitle("m_{4K} in GeV")
        frame2.GetYaxis().SetTitle("#frac{Fit var.}{Fit nom.} (shape)")
        frame2.GetXaxis().SetTitleSize(0.14)
        frame2.GetXaxis().SetTitleOffset(1)
        frame2.GetXaxis().SetLabelSize(0.11)
        frame2.GetYaxis().SetTitleSize(0.09)
        frame2.GetYaxis().SetTitleOffset(0.45)
        frame2.GetYaxis().SetLabelSize(0.11)
        frame2.Draw()
        one = ROOT.TLine(mass.getMin(), 1, mass.getMax(), 1)
        one.SetLineStyle(2)
        one.Draw()
        for g in graphs:
            g.Draw("L SAME")

        out = f"/web/jhornung/public_html/analysis_plots/{year}/ma{ma_value}/sig_syst_{label}_{year}_{scope}_ma{ma_value}"
        c.SaveAs(f"{out}.png")
        c.SaveAs(f"{out}.pdf")

    set_params("nominal")
    mass.setVal(125)


def load_dcb_tails(scope, year, ma_value):
    """Tail parameters of the double-sided CB for this mass point, from model D of compare_signal_models.py."""
    path = f"/work/jhornung/Haa/limits/{year}/{scope}/dcb_tails_{scope}.json"
    if not os.path.exists(path):
        sys.exit(f"{path} not found, run compare_signal_models.py {scope} {year} first")
    with open(path) as f:
        tails = json.load(f)["tails"]
    if ma_value not in tails:
        sys.exit(f"no DCB tails for ma = {ma_value} in {path}")
    return tails[ma_value]


def process_signal(scope, year, ma_value, signal_path, global_boundary_hits):
    os.makedirs(f"/work/jhornung/Haa/limits/{year}/{scope}/ma{ma_value}", exist_ok=True)
    os.makedirs(f"/web/jhornung/public_html/analysis_plots/{year}/ma{ma_value}", exist_ok=True)

    mc, yields, yield_vars, systs, mass = prep_data(scope, year, ma_value, signal_path)

    # Double-sided Crystal Ball. The tails are fixed to the log-linear-in-ma description fitted to all mass points
    # together (model D in compare_signal_models.py); only MH and sigma are fitted here and varied by the systematics.
    # all workspace objects carry the year: text2workspace merges same-named variables of different
    # workspaces (RecycleConflictNodes), so without it one year's pdf would silently use another year's values.
    # The nuisances ({scope}_mean_nuisance, {scope}_sigma_nuisance) keep their names on purpose.
    tail_values = load_dcb_tails(scope, year, ma_value)
    MH = ROOT.RooRealVar(f"{scope}_{year}_MH", "MH", 125.0, 120.0, 130.0)
    sigma = ROOT.RooRealVar(f"{scope}_{year}_sigma", "width of dcb core", 1.0, 0.1, 10.0)
    tails = {}
    for name in ["alphaL", "nL", "alphaR", "nR"]:
        tails[name] = ROOT.RooRealVar(f"{scope}_{year}_{name}", name, tail_values[name])
        tails[name].setConstant(True)

    model = ROOT.RooCrystalBall(f"{scope}_{year}_model", "double-sided CB", mass, MH, sigma, tails["alphaL"], tails["nL"], tails["alphaR"], tails["nR"])

    mean_vals = {}
    sigma_vals = {}

    model.fitTo(mc["nominal"], RooFit.Save(), RooFit.SumW2Error(False), RooFit.Minimizer("Minuit2", "migrad"), RooFit.PrintLevel(-1), RooFit.PrintEvalErrors(-1))
    # AsymptoticError so that the parameter errors reflect the MC statistics and not the weight normalisation
    nominal_result = model.fitTo(mc["nominal"], RooFit.Save(), RooFit.AsymptoticError(True), RooFit.Minimizer("Minuit2", "migrad"), RooFit.PrintLevel(-1), RooFit.PrintEvalErrors(-1))

    param_vars = [MH, sigma]
    snapshots = {"nominal": [v.getVal() for v in param_vars]}

    mean_vals["nominal"] = [MH.getVal(), MH.getError()]
    sigma_vals["nominal"] = [sigma.getVal(), sigma.getError()]

    # Track all boundary hits for summary at the end
    all_boundary_hits = []
    if ma_value not in global_boundary_hits:
        global_boundary_hits[ma_value] = []

    def check_boundaries(label):
        for var in param_vars:
            margin = 0.01 * (var.getMax() - var.getMin())
            if var.getVal() - var.getMin() < margin or var.getMax() - var.getVal() < margin:
                msg = f"[{label}] {var.GetName()} at its range limit: {var.getVal():.3f} [{var.getMin()}, {var.getMax()}]"
                print(f"  WARNING: {msg}")
                all_boundary_hits.append(msg)
                global_boundary_hits[ma_value].append(msg)

    check_boundaries("NOMINAL")

    for syst in systs:
        # start every variation from the nominal fit, otherwise the fits drift into different minima
        for var, val in zip(param_vars, snapshots["nominal"]):
            var.setVal(val)
        model.fitTo(mc[syst.GetName()], RooFit.Save(), RooFit.SumW2Error(False), RooFit.Minimizer("Minuit2", "migrad"), RooFit.PrintLevel(-1), RooFit.PrintEvalErrors(-1))
        model.fitTo(mc[syst.GetName()], RooFit.Save(), RooFit.SumW2Error(False), RooFit.Minimizer("Minuit2", "migrad"), RooFit.PrintLevel(-1), RooFit.PrintEvalErrors(-1))

        snapshots[syst.GetName()] = [v.getVal() for v in param_vars]
        mean_vals[syst.GetName()] = [MH.getVal(), MH.getError()]
        sigma_vals[syst.GetName()] = [sigma.getVal(), sigma.getError()]
        check_boundaries(syst.GetName().upper())

    for var, val in zip(param_vars, snapshots["nominal"]):
        var.setVal(val)

    # Compare each shift to the MC-stat uncertainty of the nominal fit; shifts well below 1 sigma are fit noise
    print(f"\nShift / fit error for {scope} ma={ma_value} (nominal fit status {nominal_result.status()}, covQual {nominal_result.covQual()}):")
    for name, vals in [("MH", mean_vals), ("sigma", sigma_vals)]:
        nom, err = vals["nominal"]
        for syst in systs:
            shift = vals[syst.GetName()][0] - nom
            # err is 0 when the error computation failed (e.g. parameter at its limit)
            ratio = f"{shift/err:+.2f}" if err > 0 else "n/a (no error)"
            print(f"  {name:6s} {syst.GetName():22s} shift={shift:+.4g}  err={err:.4g}  shift/err={ratio}")

    max_mean_var = np.max([np.abs(mean_vals[s.GetName()][0] - mean_vals["nominal"][0]) / mean_vals["nominal"][0] for s in systs])
    max_sigma_var = np.max([np.abs(sigma_vals[s.GetName()][0] - sigma_vals["nominal"][0]) / sigma_vals["nominal"][0] for s in systs])

    mean_nuisance = ROOT.RooRealVar(f"{scope}_mean_nuisance", "mean nuisance", 0, -1, 1)
    sigma_nuisance = ROOT.RooRealVar(f"{scope}_sigma_nuisance", "sigma nuisance", 0, -1, 1)
    mean_nuisance.setConstant(True)
    sigma_nuisance.setConstant(True)

    mean_formula = ROOT.RooFormulaVar(f"{scope}_{year}_mean_formula", f"mean formula", f"@0*(1+{max_mean_var}*@1)", ROOT.RooArgList(MH, mean_nuisance))
    sigma_formula = ROOT.RooFormulaVar(f"{scope}_{year}_sigma_formula", f"sigma formula", f"@0*(1+{max_sigma_var}*@1)", ROOT.RooArgList(sigma, sigma_nuisance))
    # no nuisances on the tails: they are fixed from the fit to all mass points, the systematics only move MH and sigma

    model_nom = ROOT.RooCrystalBall(f"model_Haa_H_mass_{year}_{scope}_ma{ma_value}", "double-sided CB", mass, mean_formula, sigma_formula, tails["alphaL"], tails["nL"], tails["alphaR"], tails["nR"])
    result = model_nom.fitTo(mc["nominal"], RooFit.Save(), RooFit.SumW2Error(False), RooFit.Minimizer("Minuit2", "migrad"), RooFit.PrintLevel(-1), RooFit.PrintEvalErrors(-1))

    eff_nom = float("%.2f"%(mc["nominal"].sumEntries()))
    eff_nom = ROOT.RooRealVar(f"model_Haa_H_mass_{year}_{scope}_ma{ma_value}_norm", "eff_nom", eff_nom)

    MH.setConstant(True)
    sigma.setConstant(True)

    f_out = ROOT.TFile(f"/work/jhornung/Haa/limits/{year}/{scope}/ma{ma_value}/workspace_sig_{scope}_ma{ma_value}.root", "RECREATE")
    w_sig_nom = ROOT.RooWorkspace(f"workspace_sig_nom", f"workspace_sig")
    getattr(w_sig_nom, 'import')(model_nom)
    getattr(w_sig_nom, 'import')(eff_nom)
    w_sig_nom.Print()
    w_sig_nom.Write()
    f_out.Close()

    # --- Plotting ---
    c = ROOT.TCanvas("canvas", "canvas", 1200, 1200)
    c.Divide(1, 2, 0.1, 0)

    legend = ROOT.TLegend(0.60, 0.74, 0.93, 0.90)

    xframe = mass.frame()

    mc["nominal"].plotOn(xframe)
    hist = xframe.getObject(int(xframe.numItems()-1))

    model_nom.plotOn(xframe, RooFit.LineColor(ROOT.kRed), RooFit.LineStyle(ROOT.kSolid))
    sig_fit = xframe.getObject(int(xframe.numItems()-1))

    mc["nominal"].plotOn(xframe)

    legend.AddEntry(hist, "Simulated Data", "ep")
    legend.AddEntry(sig_fit, "Signal Fit (DCB)", "l")

    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.042)

    num_params = result.floatParsFinal().getSize()
    chi2 = xframe.chiSquare(num_params)

    hpull = xframe.pullHist()
    xframe2 = mass.frame()
    xframe2.SetTitle("")
    xframe2.addPlotable(hpull, "P")

    c.cd(1)
    c.cd(1).SetTopMargin(0.08)
    c.cd(1).SetRightMargin(0.05)
    c.cd(1).SetBottomMargin(0.001)
    c.cd(1).SetLeftMargin(0.12)
    ROOT.gPad.SetBottomMargin(0.001)
    c.cd(1).SetPad(0, 0.3, 1, 1)

    xframe.GetXaxis().SetLabelSize(0.05)
    xframe.GetYaxis().SetTitleSize(0.06)
    xframe.GetYaxis().SetLabelSize(0.06)
    xframe.GetYaxis().SetTitleOffset(1.1)
    xframe.GetYaxis().SetTitle("Events / 0.5 GeV")
    xframe.SetMinimum(0)

    xframe.Draw()
    legend.Draw()

    pw_label = ROOT.TLatex()
    pw_label.SetNDC()
    pw_label.SetTextFont(42)
    pw_label.SetTextSize(0.055)
    pw_label.DrawLatex(0.12, 0.93, "#it{Private Work} #bf{(CMS Simulation)}")

    lumi_map = {
        "2018": "59.8 fb^{-1} (13 TeV)",
        "2017": "41.5 fb^{-1} (13 TeV)",
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
        (0.80, "Signal MC"),
        (0.74, scope_str),
        (0.68, "#chi^{2}/dof = " + f"{chi2:.2f}"),
        (0.62, f"m_{{a}} = {ma_value} GeV"),
    ]:
        ann_label.DrawLatex(0.15, y_pos, text)

    c.cd(2)
    c.cd(2).SetTopMargin(0.001)
    c.cd(2).SetRightMargin(0.05)
    c.cd(2).SetBottomMargin(0.3)
    c.cd(2).SetLeftMargin(0.12)
    c.cd(2).SetPad(0, 0, 1, 0.3)

    xframe2.GetXaxis().SetTitle("m_{4K} in GeV")
    xframe2.GetYaxis().SetTitle("#frac{MC - Fit}{Error}")
    xframe2.GetXaxis().SetTitleSize(0.14)
    xframe2.GetXaxis().SetTitleOffset(1)
    xframe2.GetXaxis().SetLabelSize(0.11)
    xframe2.GetYaxis().SetTitleSize(0.11)
    xframe2.GetYaxis().SetTitleOffset(0.35)
    xframe2.GetYaxis().SetLabelSize(0.11)

    xframe2.Draw()

    zero = ROOT.TLine(115, 0, 135, 0)
    zero.SetLineColor(ROOT.kRed)
    zero.SetLineWidth(2)
    zero.SetLineStyle(2)
    zero.Draw()

    model_nom.Print("t")

    c.SaveAs(f"/web/jhornung/public_html/analysis_plots/{year}/ma{ma_value}/sig_fit_{year}_{scope}_ma{ma_value}_new.png")
    c.SaveAs(f"/web/jhornung/public_html/analysis_plots/{year}/ma{ma_value}/sig_fit_{year}_{scope}_ma{ma_value}_new.pdf")

    plot_systematic_fits(model, mass, mc, systs, param_vars, snapshots, scope, year, ma_value)

    print(f"\n=== Results for ma = {ma_value} GeV ===")
    print("DCB tails (fixed): " + ", ".join(f"{name} = {var.getVal():.3f}" for name, var in tails.items()))
    for y in yields.keys():
        print(f"Yield for {y}: {yields[y]/yields['nominal']}")

    print(f"Max mean variation: {max_mean_var}")
    print(f"Max sigma variation: {max_sigma_var}")

    params = result.floatParsFinal()
    param_names = [params[i].GetName() for i in range(params.getSize())]
    print("Parameter order in covariance matrix:")
    print(param_names)

    cov = result.covarianceMatrix()
    cov.Print()

    print(eff_nom)

    # Print boundary hit summary
    if all_boundary_hits:
        print("\n" + "=" * 60)
        print(f"BOUNDARY HITS SUMMARY for {scope} ma={ma_value}")
        print("=" * 60)
        for hit in all_boundary_hits:
            print(f"  {hit}")
        print("=" * 60)
    else:
        print(f"\nNo boundary hits for {scope} ma={ma_value}")

    print("=" * 50 + "\n")

    # Dump results to JSON
    output_dict = {
        "ma_value": ma_value,
        "scope": scope,
        "year": year,
        "yield_variations": {},
        "max_parameter_variations": {
            "mean": float(max_mean_var),
            "sigma": float(max_sigma_var)
        }
    }
    for y in yields.keys():
        if y != "nominal":
            output_dict["yield_variations"][y] = float(yields[y]/yields['nominal'])

    json_output_path = f"/work/jhornung/Haa/limits/{year}/{scope}/ma{ma_value}/systematics_results_{scope}_ma{ma_value}.json"
    with open(json_output_path, 'w') as json_file:
        json.dump(output_dict, json_file, indent=4)

    print(f"Results saved to: {json_output_path}")


if __name__ == "__main__":
    scope = sys.argv[1]
    year = sys.argv[2]

    signal_dir = glob.glob(f"/ceph/jhornung/MC_2018/{year}/Haa*")

    for signal_path in signal_dir:
        signal = os.path.basename(signal_path)

        # Extract ma value
        if "ma" in signal:
            ma_value = signal.split("ma")[1].split("_")[0]
            print(f"Processing signal: {signal}, ma = {ma_value}")
            process_signal(scope, year, ma_value, signal_path, global_boundary_hits)

    # After processing all mass points, print global boundary summary
    print("\n" + "="*80)
    print(f"GLOBAL BOUNDARY HITS SUMMARY FOR ALL MASS POINTS ({scope}, {year})")
    print("="*80)

    if global_boundary_hits:
        total_hits = sum(len(hits) for hits in global_boundary_hits.values())
        print(f"Total mass points with boundary hits: {len(global_boundary_hits)}")
        print(f"Total boundary hits: {total_hits}")
        print("")
    
        for ma_val in sorted(global_boundary_hits.keys(), key=lambda x: float(x)):
            if global_boundary_hits[ma_val]:
                print(f"\nma = {ma_val} GeV ({len(global_boundary_hits[ma_val])} hits):")
                for hit in global_boundary_hits[ma_val]:
                    print(f"  {hit}")
    else:
        print("No boundary hits detected for any mass point!")

    print("="*80 + "\n")