import ROOT

f = ROOT.TFile("higgsCombine.bestfit.AsymptoticLimits.mH125.root")
w = f.Get("w")
w.Print("v")

n_bins = 65
binning = ROOT.RooFit.Binning(n_bins,70,200)

can = ROOT.TCanvas()
plot = w.var("H_mass").frame()
w.data("data_obs").plotOn( plot, binning )

# Load the S+B model
sb_model = w.pdf("model_s").getPdf("H_mass")
b_model = w.pdf("model_b").getPdf("H_mass")


# Prefit
sb_model.plotOn( plot, ROOT.RooFit.LineColor(2), ROOT.RooFit.Name("prefit") )
w.data("data_obs").plotOn( plot, binning )
# Postfit
w.loadSnapshot("MultiDimFit")
sb_model.plotOn( plot, ROOT.RooFit.LineColor(4), ROOT.RooFit.Name("postfit") )
r_bestfit = w.var("r").getVal()

#b_model.plotOn( plot, ROOT.RooFit.LineColor(1), ROOT.RooFit.Name("postfit") )

plot.Draw()

leg = ROOT.TLegend(0.55,0.6,0.85,0.85)
leg.AddEntry("prefit", "Prefit S+B model (r=1.00)", "L")
leg.AddEntry("postfit", "Postfit S+B model (r=%.2f)"%r_bestfit, "L")
leg.Draw("Same")

can.Update()
can.SaveAs("/web/jhornung/public_html/test_datacard.png")
