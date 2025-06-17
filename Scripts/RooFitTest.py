import ROOT

# Create a histogram with a Gaussian distribution
hist = ROOT.TH1F("hist", "Example Histogram", 100, -5, 5)
hist.FillRandom("gaus", 1000)  # Fill with 1000 events

# Define the fit model
x = ROOT.RooRealVar("x", "x", -5, 5)
mean = ROOT.RooRealVar("mean", "mean", 0, -1, 1)
sigma = ROOT.RooRealVar("sigma", "sigma", 1, 0.1, 2)
gauss = ROOT.RooGaussian("gauss", "gaussian PDF", x, mean, sigma)

# Convert histogram to RooDataHist
data = ROOT.RooDataHist("data", "data", ROOT.RooArgList(x), hist)

# Perform the fit without scaling factor
fit_result = gauss.fitTo(data, ROOT.RooFit.Save(), ROOT.RooFit.Minimizer("Minuit2", "migrad"))

# Print fit results
fit_result.Print()

# Apply a scaling factor to the histogram
scaling_factor = 1000
hist.Scale(scaling_factor)

# Convert scaled histogram to RooDataHist
data_scaled = ROOT.RooDataHist("data_scaled", "data_scaled", ROOT.RooArgList(x), hist)

# Perform the fit with scaling factor
fit_result_scaled = gauss.fitTo(data_scaled, ROOT.RooFit.Save(), ROOT.RooFit.Minimizer("Minuit2", "migrad"))

# Print fit results
fit_result_scaled.Print()

# Plot the results
c = ROOT.TCanvas("c", "c", 800, 600)
frame = x.frame()
data.plotOn(frame)
gauss.plotOn(frame)
frame.Draw()
c.SaveAs("/web/jhornung/public_html/fit_result.png")

c_scaled = ROOT.TCanvas("c_scaled", "c_scaled", 800, 600)
frame_scaled = x.frame()
data_scaled.plotOn(frame_scaled)
gauss.plotOn(frame_scaled)
frame_scaled.Draw()
c_scaled.SaveAs("/web/jhornung/public_html/fit_result_scaled.png")