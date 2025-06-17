import ROOT
import math
import numpy as np
import sys
from array import array 
import os

ROOT.gROOT.SetBatch(True)

ROOT.gROOT.SetStyle("Plain")
ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
    
c = ROOT.TCanvas("canvas", "canvas", 1920, 1080)
c.Divide(1,2,0.1,0)

hist = ROOT.TH1F("hist", " ; ;#events", 23, 0, 230) 
nano = ROOT.TFile.Open("/work/jhornung/Haa/samples/signal/HaaKKKK_1227059.root")
events = nano.Get("Events")

for n, event in enumerate(events):																							#event loop
	
	nmu = event.nPFCands
	pts = np.array(event.PFCands_pt)
	etas = np.array(event.PFCands_eta)
	phis = np.array(event.PFCands_phi)
	masses = np.array(event.PFCands_mass)
	pdgid = np.array(event.PFCands_pdgId)
	
	hadrpt = []
	hadreta = []
	hadrphi = []
	hadrmass = []
	nhadr = 0
	
	for i in range(nmu):
		if np.abs(pdgid[i]) == 211:
			hadrpt.append(pts[i])
			hadreta.append(etas[i])
			hadrphi.append(phis[i])
			hadrmass.append(masses[i])
			nhadr += 1
			
	hadrs = zip(hadrpt, hadreta, hadrphi, hadrmass)
	
	sortedHadrs = sorted(hadrs, reverse=True)
	
	sortedhadrpts, sortedhadretas, sortedhadrphis, sortedhadrmasses = zip(*sortedHadrs)
	
	#print(sortedhadrpts[:4], sortedhadretas[:4])
	
	if nhadr >= 4:
			
		invSyst = ROOT.TLorentzVector()	
		
		for i in range(4):		#calc some stuff
			
			hadr = ROOT.TLorentzVector() 
			hadr.SetPtEtaPhiM(sortedhadrpts[i],sortedhadretas[i],sortedhadrphis[i],sortedhadrmasses[i])
			invSyst += hadr
		
		print(invSyst.Mag())

		hist.Fill(invSyst.Mag())

#hist.Scale(1./hist.Integral())
hist.Draw("HIST")
	
c.SaveAs("/work/jhornung/Haa/Plots/testHaa4HardestHadrons.pdf","pdf")
