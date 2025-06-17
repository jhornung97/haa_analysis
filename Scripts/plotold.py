import ROOT
import math
import numpy as np
import sys
from array import array
import os

ROOT.gROOT.SetBatch(True)


histG = ROOT.TH1F("histG", " ; ;#events (normalized)", 70, 0, 140)
geant = ROOT.TFile.Open("/ceph/bmaier/Haa/2018/ttHToNonbb_M125_TuneCP5_13TeV-powheg-pythia8+RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2+MINIAODSIM/01368E94-683A-0844-ACEE-0C7770C6F46D.root")
tree = geant.Get("Events")

pt = array("f",[0,0,0,0,0,0,0])
eta = array("f",[0,0,0,0,0,0,0])
phi = array("f",[0,0,0,0,0,0,0])
mass = array("f",[0,0,0,0,0,0,0])
test = array("i",[0])

tree.SetBranchAddress("Muon_pt", pt)
tree.SetBranchAddress("Muon_eta", eta)
tree.SetBranchAddress("Muon_phi", phi)
tree.SetBranchAddress("Muon_mass", mass)
tree.SetBranchAddress("nMuon", test)

nevents = tree.GetEntries()

for i in range(nevents):
  invMass = ROOT.TLorentzVector()
  tree.GetEntry(i)
  for j in range(test[0]):
    if pt[j] <= 10 or abs(eta[j]) > 2.4:
      test[0] = 0
  if test[0] > 1:
    for j in range(test[0]):
      vec = ROOT.TLorentzVector()
      vec.SetPtEtaPhiM(pt[j],eta[j],phi[j],mass[j])
      invMass += vec
    histG.Fill(invMass.Mag())

histG.Scale(1./histG.Integral())

hist.Draw("HISTE")
	
c.SaveAs("/work/jhornung/Haa/Plots/testold.pdf","pdf")
