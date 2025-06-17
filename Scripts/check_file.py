import sys
import ROOT

ifile = ROOT.TFile(sys.argv[1])
tree = ifile.Get("ntuple")

print(tree.GetEntries())