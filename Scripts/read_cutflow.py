import sys
import ROOT

ifile = ROOT.TFile(sys.argv[1])
cutflow = ifile.Get("cutflow")
ntuple = ifile.Get("ntuple")

if "old" in sys.argv[1]:
    ec_before = cutflow.GetBinContent(1)
    ec_after = ntuple.GetEntries()
    print("event count before workflow: {}".format(ec_before))
    print("event count after first workflow: {}".format(ec_after))

else:
    ec_after = ntuple.GetEntries()
    print("event count after second workflow: {}".format(ec_after))
