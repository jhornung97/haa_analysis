import pandas as pd
import json
import ROOT
import sys
import subprocess

'''
data = pd.read_csv('/work/jhornung/Haa/crosssections.txt', sep=",") 
print(data)
'''

with open('/work/jhornung/Haa/samples/crosssections.json') as f:
    data = json.load(f)

ifile = ROOT.TFile.Open(sys.argv[1],"READ")

cutflow = ifile.Get("cutflow")

if "Haa" in sys.argv[1]:
    print("signal file")
    N_events = ifile.Get("cutflow").GetBinContent(1)
    dummy = 1.0/N_events
    weight= "%.8f"%(dummy)
else:
    #Extract the dataset name from the input file path
    dataset_name = sys.argv[1]

    #Find the cross section for the dataset
    crosssection = None
    for process, details in data.items():
        if dataset_name in details['datasets']:
            crosssection = details['xsec']
            break

    if crosssection is None:
        print(f'Error: no cross section found for {dataset_name}')
        sys.exit(1)

    #crosssection = data[sys.argv[1].split("/")[-1]]
    N_events = ifile.Get("cutflow").GetBinContent(1)
    normalizedWeight = crosssection/N_events
    print(ifile.Get("cutflow").GetBinContent(1))
    print("cross section: " + str(crosssection))
    weight="%.8f"%(normalizedWeight)


print("weight: " + str(weight))

branches = []

for branch in ifile.ntuple.GetListOfBranches():
    name = branch.GetName()
    branches.append(str(name))

ifile.Close()

ofile=sys.argv[-1].replace(".root","_final.root")
df = ROOT.RDataFrame("ntuple", sys.argv[1])
df = df.Define("evtweight",weight)
branches.append("evtweight")

df.Snapshot("ntuple", ofile, branches)

#ofile = ROOT.TFile.Open(ofile,"UPDATE")
#.cutflow.Write()

subprocess.call(f"mv {ofile} {sys.argv[-1]}",shell=True)
