import ROOT
import numpy as np

double_muon_files = np.loadtxt("/work/jhornung/Haa/filelist.txt", dtype=str)

single_muon_df = ROOT.RDataFrame("ntuple", "/ceph/jhornung/analysis/2018/single_muon_data/mm/single_muon_data.root")
single_muon_cols = single_muon_df.AsNumpy(["run", "lumi", "event"])

single_muon_runs = np.unique(single_muon_cols["run"])
single_muon_lumi = np.unique(single_muon_cols["lumi"])
single_muon_events = np.unique(single_muon_cols["event"])

#print(single_muon_runs)

for f in double_muon_files:
    double_muon_df = ROOT.RDataFrame("ntuple", f)
    print("{} has been opened".format(f))
    
    if double_muon_df.Count().GetValue() == 0:
        print("File is empty, loop will continue.")
    else:    
        print("Tree is full.")
        double_muon_cols = double_muon_df.AsNumpy(["run", "lumi", "event"])
        
        double_muon_runs = np.unique(double_muon_cols["run"])
        double_muon_lumi = np.unique(double_muon_cols["lumi"])
        double_muon_events = np.unique(double_muon_cols["event"])

        runs_overlap = np.all(np.isin(double_muon_runs, single_muon_runs))
        lumi_overlap = np.all(np.isin(double_muon_lumi, single_muon_lumi))
        events_overlap = np.all(np.isin(double_muon_events, single_muon_events))

        if (runs_overlap == True): 
            print("File {} only has matching runs.".format(f))
            if (lumi_overlap == True):
                print("File {} only has matching lumi sections.".format(f))
                if (events_overlap == True):
                    print("File {} only has matching events.".format(f))
   
'''
file2 = ROOT.TFile("/ceph/jhornung/analysis/2018/double_muon_data/mm/double_muon_data.root")

# Lese die Bäume aus den Dateien
tree1 = file1.Get("ntuple")  # "tree_name" durch den tatsächlichen Namen des Baumes ersetzen
tree2 = file2.Get("ntuple")  # "tree_name" durch den tatsächlichen Namen des Baumes ersetzen

for i, event_one in enumerate(tree1):
    has_match = False
    for j, event_two in enumerate(tree2):
        print([i, j])
        if (event_one.lumi == event_two.lumi) and (event_one.run == event_two.run) and (event_one.event == event_two.event):
            print("Event {} matches with event {}.".format(i, j))
            has_match = True
            break
    if has_match == False:
        print("event {} would be appended to a cloned double muon file.".format(i))

# Klonen Sie den Baum aus Datei 1
output_tree = tree1.Clone("ntuple")  # Das Argument 0 bedeutet, dass nur die Struktur des Baumes kopiert wird
# Erstelle eine Ausgabedatei im Schreibmodus
output_file = ROOT.TFile("/ceph/jhornung/analysis/2018/full_data/mm/full_mm_data.root", "RECREATE")

# Schreibe den Ausgabebaum in die Ausgabedatei
output_file.cd()
output_tree.Write()
output_file.Close()

# Schließe die Dateien
file1.Close()
file2.Close()
'''
