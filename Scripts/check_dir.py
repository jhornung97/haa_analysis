import uproot
import sys
import glob
import os

dir_to_be_checked = sys.argv[1]

files = glob.glob(os.path.join(dir_to_be_checked, "Single*"))

event_sum = 0

for file in files:
    ifile = uproot.open(file)
    tree = ifile["ntuple"]
    nevents = tree.num_entries
    print(f"Number of events in {file}: {nevents/1e6}M")
    event_sum += nevents

print(f"Total number of events in {dir_to_be_checked}: {event_sum}")

hadded_file = os.path.join(dir_to_be_checked, "single_muon_data.root")
hadded_tree = uproot.open(hadded_file)["ntuple"]
nevets = hadded_tree.num_entries

print(f"Total number of events in {hadded_file}: {nevets}")