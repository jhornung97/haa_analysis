import json
import re

ifile = '/work/jhornung/Haa/crosssections.txt'
ofile = '/work/jhornung/Haa/crosssections.json'

with open(ifile, 'r') as f:
    lines = f.readlines()

names = lines[0].strip().split(',')
crosssections = lines[1].strip().split(',')

datasets = {}
for name, xsec in zip(names, crosssections):
    match = re.match(r'([A-Za-z0-9_-]+)_Tune', name)
    if match:
        process_name = match.group(1)
        if process_name not in datasets:
            datasets[process_name] = {'datasets': [], 'xsec': xsec}
        datasets[process_name]['datasets'].append(name)
        datasets[process_name]['xsec'] = xsec
    else:
        print(f'No match for {name}')

#Check if all xsecs for a given process are the same
for process, data in datasets.items():
    xsecs = [data['xsec']]
    if len(set(xsecs)) > 1:
        print(f'Warning: {process} has different xsecs: {xsecs}')
        #Choose the xsec with fewer decimal places
        chosen_xsec = min(xsecs, key=lambda x: len(x.split('.')[1]) if '.' in x else 0)
        print(f'Choosing {chosen_xsec}')
        #Update the datasets with the chosen xsec
        data['xsec'] = chosen_xsec
    else:
        data['xsec'] = xsecs[0]

#Convert the xsecs to floats
for process, data in datasets.items():
    data['xsec'] = float(data['xsec'])

with open(ofile, 'w') as f:
    json.dump(datasets, f, indent=4)

print(f'Wrote {ofile}')