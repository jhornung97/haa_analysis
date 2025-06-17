import json
import re
import sys
import os
import glob

#Paths to json file and sample directory

json_file = '/work/jhornung/Haa/crosssections.json'
sample_dir = sys.argv[1]

with open(json_file) as f:
    data = json.load(f)

sample_list = glob.glob(os.path.join(sample_dir, '*'))

#Update json file with new samples
for sample in sample_list:
    sample_name = sample.split('/')[-1] + '.root'
    match = re.match(r'([A-Za-z0-9_-]+)_Tune', sample_name)
    if match:
        process_name = match.group(1)
        if process_name in data:
            if sample_name not in data[process_name]['datasets']:
                data[process_name]['datasets'].append(sample_name)
        else:
            print(f'Process {process_name} not found in json file. Skipping...')
    else:
        print(f'Could not extract process name from {sample_name}. Skipping...')

#Write updated json file
with open(json_file, 'w') as f:
    json.dump(data, f, indent=4)

print('Updated json file with new samples')
