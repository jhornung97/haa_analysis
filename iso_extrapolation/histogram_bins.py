import numpy as np
import csv

#years = ["2016preVFP", "2016postVFP", "2017", "2018"]
years = ["2024"]

#scopes = ["mm", "ee"]
scopes = ["mm"]

cuts_2017 = np.arange(.5, 10.5, .5)
cuts_2016_2018_mm = np.arange(1, 10.5, .5)
cuts_2018_ee = np.arange(1.5, 10.5, .5)
cuts_2024 = np.arange(1.5, 10.5, .5)


with open('combinations.txt', 'w') as csvfile:
    for year in years:
        for scope in scopes:
            if year == "2017":
                cuts = cuts_2017
            elif year == "2016preVFP" or year == "2016postVFP":
                cuts = cuts_2016_2018_mm
            else:  
                if scope == "mm":
                    cuts = cuts_2016_2018_mm
                else:
                    cuts = cuts_2018_ee
            
            for cut in cuts:
                csvfile.write(f"{scope} {year} {cut}\n")