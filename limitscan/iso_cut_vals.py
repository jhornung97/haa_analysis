import numpy as np
import csv

numbers = np.linspace(0.1, 20., 200)  # Generate 200 numbers from 0.1 to 20.0
numbers = np.round(numbers, 1).tolist()  # Round to 2 decimal places and convert to list
with open('iso_cut_vals.txt', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile, delimiter='\n')
    writer.writerow(numbers)