import numpy as np
import ROOT
import sys
import subprocess

ifile = np.genfromtxt('/work/jhornung/Haa/crosssections.txt', delimiter=",", dtype=str)

samplenames = ifile[0]
crosssections = ifile[1]

compared_pairs = set()

for i, sample1 in enumerate(samplenames):
    if ("TT" in sample1 or "ST" in sample1) and "20UL18" in sample1:
        for j, sample2 in enumerate(samplenames):
            if "20UL17" in sample2 and sample1[:15] == sample2[:15]:
                pair = tuple(sorted((i, j)))
                if pair not in compared_pairs:
                    compared_pairs.add(pair)
                    print(f"Comparing {sample1[:15]} from 18 and {sample2[:15]} from 17")
                    print(f"Cross sections are different: {crosssections[i]} and {crosssections[j]}")

for i, sample1 in enumerate(samplenames):
    print(f"{sample1}: {crosssections[i]}")

