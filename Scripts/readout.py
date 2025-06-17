import numpy as np

ifile = "/work/jhornung/Haa/crosssections.txt"

open(ifile, 'r')
names = open(ifile, 'r').readline().strip().split(',')
css = np.genfromtxt("/work/jhornung/Haa/crosssections.txt", delimiter=',', skip_header=1)

print("sample name css")
for i, name in enumerate(names):
    print("{} {}".format(name, css[i]))
