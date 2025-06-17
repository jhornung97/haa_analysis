import numpy as np
from itertools import combinations

def CheckPhiDiff(delPhi):

    delPhi = np.where(delPhi < -np.pi, delPhi + 2*np.pi, delPhi)
    delPhi = np.where(delPhi > np.pi, delPhi - 2*np.pi, delPhi)

    return delPhi

def PairDiff(it, kin1, kin2):

    Pairs = np.array(np.meshgrid(kin1[it,:],kin2[it,:])).T.reshape(-1,2)
#    Pairs = np.array(list(combinations(kin1[it,:], 2)))
#    print(Pairs)

    Deltas = []

    for pair in Pairs:

        delta = pair[0] - pair[1]
        Deltas.append(delta)

    return Deltas

def PairDiffAk(it, kin1, kin2):

    Pairs = np.array(np.meshgrid(kin1[it,:],kin2)).T.reshape(-1,2)

    #print(Pairs)                                                                                        

    Deltas = []

    for pair in Pairs:

        delta = pair[0] - pair[1]
        Deltas.append(delta)

    return Deltas

def weights(weight):

        if weight > 0:
                weight = 1
        elif weight < 0:
                weight = -1

        return weight

def duplicate_chars(strings):
    seen = set()
    for string in strings:
        for char in string:
            if char in seen:
                return True
            seen.add(char)
    return False

if __name__ == "__main__":

    combinations = list(combinations(["12", "14", "23", "34"], 2))
    valid = [c for c in combinations if not duplicate_chars(c)]
    print(valid)
