import numpy as np
from math import pi

# This is a script for determining the physical magnitude of a given lattice
# momentum vector along with finding triplets of non-exceptional momenta

def get_triplet(latt, p):
    l = 1
    for L in latt:
        l = np.lcm(L, l)
    weights = np.array(l / latt, dtype=np.int)
    p_mag = np.sum((p * weights)**2)
    allowed_p = np.arange(1, np.sqrt(p_mag / np.min(weights) + 1))
    allowed_p = np.vstack((allowed_p, -allowed_p)).reshape((-1,), order='F')
    allowed_p = np.concatenate(([0], allowed_p))
    for px in allowed_p:
        for py in allowed_p:
            for pz in allowed_p:
                for pt in allowed_p:
                    p2 = np.array([px,py,pz,pt], dtype=np.int)
                    p2_mag = np.sum((p2 * weights)**2)
                    if p2_mag == p_mag:
                        q = p - p2
                        q_mag = np.sum((q * weights)**2)
                        if q_mag == p_mag:
                            return p2
                    elif p2_mag > p_mag:
                        break
    return None

if __name__ == '__main__':
    latt = input("Lattice dimensions? ")
    latt = np.array([int(x) for x in latt.strip().split()])
    a_inv = float(input("a^{-1} (GeV) "))

    while True:
        p = input("Momentum (lattice units)? ")
        p = np.array([float(x) for x in p.strip().split()])
        p_physical = 2 * pi / latt * p * a_inv
        print("Physical magnitude:", np.linalg.norm(p_physical))
        print("Triplet:", get_triplet(latt, p))

