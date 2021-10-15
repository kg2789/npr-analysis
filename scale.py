#!/bin/python3

# Small script for computing the scale of a given momentum on a given lattice

import numpy as np
from math import pi

latt = np.array(input("Lattice size? ").split(), dtype=np.int)
p = np.array(input("Momentum? ").split(), dtype=np.float)
ainv = float(input("a^{-1}? "))

p_phys = 2 * pi / latt * p * ainv

print("Scale:", np.sqrt(np.sum(p_phys * p_phys)))
