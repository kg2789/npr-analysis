from momentum_scale import get_triplet
from math import pi
import numpy as np

L = 8

latt = np.array([32,32,32,64])

results = []
for x in range(L):
    for y in range(L):
        for z in range(L):
            for t in range(L):
                p = np.array([x,y,z,t])
                p2 = get_triplet(latt, p)
                if p2 is None:
                    continue
                p_physical = 2 * pi / latt * p
                mu = np.linalg.norm(p_physical)
                results.append((mu,p,p2,p-p2))
                
results = sorted(results, key=lambda x: x[0])
for res in results:
    mu, p1, p2, q = res
    print("{}\t{}\t{}\t{}".format(mu,p1,p2,q))
