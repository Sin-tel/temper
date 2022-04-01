import numpy as np
import re
from lib_temper import *

subgroup = [2,3,5]

m1 = patent_map(19,subgroup)
m2 = patent_map(31,subgroup)

T = np.vstack([m1,m2])
print(T)
H, U = hnf(T, transformation = True)
print(H)
print(U)

print(preimage(T))

subgroup = p_limit(31)

m3 = patent_map(22.8,subgroup)
m4 = np.array([[23, 36, 53, 64, 79, 84, 93, 97, 103, 111, 113]], dtype = np.int64)

print(np.average(m4 / log_subgroup(subgroup)))
print(m3)
print(m4)