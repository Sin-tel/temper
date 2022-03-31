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