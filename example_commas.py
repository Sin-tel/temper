import numpy as np
from fractions import Fraction

from lib_temper import *

subgroup = p_limit(7)

# make list of commas 81/80 and 225/224
# you have to use Fraction instead of e.g. `81/80` as python will make that into a float
commas = [Fraction(81, 80), Fraction(225, 224)]

# construct comma basis matrix by factorizing and stacking
c_basis = np.hstack([factors(c, subgroup) for c in commas])
print("commas basis: ")
print(c_basis)

# find the nullspace to get the temperament map
# you don't need to call canonical form because cokernel is already guaranteed to be saturated
t_map = cokernel(c_basis)

# check that we got septimal meantone
print("canonical form:")
print(t_map)
