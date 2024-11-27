import numpy as np
from lib_temper import *

# a non-trivial example of solving for optimal tuning

# 7-limit
subgroup = p_limit(7)

# construct matrix of commas
comma_1 = factors((81, 80), subgroup)
comma_2 = factors((225, 224), subgroup)
commas = np.hstack([comma_1, comma_2])

# solve for the mapping
mapping = cokernel(commas)

# mapping will be in HNF
print("Temperament matrix:")
print(mapping)

# constrainted optimization will use octaves by default
solution, _ = cte((mapping, subgroup), weight='tenney')

# solution is in octaves, multiply by 1200 to get cents
print("Generators:")
for s in solution:
    print(f"{s*1200:.3f}")
