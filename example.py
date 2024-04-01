# basic examples on how to use the library

import numpy as np
from lib_temper import *

# temperament maps are represented as integer matrices
# subgroups are lists of integers / fractions

# set subgroup to 7-limit
# you can also do `subgroup = [2, 3, 5, 7]` etc.
subgroup = p_limit(7)

# get patent map for 31 edo
p_map = patent_map(31, subgroup)
# calculate badness metric
badness = temp_badness((p_map, subgroup))
print(f"31 edo badness: {badness}")

# calculate the map for the join of 19 and 31, by stacking them:
print("31 & 19:")
t_map = np.vstack([patent_map(19, subgroup), patent_map(31, subgroup)])
print(t_map)
# calculate the canonical form (HNF + saturation)
t_map = canonical(t_map)
print("canonical form:")
print(t_map)

# get the nullspace (commas)
# this is a mtrix of column vectors
commas = kernel(t_map)
print("commas: ")
for comma in commas.T:
    # ratio just converts a vector to a ratio
    print(f"\t{ratio(comma, subgroup)}")

# these commas are bad, we want simple ones!
# you can use LLL reduction to find a 'good' basis:
commas = LLL(commas)
print("commas (LLL reduced):")
for comma in commas.T:
    print(f"\t{ratio(comma, subgroup)}")

# finding generators
gens = preimage(t_map)
# again, we want to simplify them
# because we already have a reduced comma basis, this is easy:
gens = simplify(gens, commas)

print("generators:")
for gen in gens.T:
    print(f"\t{ratio(gen, subgroup)}")
