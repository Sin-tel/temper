from util import *

# this trigger an overflow exception in diophantine hermite solve!

subgroup = p_limit(41)

print(subgroup)

M = patent_map(311, subgroup)

print(M)

k = kernel(M)

print(k)
