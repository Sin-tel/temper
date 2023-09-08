from lib_temper import *
import numpy as np
from typing import List


def badness(T):
    q = T[0][0]

    err = np.max(np.abs(T - q * s))
    com = (q) ** (1 / (d - 1))

    bad = err * com  # * 5 ** (1 / d)
    return bad, err


# subgroup = [2, 3]
subgroup = [2, 3, 5, 7]

s = log_subgroup(subgroup)[np.newaxis, :]

d = s.shape[1]
# print(d)

search_range = (0.5, 1000.0)
# search_range = (7.5, 20.5)

count_1 = 0

l = []

for m1 in Pmaps(search_range, subgroup):
    count_1 += 1
    b, e = badness(m1)
    if np.gcd.reduce(m1.flatten()) == 1:
        if b < 1.0:
            l.append((b, m1.copy()))

# print(len(l))
# print(count_1)

for m in l:
    t = m[1]
    print(f"{m[0]:.4f}, " + " ".join(map(str, t[0])) + ", ")
