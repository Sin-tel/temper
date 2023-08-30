from lib_temper import *
import numpy as np
from typing import List


def h(a):
    v = a @ w
    return np.sqrt(np.linalg.det(v @ v.T))


def random_map(subgroup):
    n = np.random.random() * 195 + 5
    return patent_map(n, subgroup)


def badness(T):
    r, d = T.shape

    b = h(np.vstack([T, s]))
    c = h(T) ** (r / (d - r))

    return b * c / h(s)


# subgroup = [2, 3 / 2, 5 / 4, 7 / 4, 11 / 8, 13 / 8, 17 / 16, 19 / 16]
subgroup = [2, 3, 5]


s = log_subgroup(subgroup)[np.newaxis, :]

w = np.diagflat(1.0 / s)
# normalize w so it has det(w) = 1
w = w / np.power(np.linalg.det(w), 1 / w.shape[0])


# search_range = (4.5, 1999.5)
search_range = (4.5, 53.5)
for m1 in Pmaps(search_range, subgroup):
    b = badness(m1)
    if np.gcd.reduce(m1.flatten()) == 1:
        if b < 2.0:
            # m1 = defactored_hnf(m1)

            print(m1)
