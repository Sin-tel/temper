from lib_temper import *
import numpy as np
from typing import List

# https://stackoverflow.com/questions/35736823/how-to-save-a-dictionary-of-arrays-to-file-in-numpy


def h(a):
    v = a @ w
    return np.sqrt(np.linalg.det(v @ v.T))


def badness(T):
    r, d = T.shape

    b = h(np.vstack([T, s]))
    c = h(T) ** (r / (d - r))

    return -b * c / h(s), -b


subgroup = [2, 3, 5, 7, 11, 13, 17, 19]

s = log_subgroup(subgroup)[np.newaxis, :]
d = s.shape[1]

w = np.diagflat(1.0 / s)
w = w / np.power(np.linalg.det(w), 1 / w.shape[0])

# w = np.eye(d)


search_range = (5.5, 100.5)
# search_range = (0.5, 800.5)

count_1 = 0

l = []

# for m1 in Pmaps(search_range, subgroup):
for x in range(1, 200):
    m1 = patent_map(x, subgroup)
    count_1 += 1
    b, err = badness(m1)
    if np.gcd.reduce(m1.flatten()) == 1:
        if b < 1.0:
            l.append((err, m1.copy()))

l = sorted(l)

for m in l:
    t = m[1]
    # print(f"{m[0]:.4f}, " + " ".join(map(str, t[0])) + ", ")
    # print(np.average(t / s))
    print(t[0, 0])
