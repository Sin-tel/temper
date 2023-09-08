from lib_temper import *
import numpy as np
from typing import List

# https://stackoverflow.com/questions/35736823/how-to-save-a-dictionary-of-arrays-to-file-in-numpy


def h(a):
    v = a @ w
    return np.sqrt(np.linalg.det(v @ v.T))


def rel_err(T):
    return h(np.vstack([T, s]))


def badness(T):
    r, d = T.shape

    b = h(np.vstack([T, s]))
    c = (h(T)) ** (r / (d - r))

    return b * c / h(s)


subgroup = [2, 3, 5, 7, 11, 13]

s = log_subgroup(subgroup)[np.newaxis, :]

w = np.diagflat(1.0 / s)
# normalize w so it has det(w) = 1
w = w / np.power(np.linalg.det(w), 1 / w.shape[0])


# search_range = (4.5, 1999.5)
search_range = (4.5, 199.5)

count_1 = 0

l = []

for m1 in Pmaps(search_range, subgroup):
    count_1 += 1
    b = badness(m1)
    if np.gcd.reduce(m1.flatten()) == 1:
        if b < 2.5:
            l.append((b, m1.copy()))

print(len(l))
print(count_1)

print("==========")
checked = set()
l = sorted(l)

count_2 = 0

l2 = []
for i1 in range(len(l)):
    for i2 in range(i1 + 1, len(l)):
        m1 = l[i1][1]
        m2 = l[i2][1]
        t = np.vstack([m1, m2])
        b = badness(t)
        if b < 1.0:
            t = defactored_hnf(t)
            key = t.tobytes()
            if key not in checked:
                checked.add(key)
                l2.append((h(t), t, m1, m2))


best_rel = 10000000

l2 = sorted(l2)
print(len(l2))
print("==========")
for m in l2:
    t = m[1]
    b = m[0]

    r, d = t.shape

    u = rel_err(t)
    if u < best_rel:
        best_rel = u

        print("====")
        print(b)
        print(rel_err(t))
        print(badness(t))

        print(m[2][0])
        print(m[3][0])

        print(t)
        # print(kernel(t).T)
        # print(ratio(kernel(t), subgroup))

        # print(f"{b:.4f}, " + str(t).replace("\n", "") + ", ")
