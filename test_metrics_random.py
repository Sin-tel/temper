from lib_temper import *
import numpy as np
from typing import List


def h(a):
    v = a @ w
    return np.sqrt(np.linalg.det(v @ v.T))


def random_map(subgroup):
    n = np.random.random() * 195 + 5
    return patent_map(n, subgroup)


# subgroup = [2, 3 / 2, 5 / 4, 7 / 4, 11 / 8, 13 / 8, 17 / 16, 19 / 16]
subgroup = [2, 3, 5, 7, 11, 13]


s = log_subgroup(subgroup)[np.newaxis, :]

w = np.diagflat(1.0 / s)
# normalize w so it has det(w) = 1
w = w / np.power(np.linalg.det(w), 1 / w.shape[0])

b_list = []

n_sample = 10000

for i in range(n_sample):
    t = np.vstack(
        [
            random_map(subgroup),
            random_map(subgroup),
            # random_map(subgroup),
            # random_map(subgroup),
            # random_map(subgroup),
            # random_map(subgroup),
        ]
    )

    r, d = t.shape

    b = h(np.vstack([t, s]))
    c = h(t) ** (r / (d - r))

    badness = b * c / h(s)
    if badness > 1e-7:
        # if badness < 1.0:
        #     print(hnf(t))
        b_list.append(b * c / h(s))

b_list = np.array(sorted(b_list)).astype(np.float64)
l = b_list.shape[0]
# print(l)

count = np.sum(b_list < 1.0)
print(count)

print(100 * count / l)
