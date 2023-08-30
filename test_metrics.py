from lib_temper import *
import numpy as np


def h(a):
    return np.sqrt(np.linalg.det(a @ a.T))


subgroup = [2, 3 / 2, 5 / 4]

s = log_subgroup(subgroup)[np.newaxis, :]

t = np.vstack(
    [
        patent_map(11, subgroup),
        patent_map(13, subgroup),
    ]
)

r, d = t.shape

print(hnf(t))

b = h(np.vstack([t, s]))
c = h(t) ** (r / (d - r))

print(b * c / h(s))

print("efficiency:", "%.3f" % (h(s) / (b * c)))
