from lib_temper import *
import numpy as np


def h(a, wm):
    v = a @ wm
    return np.sqrt(np.linalg.det(v @ v.T))


def weight(a):
    return np.diagflat(1 / a)


subgroup = [2, 3, 5]

s = log_subgroup(subgroup)[np.newaxis, :]

w = weight(s)
# w = np.diagflat(np.ones_like(s))

# normalize so det(w) == 1
w = w / np.power(np.linalg.det(w), 1 / w.shape[0])
print(np.linalg.det(w))

t = np.vstack(
    [
        # patent_map(54, subgroup),
        patent_map(19, subgroup),
        patent_map(31, subgroup),
        # patent_map(13, subgroup),
    ]
)

r, d = t.shape

print(hnf(t))
# print(hnf(t) @ w)
print(h(s, w))

omega = r / (d - r)
b = h(np.vstack([t, s]), w)
c = h(t, w) ** omega

badness = b * c / h(s, w)
print("badness:", badness)

# print("efficiency:", "%.3f" % badness)
