from lib_temper import *
import numpy as np

subgroup = [2, 3, 5, 7, 11]
s = log_subgroup(subgroup)
# s = np.array([1, np.log2(5) / 4])
d = s.shape[0]

evalues, evectors = np.linalg.eigh(metric_weil(subgroup))
w = evectors * np.sqrt(evalues) @ np.linalg.inv(evectors)

# w = np.diagflat(1.0 / s)
w = w / np.power(np.linalg.det(w), 1 / w.shape[0])

# w = np.eye(d)


def err(a):
    a = a[np.newaxis, :]
    return h(np.vstack([a, s]))


def h(a):
    v = a @ w
    return np.sqrt(np.linalg.det(v @ v.T))


def badness(a):
    com = h(a[np.newaxis, :]) ** (1 / (d - 1))
    return err(a) * com


def dist(a, b):
    a = a[np.newaxis, :]
    b = b[np.newaxis, :]
    return 1200 * h(np.vstack([a, b])) / (h(a) * h(b))


def bary(r, w):
    x = r @ np.linalg.inv(w)
    # x = r @ solve_diophantine(w, np.eye(d).astype(np.int64))

    x = x / np.sum(x)
    return x


w_mat = np.eye(d, dtype=np.int64)
# w_mat = np.tri(d, d, dtype=np.int64).T
# w_mat = np.array([[15, 24, 35], [12, 19, 28], [7, 11, 16]])

# print(np.linalg.det(w_mat))

for _ in range(100):
    b = bary(s, w_mat)

    i1 = 0
    i2 = 1
    best = 1000000000.0
    for i in range(d):
        for j in range(i + 1, d):
            # this_d = -dist(w_mat[i], w_mat[j])
            # this_d = badness(w_mat[i] + w_mat[j])
            # this_d = err(w_mat[i] + w_mat[j])

            k1, k2 = i, j
            if b[k1] < b[k2]:
                k1, k2 = k2, k1

            this_d = err(w_mat[k1] + w_mat[k2]) - err(w_mat[k2])
            # this_d = badness(w_mat[k1] + w_mat[k2]) - badness(w_mat[k2])

            if this_d < best:
                best = this_d
                i1 = k1
                i2 = k2

    new = w_mat[i2] + w_mat[i1]
    w_mat[i2] = new

    # print(new)
    # print(badness(new))
    print(w_mat)

    # print(dist(new, s))
    # print(ratio(kernel(hnf(np.vstack([w_mat[i2], w_mat[i1]]))), subgroup))
    # print(hnf(np.vstack([w_mat[i2], w_mat[i1]])))
    for i in range(d):
        sub = hnf(np.delete(w_mat, i, axis=0))
        print(ratio(kernel(sub), subgroup))

    # print(solve_diophantine(w_mat, np.eye(d).astype(np.int64)))
    if np.any(w_mat[:, 0] > 500):
        break
