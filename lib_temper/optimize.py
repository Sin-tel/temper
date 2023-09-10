# generator optimization and metrics

import numpy as np

from .interval import *
from .farey import farey


def lstsq(temp, weight="tenney", V=None):
    M = temp[0]
    s = temp[1]

    M = np.atleast_2d(M)

    j = log_subgroup(s)

    if V is None:
        V = np.eye(len(s))

    if weight == "unweighted":
        G = np.eye(len(s))
    elif weight == "tenney":
        G = metric_tenney(s)
    elif weight == "weil":
        G = metric_weil(s)
    else:
        raise ValueError("unknown weight parameter")

    j = np.atleast_2d(j)

    sol = j @ G @ M.T @ np.linalg.inv(M @ G @ M.T)
    sol = sol.T

    tun = sol.T @ M
    err = tun - j

    return sol, err.flatten()


def cte(temp, weight="tenney", V=None):
    M = temp[0]
    s = temp[1]
    M = np.atleast_2d(M)

    j = log_subgroup(s)

    if weight == "unweighted":
        G = np.eye(len(s))
    elif weight == "tenney":
        G = metric_tenney(s)
    elif weight == "weil":
        G = metric_weil(s)
    else:
        raise ValueError("unknown weight parameter")

    j = np.atleast_2d(j)

    A = M.T
    b = j.T

    r, d = M.shape

    if V is None:
        V = np.zeros((d, 1))
        V[0, 0] = 1

    C = (M @ V).T
    d = (j @ V).T

    Z = np.zeros((C.shape[0], C.shape[0]))

    sol = np.linalg.solve(np.block([[A.T @ G @ A, C.T], [C, Z]]), np.vstack([A.T @ G @ b, d]))

    sol = sol[:r]

    tun = sol.T @ M
    err = tun - j

    return sol.flatten(), err.flatten()


# these metrics are defined on dual space!
# invert them to get the interval norm


# diagonal inverse log primes
def metric_tenney(s):
    j = log_subgroup(s)
    W = np.diag(1.0 / j)
    G = W @ W.T

    return G


# diagonal inverse primes
def metric_wilson(s):
    j = np.array(s).astype(np.double)
    W = np.diag(1.0 / j)
    G = W @ W.T

    return G


# weil norm
def metric_weil(s):
    n = len(s)
    j = log_subgroup(s)
    Bw = np.block([[np.eye(n)], [np.ones((1, n))]])
    Bw2 = Bw @ np.diag(j.flatten())
    Gd = Bw2.T @ Bw2
    G = np.linalg.inv(Gd)
    G = G / G[0, 0]

    return G


# k-weil norm
def metric_weil_k(s, k):
    n = len(s)
    j = log_subgroup(s)
    Bw = np.block([[np.eye(n)], [k * np.ones((1, n))]])
    Bw2 = Bw @ np.diag(j.flatten())
    Gd = Bw2.T @ Bw2
    G = np.linalg.inv(Gd)
    G = G / G[0, 0]

    return G


# farey series (integer limit) norm
def metric_farey(n, s):
    l = farey(n)  # integer limit
    lv = []
    for i in l:
        f = factors_unchecked(i, s)
        if f is not None:
            lv.append(f)

    lv = np.hstack(lv)
    G = lv @ lv.T
    G = G / G[0, 0]

    return G
