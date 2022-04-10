from util import *
import numpy as np
import math
from farey import farey

### http://x31eq.com/composite.pdf

s = p_limit(5)

# print(s)

l = farey(11)

lv = []
# print(list(map(str, l)))
for i in l:
	f = factors_unchecked(i, s)
	if f is not None:
		lv.append(f)

lv = np.hstack(lv)

# print(lv)

G = lv @ lv.T

print(G)
# G = G / np.sqrt(np.linalg.det(G))
G = G / G[0, 0]

Gd = np.linalg.inv(G)

j = log_subgroup(s)
W = np.diag(1.0 / j)
j = np.atleast_2d(j).T
# print(j)
G2 = W @ W.T

# G2 = G2 /  np.sqrt(np.linalg.det(G2))

Gd2 = np.linalg.inv(G2)

np.set_printoptions(precision=3, suppress=True)

# U = (G2 @ j @ j.T @ G2)/(j.T @ G2 @ j)
U = (G2 @ j @ j.T @ G2)

U = G2 - (U - G2) / (j.T @ G2 @ j)

j_i = 1/j
n = len(s)

U2 = G2 - (j_i @ j_i.T - G2) / n


# print(G - G2)
print(G)
print(U)

Gd_u = np.linalg.inv(U)

# print(G - U)

# weil norm 
n = len(s)
Bw = np.block([[np.eye(n)],[np.ones((1,n))]])
# Bw = np.ones((1,n))
Bw2 = Bw @ np.diag(j.flatten())
Gd3 = Bw2.T @ Bw2
G3 = np.linalg.inv(Gd3)
Gd3 = np.linalg.inv(G3 / G3[0,0])
print("===================")

print(Gd)

print(Gd_u)
print(Gd3)
print(Gd2)


def norm(v, g):
	return np.sqrt(v.T @ g @ v)


# v1 = factors("11/10", s)
# v2 = factors("55/2", s)

v1 = factors("8/5", s)
v2 = factors("5/3", s)

# print(1200*np.log2(ratio(v1, s)))
print(1200*np.sum(j*v1))
print(1200*np.sum(j*v2))
print("======")
print("farey")
print(norm(v1, Gd))
print(norm(v2, Gd))
print(norm(v1, Gd) / norm(v2, Gd))
print("farey asymptotic")
print(norm(v1, Gd_u))
print(norm(v2, Gd_u))
print(norm(v1, Gd_u) / norm(v2, Gd_u))
print("TE")
print(norm(v1, Gd2))
print(norm(v2, Gd2))
print(norm(v1, Gd2) / norm(v2, Gd2))
print("Weil")
print((norm(v1, Gd3)))
print((norm(v2, Gd3)))
print(norm(v1, Gd3) / norm(v2, Gd3))

print("===========")


