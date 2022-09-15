from util import *
import numpy as np
import math
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt

def nullspace(A, atol=1e-13, rtol=0):
    A = np.atleast_2d(A)
    u, s, vh = np.linalg.svd(A)
    tol = max(atol, rtol * s[0])
    nnz = (s >= tol).sum()
    ns = vh[nnz:].conj().T
    return ns


def innerprod(a, b, W):
	return np.dot(a, W @ b)


# number of points in each dimension
nx, ny, nz = (5, 5, 5)

# create points
xv = np.arange(-nx, nx+1)
yv = np.arange(-ny, ny+1)
zv = np.arange(-nz, nz+1)

x, y, z = np.meshgrid(xv, yv, zv)

pts = np.vstack([x.flatten(),y.flatten(),z.flatten()])

####################

# s = p_limit(5)
s = [2,3,5]

n = len(s)

Gd1 = np.linalg.inv(metric_farey(7,s))
Gd2 = np.linalg.inv(metric_tenney(s))
Gd3 = np.linalg.inv(metric_weil(s))
Gd4 = np.linalg.inv(metric_kweil(5,s))
Gd5 = np.eye(n)


Gd = Gd4
G = np.linalg.inv(Gd)


L = np.linalg.cholesky(Gd)

basis = L.T

# aabaaab
# abbabab

# LLsLLLs
# LssLsLs
# LmsLmLs

# LLsLLLs
# sLsLsLs
# mLsLmLs

# LLsLLLs
# ssLsLsL
# MMsMLMs

# LLsLLsL
# LsLLsLs
# LMsLMsM

# MV3a
# C  D  E- F  G  A- B- C
# MV3b
# C  D- E- F  G  A- B- C

# D  E  F+ G+ A  B  C+ D


# Lb = np.eye(len(s), dtype = np.int64)
# Lb = LLL(Lb, Gd)
# print(Lb)

## comma
prat = factors("250/243",s)

# porcupine  250/243
# meantone   81/80
# schisma    32805/32768
# diaschisma 2048/2025
# magic      3125/3072
# wurschmidt 393216/390625

M = cokernel(prat)

# M = LLL(M.T, G).T
# print(M)
proj = nullspace((basis @ prat).T).T

proj = proj @ basis

print(proj @ prat)

g = preimage(M)

print(g)

Bm = proj @ g

Gm = Bm.T @ Bm
print(Gm)


Gm = np.linalg.inv(M @ G @ M.T)
print(Gm)

gens = np.eye(2, dtype = np.int64)

print(np.dot(gens.T[0], Gm @ gens.T[0]))
print(np.dot(gens.T[1], Gm @ gens.T[1]))
print(np.dot(gens.T[1], Gm @ gens.T[0]))
normdot = np.dot(gens.T[1], Gm @ gens.T[0])
normdot = normdot / np.sqrt(np.dot(gens.T[0], Gm @ gens.T[0]) * np.dot(gens.T[1], Gm @ gens.T[1]))
# print(normdot)
print("angle: ", np.degrees(np.arccos(normdot)))
gens = LLL(gens, Gm)

# should be 1
print(integer_det(gens))

# M = 

print(gens)
print(gens.T[0])

print(np.dot(gens.T[0], Gm @ gens.T[0]))
print(np.dot(gens.T[1], Gm @ gens.T[1]))
print(np.dot(gens.T[1], Gm @ gens.T[0]))
normdot = np.dot(gens.T[1], Gm @ gens.T[0])
normdot = normdot / np.sqrt(np.dot(gens.T[0], Gm @ gens.T[0]) * np.dot(gens.T[1], Gm @ gens.T[1]))
# print(normdot)
print("angle: ", np.degrees(np.arccos(normdot)))

# get inverse
H, U = hnf(gens, transformation = True)
print(gens)
print(H)
print(U @ gens)
print(M)
print(U @ M)
g = preimage(U @ M)
print(g)
print(ratio(g.T[0],s))
print(ratio(g.T[1],s))

# exit()
nb = proj @ g
# exit()
pts2 = proj @ pts



# transform the points
prat = basis @ prat
pts = basis @ pts
 

bx = basis.T[0]
by = basis.T[1]
bz = basis.T[2]
o = np.zeros_like(bx)
bx = np.vstack([o,bx]).T
by = np.vstack([o,by]).T
bz = np.vstack([o,bz]).T

prat = np.vstack([o,prat.T]).T

bx2 = proj.T[0]
by2 = proj.T[1]
bz2 = proj.T[2]
o2 = np.zeros_like(bx2)
bx2 = np.vstack([o2,bx2]).T
by2 = np.vstack([o2,by2]).T
bz2 = np.vstack([o2,bz2]).T

# cx = 

# Creating figure
fig = plt.figure(figsize = (10, 7))
ax = plt.axes(projection ="3d", proj_type = 'ortho')

fig2 = plt.figure(figsize = (10, 7))
ax2 = plt.axes()
ax2.scatter(pts2[0], pts2[1], color = "green")
ax2.plot(bx2[0], bx2[1], color='red', alpha=0.8, lw=3)
ax2.plot(by2[0], by2[1], color='red', alpha=0.8, lw=3)
ax2.plot(bz2[0], bz2[1], color='red', alpha=0.8, lw=3)

lim = (-10,10)
ax2.set_xlim(lim)
ax2.set_ylim(lim)

ax.scatter3D(pts[0], pts[1], pts[2], color = "green")
ax.plot(bx[0], bx[1], bx[2], color='red', alpha=0.8, lw=3)
ax.plot(by[0], by[1], by[2], color='red', alpha=0.8, lw=3)
ax.plot(bz[0], bz[1], bz[2], color='red', alpha=0.8, lw=3)
ax.plot(prat[0], prat[1], prat[2], color='blue', alpha=0.8, lw=3)

lim = (-10,10)
ax.set_xlim(lim)
ax.set_ylim(lim)
ax.set_zlim(lim)
ax.set_xlabel('2')

plt.show()