from util import *
import numpy as np
import math
from farey import farey
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt

# number of points in each dimension
nx, ny, nz = (2, 1, 1)

# create points
xv = np.arange(-nx, nx+1)
yv = np.arange(-ny, ny+1)
zv = np.arange(-nz, nz+1)

x, y, z = np.meshgrid(xv, yv, zv)

pts = np.vstack([x.flatten(),y.flatten(),z.flatten()])

####################

s = p_limit(5)

## farey norm
l = farey(6) # integer limit
lv = []
for i in l:
	f = factors_unchecked(i, s)
	if f is not None:
		lv.append(f)

lv = np.hstack(lv)
G = lv @ lv.T
G = G / G[0, 0]
Gd = np.linalg.inv(G)

## TE norm
j = log_subgroup(s)
W = np.diag(1.0 / j)
j = np.atleast_2d(j).T
G2 = W @ W.T
Gd2 = np.linalg.inv(G2)

# weil norm 
n = len(s)
Bw = np.block([[np.eye(n)],[np.ones((1,n))]])
# Bw = np.ones((1,n))
Bw2 = Bw @ np.diag(j.flatten())
Gd3 = Bw2.T @ Bw2
G3 = np.linalg.inv(Gd3)
Gd3 = np.linalg.inv(G3 / G3[0,0])

### HERE you can set the norm used to draw the image
# farey norm = Gd
# TE norm = Gd2
# weil norm = Gd3

L = np.linalg.cholesky(Gd3)

basis = L.T

## extra vector to draw in blue
prat = factors("5/4",s)

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
print(bx[0])

prat = np.vstack([o,prat.T]).T


# Creating figure
fig = plt.figure(figsize = (10, 7))
ax = plt.axes(projection ="3d")
 
# Creating plot
ax.scatter3D(pts[0], pts[1], pts[2], color = "green")
ax.plot(bx[0], bx[1], bx[2], color='red', alpha=0.8, lw=3)
ax.plot(by[0], by[1], by[2], color='red', alpha=0.8, lw=3)
ax.plot(bz[0], bz[1], bz[2], color='red', alpha=0.8, lw=3)
ax.plot(prat[0], prat[1], prat[2], color='blue', alpha=0.8, lw=3)

lim = (-3,3)
ax.set_xlim(lim)
ax.set_ylim(lim)
ax.set_zlim(lim)
ax.set_xlabel('2')

plt.show()