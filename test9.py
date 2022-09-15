from util import *
import numpy as np
import math

s = p_limit(7)
# s = [2,3,7,11]

# complexity limit
# Q = 1*1200

c = 20
Q = 1200/c

G = metric_tenney(s)
Gd = np.linalg.inv(G)
# print(Gd)

n = len(s)
Gd = np.block([[Gd, np.zeros((n,1))],[np.zeros((1,n)), np.array([[1]])]])
# G = np.block([[G, np.zeros((n,1))],[np.zeros((1,n)), np.array([[1]])]])

# print(Gd)
# exit()
slog = log_subgroup(s)

# print(slog)

slog = slog / slog[0]
# slog = slog[1:]

# print(slog)
# slog = np.mod(slog,1.0)
# print(slog)
# exit()
n = slog.size

m = 1

# print(2 ** ((m+n)/m))
# exit()

print(Q**(-m/n) * (2**((m+n+3)*(m+n)/(4*n))))
# c = 200.0
# c = c / 1200
# Q = c**(-n/m)
# print(Q)
# exit()
P = np.atleast_2d(slog).T

# P = P @ np.ones((1,m))
P = P @ P.T
P = P[:,0:m]
print()
# exit()

# L = np.block([[np.array([[1]]),np.zeros((1,n))],[Q*np.atleast_2d(slog).T,Q*np.eye(n)]])

# L = np.block([[np.eye(m),np.zeros((m,n))],[Q*P,Q*np.eye(n)]])
# L = np.block([[Q*np.eye(n),Q*P],[np.zeros((m,n)),np.eye(m)]])

# L = np.block([[np.eye(n),np.zeros((n,1))],[Q*np.atleast_2d(slog),np.array([[Q]])]])

## find edos
# L = np.block([[np.eye(m),np.zeros((m,n))],[Q*P,Q*np.eye(n)]])
## find commas
L = np.block([[np.eye(n)],[Q*np.atleast_2d(slog)]])

L = np.round(L).astype(np.int64)

# L = antitranspose(L.T)
print(L)
# exit()
R = LLL(L,Gd,delta = 0.99)  # commas
# R = LLL(L,G)
# R = LLL(L,np.eye(n+m))
print("==========")
print(R[0,0])
print(R.T)
print("==========")
print(np.abs(R[0:m,0]))
# print(R[0:m,0])

for c in R.T:
	print(ratio(c[:-1],s))

for c in R.T:
	print(c[:-1])