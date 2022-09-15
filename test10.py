# algorithm for finding nice comma bases

from util import *
import numpy as np
import math
from time import perf_counter


s = p_limit(31)
# s = [2,5,7,11]

# complexity limit
# Q = 1*1200

c = 0.1
Q = 1200/c

G = metric_tenney(s)
# G = metric_weil(s)
Gd = np.linalg.inv(G)
# print(Gd)
# exit()

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

P = np.atleast_2d(slog).T


L = np.block([[np.eye(n)],[Q*np.atleast_2d(slog)]])

L = np.round(L).astype(np.int64)

# L = antitranspose(L.T)
print(L)
# exit()
t1_start = perf_counter()

R = LLL(L,Gd,delta = 0.99)  # commas
t1_stop = perf_counter()
print("Elapsed time:", 1000*(t1_stop - t1_start))
print("==========")
print(R.T)

for c in R.T:
	print(ratio(c[:-1],s))
	print(1200*log_interval(c[:-1],s))

for c in R.T:
	print(c[:-1])

# print(R)
R = R[:-1] 
# print(integer_det(R))
H, U = hnf(R, transformation = True)

print(np.abs(U))