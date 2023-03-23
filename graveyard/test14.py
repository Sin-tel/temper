from util import *
import numpy as np
import math

s = p_limit(7)
# s = [2,3,5,11]

n = len(s)

Gd1 = np.linalg.inv(metric_farey(15, s))
Gd2 = np.linalg.inv(metric_tenney(s))
Gd3 = np.linalg.inv(metric_weil(s))
Gd4 = np.linalg.inv(metric_kweil(7, s))
Gd5 = np.eye(n)
Gd6 = np.linalg.inv(metric_wilson(s))

Gd = Gd4
G = np.linalg.inv(Gd)


# M = np.vstack([patent_map(19,s),patent_map(12,s)])
M = np.vstack([patent_map(x, s) for x in [19, 31]])
# M = np.vstack([patent_map(x,s) for x in [1578,342,270]])

print(M)
M = hnf(M)

print(M)

print("===========")

g = preimage(M)

# print(g)

Gm = np.linalg.inv(M @ G @ M.T)
# print(Gm)

# print()
gens = np.eye(M.shape[0], dtype=np.int64)
gens = LLL(gens, Gm)

tun = lstsq((M, s))[0].T

# print(1200*tun)

# get inverse
H, U = hnf(gens, transformation=True)

Mn = U @ M

g = preimage(Mn)

for i in range(Mn.shape[0]):
    if log_interval(g[:, i], s) < 0:
        Mn[i, :] = -Mn[i, :]
        g[:, i] = -g[:, i]

print(Mn)
print("generators:")

# simplify
commas = LLL(kernel(M), Gd2)
# print(commas)
# print(g)
# for k in g.T:
#     print(ratio(k,s), M @ k, 1200*tun @ M @ k)

for i in range(len(commas.T)):
    commas[:, i] = make_positive(commas[:, i], s)
g = simplify(g, commas, Gd2)

# print(g)
for k in g.T:
    print(ratio(k, s), M @ k, cents((tun @ M @ k)[0]) + "c")
