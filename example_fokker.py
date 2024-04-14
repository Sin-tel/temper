# example: finding fokker blocks
# this is equivalent to finding a "nice" comma basis that is unimodular

import numpy as np
from lib_temper import *

s = p_limit(11)

# search parameter, higher = smaller commas / larger edos
k = 800.0

B = np.eye(len(s), dtype=np.int64)
W = metric_weil_k(s, k)
W_inv = np.linalg.inv(W)

B = LLL(B, W_inv, delta=0.99)

print("comma basis")
print(B)

print("commas:")
for c in B.T:
    print(f"\t{ratio(c, s)}")


def inv(T):
    # B is guaranteed invertible since the basis we find is unimodular
    K = np.linalg.inv(T)
    return np.round(K).astype(np.int64)


print("dual basis:")
B_inv = np.abs(inv(B))
print(B_inv)
