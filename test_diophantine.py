# test diophantine solver with random matrices

import numpy as np
from lib_temper import solve_diophantine, hnf

count = 0
for i in range(1, 10000):
    r1 = np.random.randint(1, 10)
    r2 = np.random.randint(1, 10)
    r3 = np.random.randint(1, 10)

    if r1 < r2:
        continue

    A = np.random.randint(-5, 5, size=(r1, r2))
    X = np.random.randint(-5, 5, size=(r2, r3))

    # system should be solvable
    A_hnf = hnf(A.T)
    idx = np.argwhere(np.all(A_hnf[:] == 0, axis=1))
    if idx.shape[0] > 0:
        continue

    B = A @ X

    X2 = solve_diophantine(A, B)
    assert np.all((X - X2) == 0)
    print(r1, r2, r3, "ok")
    count += 1

print(f"checked {count}")
