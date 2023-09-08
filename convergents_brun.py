import numpy as np


def f(x, n):
    a = 1200 * np.log2(x)
    X = np.identity(len(x), dtype=np.int64)

    for _ in range(n):
        i = np.argsort(-a)
        a = a[i]
        X = X[i]
        # print(f"{a}\n{X}\n")
        a[0] -= a[-1]
        X[-1] += X[0]
        print(X[1])


f([2, 3, 5, 7, 11], n=50)
