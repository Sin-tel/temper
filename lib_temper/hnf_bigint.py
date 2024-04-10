# Algorithm adapted from https://github.com/lan496/hsnf
# Works with lists native python integers, as this turns out to be faster for small matrices!
# Also guaranteed to not overflow

import numpy as np

from .util_types import IntMat


def get_pivot(A, i1, j):
    idx = None
    valmin = None

    for i in range(i1, len(A)):
        if A[i][j] == 0:
            continue
        if (valmin is None) or (abs(A[i][j]) < valmin):
            idx = i
            valmin = abs(A[i][j])
    return idx


def _hnf_row(A):
    num_row, num_column = len(A), len(A[0])
    si, sj = 0, 0

    while True:
        if (si == num_row) or (sj == num_column):
            return A

        # choose a pivot
        row = get_pivot(A, si, sj)

        if row is None:
            # if there does not remain non-zero elements, go to a next column
            sj = sj + 1
            continue
        # swap
        A[si], A[row] = A[row], A[si]

        # eliminate the s-th column entries
        for i in range(si + 1, num_row):
            if A[i][sj] != 0:
                k = A[i][sj] // A[si][sj]
                for j in range(0, num_column):
                    A[i][j] -= k * A[si][j]

        # if there does not remain non-zero element in s-th column, find a next entry
        row_done = True
        for i in range(si + 1, num_row):
            if A[i][sj] != 0:
                row_done = False
        if row_done:
            if A[si][sj] < 0:
                for j in range(0, num_column):
                    A[si][j] *= -1

            if A[si][sj] != 0:
                for i in range(si):
                    k = A[i][sj] // A[si][sj]
                    if k != 0:
                        for j in range(0, num_column):
                            A[i][j] -= k * A[si][j]

            si += 1
            sj += 1


def hnf_bigint(A: IntMat) -> IntMat:
    return np.array(_hnf_row(A.tolist()), dtype=np.int64)
