from fractions import Fraction
from typing import TypeAlias, Any, Sequence
import numpy as np

FloatArray: TypeAlias = np.ndarray[Any, np.dtype[np.int64]]
FloatVec: TypeAlias = np.ndarray[tuple[int], np.dtype[np.float64]]
FloatMat: TypeAlias = np.ndarray[tuple[int, int], np.dtype[np.float64]]

IntArray: TypeAlias = np.ndarray[Any, np.dtype[np.int64]]
IntVec: TypeAlias = np.ndarray[tuple[int], np.dtype[np.int64]]
IntMat: TypeAlias = np.ndarray[tuple[int, int], np.dtype[np.int64]]

Subgroup: TypeAlias = Sequence[int | Fraction]
SubgroupFr: TypeAlias = Sequence[Fraction]
SubgroupInt: TypeAlias = Sequence[int]
