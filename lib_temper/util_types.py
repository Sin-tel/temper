from fractions import Fraction
from typing import TypeAlias, Any, Sequence
import numpy as np

FloatArray: TypeAlias = np.ndarray[Any, np.dtype[np.int64]]

IntArray: TypeAlias = np.ndarray[Any, np.dtype[np.int64]]

Subgroup: TypeAlias = Sequence[int | Fraction]
SubgroupFrac: TypeAlias = Sequence[Fraction]
SubgroupInt: TypeAlias = Sequence[int]
