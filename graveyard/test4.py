import numpy as np
import re
from lib_temper import *



for i in range(100):
	n = np.random.randint(2,20)
	M = np.random.randint(-10,10,size=(n,n)).astype(np.int64)
	D1 = integer_det(M)
	D2 = int(np.round(np.linalg.det(M)))
	print(D1, D2)
	print(M)
	assert (D1 == D2)

