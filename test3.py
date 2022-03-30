import numpy as np
import re
from lib_temper import *

subgroup = [2,3,5,7,11]
# M = np.array([[ 2, 1,3 , 4 , 8],[ 0 , 4, 3 , 3 ,-2]],dtype = np.int64)
M = np.array([[ 2, 1 , 4 , 8],[ 0 , 4 , 3 ,-2]],dtype = np.int64)
# M = np.array([[ 2, 1 , 3 , 4],[ 0 , 4 , 3 ,3]],dtype = np.int64)

it = np.array([[1,0],[0,1]], dtype = np.int64)



solve_diophantine(M,it)