from itertools import combinations
import numpy as np
# l = list(range(10))




l = [12,7,31,19,50,20,5]

for combo in comboBySum(3,0,len(l)):
    print([l[i] for i in combo],combo,sum(combo))

# B = list(sorted(combinations(range(minVal,maxVal+1),size),key=sum))