# generate all combinations (n choose p) ordered by sum of indices


# generate offsets in increasing order (>=)
# to produce a total value
def getOffsets(size,total,maxValue):
    #print(size,total,maxValue)
    if not total: yield [0]*size; return
    if size == 1 and total==maxValue: yield [maxValue]; return
    while total>=0 and size*maxValue>=total:
        for prefix in getOffsets(size-1,total-maxValue,maxValue):
            yield prefix + [maxValue]
        maxValue -= 1

# generate all combinations of a range of values
# that produce a given total
def comboOfSum(total,size,minValue,maxValue):
    if size == 1: yield (total,); return
    base        = list(range(minValue,minValue+size)) # start with smallest(s)
    base[-1]    = min(total-sum(base[:-1]),maxValue)  # end with largest
    maxOffset   = base[-1]-base[-2]-1 # freedom of moving smaller values
    totalOffset = total-sum(base)     # compensate decreasing last
    minLast     = (total + size*(size-1)//2)//size # minimum to reach total
    while base[-1]>base[-2] and base[-1] >= minLast:
        for offsets in getOffsets(size-1,totalOffset,maxOffset):
            yield tuple(b+o for b,o in zip(base,offsets+[0])) # apply offsets
        base[-1]    -= 1 # decrease last value
        totalOffset += 1 # increase total to compensate for decrease
        maxOffset   -= 1 # decrease small values' freedom of movement

# generate combinations in order of target sum  
def comboBySum(size,minValue,maxValue):
    minTotal = minValue*size + size*(size-1)//2
    maxTotal = maxValue*size - size*(size-1)//2
    for total in range(minTotal,maxTotal+1):
        yield from comboOfSum(total,size,minValue,maxValue)