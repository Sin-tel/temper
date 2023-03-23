import math
from fractions import Fraction


# "inverted" farey sequence, not including 0/1 and 1/1
# this should be the same as the "integer limit"
def farey(n):
    l = []
    # We know first two terms are
    # 0/1 and 1/n
    x1 = 0
    y1 = 1
    x2 = 1
    y2 = n

    l.append(Fraction(y2, x2))

    # For next terms to be evaluated
    x = 0
    y = 0
    while True:
        # Using recurrence relation to
        # find the next term
        x = math.floor((y1 + n) / y2) * x2 - x1
        y = math.floor((y1 + n) / y2) * y2 - y1

        if y == 1:
            break

        l.append(Fraction(y, x))

        # Update x1, y1, x2 and y2 for
        # next iteration
        x1 = x2
        x2 = x
        y1 = y2
        y2 = y

    return l
