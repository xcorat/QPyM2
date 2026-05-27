import math
def sant_ceil(val, pow10=5):
    """ return ceil, but ignore the digits beyond `pow10` """
    v = math.floor(val*math.pow(10, pow10))
    return math.ceil(v/math.pow(10, pow10))

