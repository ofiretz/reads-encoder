from __future__ import division, print_function
from math import log, ceil, floor
from mzpy_math import ispower2, isinteger
from sys import stderr

"""
Variable Length Codes
"""

# some codes only encode positive ints, therefore we define
# argument "off" to offset the actual value.
# this is also useful when we want to reserve the smallest code
# as a comma
def unary_encode(n, sep=0, off=0):
    assert n+off > 0
    assert sep in [0, 1]
    n += off
    digit = 1 if sep == 0 else 0
    return "{0}{1}".format( (n-1)*str(digit), sep )

def unary_decode(s, sep=0, off=0):
    assert sep in [0, 1]
    digit = 1 if sep == 0 else 0
    for i, c in enumerate(s):
        if c == str(sep):
            return i+1-off, i+1 # value, num bits read

def unary_length(n, off=0):
    return n+off

# binary representation
def beta_encode(n):
    return "{0:b}".format(n)

def beta_decode(s):
    return int(s, 2)

def beta_length(n):
    return 1 + int(floor(log(n, 2)))

def gamma_encode(n, off=0):
    assert n+off > 0
    n += off
    beta_n = beta_encode(n)
    beta_n_len = len(beta_n)
    return "{0}{1}".format( "0"*(beta_n_len-1), beta_n )

def gamma_decode(s, off=0):
    beta_n_len, unary_len = unary_decode(s, sep=1)
    return beta_decode(s[beta_n_len-1:2*beta_n_len-1]) - off, 2*beta_n_len-1

def gamma_length(n, off=0):
    assert n+off > 0
    n += off
    return 2 * beta_length(n) - 1

def rice_encode(n, m, off=0):
    assert ispower2(m)
    n += off
    k = int(log(m, 2))
    q = int(floor(n / m))
    r = n % m
    return ( "1" * q ) + "0" + "{0:0{1}b}".format(r, k)

def rice_decode(s, m, off=0):
    assert ispower2(m)
    k = int(log(m, 2))
    q = 0
    for i, c in enumerate(s):
        if c == "0":
            break
        q += 1
    r = beta_decode(s[q+1:q+1+k])
    return (q * m + r) - off, q + 1 + k

def rice_length(n, m, off=0):
    assert ispower2(m)
    n += off
    k = int(log(m, 2))
    return 1 + k + int(floor(n/(2**k)))

def golomb_encode(n, m, off=0):
    if ispower2(m):
        return rice_encode(n, m)

    n += off

    q = int(floor(n / m))
    r = n % m
    c = int(ceil(log(m, 2)))

    res = "1" * q + "0"

    if r < 2**c - m:
        res += "{0:0{1}b}".format(r, c-1)
    else:
        res += "{0:0{1}b}".format(r+2**c-m, c)

    return res

def golomb_decode(s, m, off=0):
    if ispower2(m):
        return rice_decode(s, m)
    q = 0
    for i, c in enumerate(s):
        if c == "0":
            break
        q += 1
    
    numbits = q + 1

    c = int(ceil(log(m, 2)))
    # read c-1 next bits
    r = beta_decode(s[q+1:q+c])
    numbits += c-1

    if r >= 2**c - m:
        r = beta_decode(s[q+1:q+c+1]) - (2**c-m)
        numbits += 1

    return (q * m + r) - off, numbits

if __name__ == "__main__":
    print()
    print("running unit tests")

    for i in xrange(1000):
        for off in [1, 2, 3]:
            gamma_i = gamma_encode(i, off=off)
            assert len(gamma_i) == gamma_length(i, off=off)
            decoded, c_len = gamma_decode(gamma_i, off=off)
            assert i == decoded
            assert len(gamma_i) == c_len

        unary_i = unary_encode(i, off=1)
        assert len(unary_i) == unary_length(i, off=1)
        decoded, c_len = unary_decode(unary_i, off=1)
        assert i == decoded
        assert len(unary_i) == c_len

        for m in xrange(2,100):
            for off in [0, 1, 2]:
                if ispower2(m):
                    rice_i = rice_encode(i, m, off)
                    assert len(rice_i) == rice_length(i, m, off)
                    decoded, c_len = rice_decode(rice_i, m, off)
                    assert i == decoded
                    assert len(rice_i) == c_len
                else:
                    golomb_i = golomb_encode(i, m, off)
                    #assert len(golomb_i) == golomb_length(i, off)
                    decoded, c_len = golomb_decode(golomb_i, m, off)
                    assert i == decoded

    print("all passed")
    print()
