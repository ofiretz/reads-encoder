from math import log

def isinteger(n):
    return n == int(n)
    
def ispower2(n):
    return isinteger(log(n, 2))
