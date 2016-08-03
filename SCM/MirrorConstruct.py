# MirrorConstruct.py
# Tests the behavior of MirrorConstruct functions and their divisors
# Written for SUMRY 2016
# By Connor Halleck-Dube

from SCM.Complex import (
    Complex,
    Kcomp,
    Union
)

from SCM.Divisor import (
    Divisor,
    compToDiv,
    Boundary
)
    

# Creates the mirror d-complex on two Kn's, sharing k vertices
def MirCon(n, d, k):
    c1 = Kcomp(n, d, range(1, n+1))
    c2 = Kcomp(n, d, range(n+1-k, 2*n+1-k))
    return Union(c1, c2)

# Tests all the mirror comps for n = 2d+1, k=d, d<=dmax
def MirConSearch(dmax):
    for d in range(1, dmax+1):
        m = MirCon(2*d+1, d, d)
        if (m.canSimBalance()):
            print((2*d+1,d,d))

# Generates the stevelike being
def Stevelike(d, u, v, smallcol, mid, w, bigcol):
    comp = MirCon(2*d+1, d, d)
    for i in range(3, (d+1)+1):
        comp = comp.vertexIdentify(i, 3)
    for i in range(d+2, (2*d+1)+1):
        comp = comp.vertexIdentify(i, d+2)
    for i in range(2*d+3, (3*d+2)+1):
        comp = comp.vertexIdentify(i, 2*d+3)
    comp = comp.vertexSwap(1, u)
    comp = comp.vertexSwap(2, v)
    comp = comp.vertexSwap(3, smallcol)
    comp = comp.vertexSwap(d+2, mid)
    comp = comp.vertexSwap(2*d+2, w)
    comp = comp.vertexSwap(2*d+3, bigcol)
    return comp