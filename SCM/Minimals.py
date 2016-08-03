# Minimals.py
# Generates all the known minimal divisors on M

from SCM.Complex import (
    Complex, 
    Stevelike, 
    Kcomp, 
    KcompSing,
    MirCon,
    Exodia
)

from SCM.Divisors import (
    Divisor
)

from SCM.Helpers import (
    exceptionalCount,
    powerset,
    partitions,
    is_sorted,
    sumCount
)

from itertools import (
    combinations_with_replacement, 
    combinations, 
    permutations,
    product
)

from sympy.matrices import *
from ast import literal_eval
from networkx import Graph, MultiGraph, is_isomorphic
from collections import deque
from copy import copy, deepcopy
import operator

## Individual Generating functions

# generates the divisor for the boundary divisor on v and w in Mbar(0, n)
def Boundary(n, v, w, model = None):
    if model != None:
        raise NotImplementedError('Only do this if necessary')
    d = Divisor(n, [1] + [-1 for i in range(exceptionalCount(n)-1)], model)
    d.incr((v,))
    d.incr((w,))
    for key in d:
        if (d[key] > 0):
            d[key] = 0
    d[()] = 1
    return d
            

# generates the bowtie with mid in the center and l1 and l2 on the same side
# and r1 and r2 on the same side
def Hypertree5(n, mid, l1, l2, r1, r2, model = None):
    if (model != None):
        raise NotImplementedError('Only implement this if necessary')
    c = Complex([[mid, l1], [mid, l2], [mid, r1], [mid, r2], [l1, l2], [r1, r2]])
    return compToDiv(n, c)
    
# generates the 6 bowtie
def Hypertree6(n, ml, mr, l1, l2, r1, r2, model = None):
    c = Complex([[ml, l1], [ml, l2], [l1, l2], [ml, mr], [mr, r1], [mr, r2], [r1, r2]])
    return compToDiv(n, c)

# Generates the chain of tetrahedra
def chain(n, i, a, b, c, d, e):
    c = Complex([[i, a, c], [i, c, d], [i, a, d], [a, c, d], [a, b, d], [a, d, e],
            [a, b, e], [b, d, e], [i, b, e], [b, c, e], [i, b, c], [i, c, e]])
    return compToDiv(n, c)

# generates the opie divisor for n=7
def OP(s, k1, k2):
    d = Divisor(7, [7-4])
    
    # specials
    for l in range(1, 7-4):
        for index in d:
            if ((len(index) == l) and (not ((k1 in index) or (k2 in index)))):
                d[index] = -(7-4) + l
        
    # k1s
    for index in d:
        if ((k1 in index) and (not(s in index))):
            d[index] = -1
    
    # k2s
    for index in d:
        if ((k2 in index) and (not(s in index))):
            d[index] = -1
    
    
    #k1 and k2 
    for index in d:
        if (d[index] != 0) and (k1 in index) and (k2 in index):
            d[index] = 0
    
    return d


## Generate all the known minimal divisors on M0,n

# The tiered list of Minimal divisors
# Minimals[i] contains degree i+1 divisors
Minimals = list()

# Generate boundary divisors
B = dict()
for pair in combinations(range(1, 7), 2):
    B[pair] = Boundary(7, pair[0], pair[1])
Minimals.append(list(B.values()))
print("Boundary done!")


# Checks that the permutations occur in a reasonable order 
# (makes addition of same degree minimals cognizant of commutivity)
def orderCheck(comb, mlist):
    for d in range(max(mlist)+1):
        if not (is_sorted([comb[i] for i in range(len(mlist)) if (mlist[i] == d)])):
            return False
    return True
        
# Checks for repeats in a list of divisors
def repeatCheck(l):
    for i in combinations(l, 2):
        if (i[0]==i[1]):
            print("Repeats occurred:")
            i[0].display()
            i[1].display()
            return False
    return True
        
# Generates all combinations of elements from Minimals[i] for i in mlist
def sums(mlist, Minimals):
    l = []
    for comb in product(*[range(len(Minimals[i])) for i in mlist]):
        try:
            if (orderCheck(comb, mlist)):
                l.append(sum([Minimals[mlist[i[1]]][i[0]] for i in zip(comb, range(len(comb)))], Divisor(7)))
        except KeyError:
            print(comb, mlist, Minimals[mlist[0]][comb[0]])
            Minimals[mlist[0]][comb[0]].display()
            exit()
    return l



## d = 1

# Generate 5-bowties (hypertrees)
H5 = list()
for i in permutations(range(1,7), 5):
    if (i[1] < i[2]) and (i[3] < i[4]) and (i[1] < i[3]):
        H5.append(Hypertree5(7, *i))
repeatCheck(H5)
print("5-ties done:", len(H5))

# Generate 6-bowties (hypertrees)
H6 = list()
for i in permutations(range(1, 7)):
    if (i[0] < i[1]) and (i[2] < i[3]) and (i[4] < i[5]):
        H6.append(Hypertree6(7, i[0], i[1], i[2], i[3], i[4], i[5]))
repeatCheck(H6)
print("6-ties done:", len(H6))

# Generate hexagons
Hex = list()
for i in permutations(range(2, 7)):
    if (i[0] < i[1]):
        Hex.append(compToDiv(7, Complex([[1, i[0]], [1, i[1]], [i[1], i[2]], [i[2], i[3]], [i[3], i[4]], [i[0], i[4]]])))
repeatCheck(Hex)
print("Hexagons done:", len(Hex))

Minimals.append(list())
Minimals[1] += H5
Minimals[1] += H6
Minimals[1] += Hex
repeatCheck(Minimals[1])


## d = 2

# Generate the chain of tetrahedras
Chains = list()
for i in permutations(range(2,7)):
    if (i[0] < i[1]) and (i[0] < i[3]) and (i[1] < i[4]):
        Chains.append(chain(7, 1, i[0], i[1], i[2], i[3], i[4]))
print("Chains done:", len(Chains))

# Generate the stevens
Stevens = list()
for i in permutations(range(1, 7)):
    if (i[0] < i[1]) and (i[1] < i[2]):
        Stevens.append(compToDiv(7, Stevelike(2, *i)))
print("Stevens done:", len(Stevens))
repeatCheck(Stevens)

# Generate the opies
opies = list()
for i in permutations(range(1, 7), 3):
    if (i[1] < i[2]):
        opies.append(OP(*i))
print("Opies done:", len(opies))

# Generate the exodias
exo = list()
for i in permutations(range(1, 7)):
    if (i[0] < i[1]) and (i[3] < i[4]) and (i[1] < i[2]):
        exo.append(compToDiv(7, Exodia(2, *i)))
print("Exodias done:", len(exo))
repeatCheck(exo)
        
Minimals.append(list())
Minimals[2] += Chains
Minimals[2] += Stevens
Minimals[2] += opies
Minimals[2] += exo
repeatCheck(Minimals[2])

## d = 3

Minimals.append(list())


## Combinations thereof

print("Minimal lengths")
for d in range(2+1):
    print(len(Minimals[d]))

print("\n****Full lengths****")
BigList.append(Minimals[0])
print(len(BigList[0]))
for d in range(1, 2+1):
    BigList.append(list())
    for part in partitions(d+1):
        print([i-1 for i in part])
        BigList[d] += sums([i-1 for i in part], Minimals)
        print("\t", len(BigList[d]))
    print(len(BigList[d]))
    print("\n")



# Checks that incr and decr are functioning as inverses.
def incrCheck():
    print("Beginning check!")
    for i in range(0, 3+1):
        print(i, ":")
        for div in Minimals[i]:
            for I in div:
                div2 = deepcopy(div)
                div.incr(I)
                div.decr(I)
                if (div2 != div):
                    print("\tNot Clear A!")
                    div2.display()
                    div.display()
                    exit(1)
                
                div.decr(I)
                div.incr(I)
                if (div2 != div):
                    print("\tNot Clear B!")
                    div2.display()
                    div.display()
                    exit(1)
                
        print("\tClear!")