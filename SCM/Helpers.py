# Simplicial Complex Helper Functions
# written for SUMRY 2016
# by Connor Halleck-Dube


from functools import reduce
from itertools import (
    chain,
    combinations,
    combinations_with_replacement,
    tee
)

from sympy.matrices import *
from ast import literal_eval
from networkx import Graph, MultiGraph, is_isomorphic
from collections import deque
from copy import copy, deepcopy
import operator

# Tests a row vector for being equal to a zero vector
def is_zero(vec):
    for i in vec:
        if ((i >= 0.01) or (i <= -0.01)):
            return False
    return True

# Tests an iterable for being sorted
def is_sorted(iterable, compare=operator.le):
  a, b = tee(iterable)
  next(b, None)
  return all(map(compare, a, b))


# Integer partitions
def partitions(n):
    a = [0 for i in range(n + 1)]
    k = 1
    y = n - 1
    while k != 0:
        x = a[k - 1] + 1
        k -= 1
        while 2 * x <= y:
            a[k] = x
            y -= x
            k += 1
        l = k + 1
        while x <= y:
            a[k] = x
            a[l] = y
            yield a[:k + 2]
            x += 1
            y -= 1
        a[k] = x + y
        y = x + y - 1
        yield a[:k + 1]

# Product of list
def prod(factors):
    return int(reduce(operator.mul, factors, 1))

# combinations function (n Choose r)
def nCr(n, r):
    return reduce(lambda x, y: x * y[0] / y[1], zip(range(n - r + 1, n+1), range(1, r+1)), 1)

# returns iterable of all subsets of an iterable
def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))

# The number of exceptional divisors in Mbar0, n
def exceptionalCount(n):
    return int(sum([nCr(n-1, i) for i in range(n-3)]))

# Returns the summed maximality of a simplex in an index
def sumCount(simp, index):
    return sum([simp.count(v) for v in index])

# Computes degree list for a complex comp
def degList(comp):
    dl = []
    for vert in comp.vertices:
        dl.append(comp.simpsCount([vert]))
    return tuple(dl)

# Computes 2 degrees of separation degList for a complex comp
def deg2List(comp):
    dl = []
    for vert in comp.vertices:
        l = []
        print(vert, set([i for i in comp.vertices if (i != vert)]))
        for neighbor in set([i for i in comp.vertices if (i != vert)]):
            l.append(comp.simpsCount([neighbor]))
        dl.append(tuple(l))
    return tuple(dl)
            

# Computes the Laplacian matrix for a simplicial complex
## PRESUMES STANDARD VERTEX NUMBERING
def Laplacian(comp):
    D = diag(*degList(comp)) # degree matrix
    A = zeros(len(comp.vertices))
    for i in range(1, len(comp.vertices)+1):
        for j in range(i):
            aij = -comp.isAdjacent(i, j)
            A[i,j] = aij
            A[j,i] = aij
    print(D+A)
    return set((D+A).eigenvals())