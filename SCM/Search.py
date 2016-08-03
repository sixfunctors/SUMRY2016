# Simplicial Complex Search Module
# Implements brute force searches on simplicial complexes and divisors
# written for SUMRY 2016
# by Connor Halleck-Dube

from SCM.Complex import (
    Simplex, 
    Complex
)

from SCM.Divisor import (
    Divisor,
    Boundary,
    BigList
)

from SCM.IsoHash import (
    CompHT
)

from SCM.Helpers import (
    exceptionalCount,
    powerset,
    Laplacian
)

from itertools import combinations_with_replacement
from sympy.matrices import *
from ast import literal_eval
from networkx import Graph, MultiGraph, is_isomorphic, is_connected
from collections import deque
from copy import copy, deepcopy
import operator
import time
import multiprocessing as mp
from queue import Empty, Full

## Brute Force Searches
def Search(d, numverts, outname):
    verts = [v for v in range(numverts)]
    simps = [s for s in combinations_with_replacement(verts, d+1)]
    for i in range(numverts):
        simps.remove(tuple([i for x in range(d+1)]))
        
    def CouldBeBalanceable(c):
        for multiset in combinations_with_replacement(verts, d):
            if (c.simpsCount(multiset) == 1):
                return False
        return True
    
    f = open(outname, "w")
    comps = []
    for comp in powerset(simps):
        if (len(comp) >= 4):
            c = Complex(comp)
            if (c.vertices != verts) or (not CouldBeBalanceable(c)):
                continue
            w = c.canSimBalance()
            if (w != Matrix([[0]])):
                comps.append(c)
                f.write(str(item))
                f.write("\n")
                print(c)
    
    # remove duplicates
    def mini(comp):
        i = 0
        c = []
        dict = {}
        for simp in comp:
            l = []
            for vert in simp:
                if not(vert in dict.keys()):
                    dict[vert] = i
                    i += 1
                l.append(dict[vert])
            c.append(sorted(l))
        return sorted(c)
    
    # Checks isomorphism against the current list
    def isocheck(testcomp, list):
        for comp in list:
            if Complex(comp).is_isomorphic(testcomp):
                return False
        return True

    f.write("****************************************************\n")
    
    newlist = []
    for comp in comps:
        if isocheck(comp, newlist):
            newlist.append(mini(comp))
    
    for item in newlist:
        f.write(str(item))
        f.write("\n")

# Search based on subkernels of complete d-complex on n verts
def LASearch(d, n, outname):
    verts = [v for v in range(n)]
    simps = [s for s in combinations_with_replacement(verts, d+1)]
    f = open(outname, "w")
    
    # Generate complete multiplicity
    def multMatrixK(d, n):
        for i in range(n):
            simps.remove(tuple(i for x in range(d+1)))
        f.write(str(simps) + "\n")
        return Complex(simps)._multMatrix()
    
    # Find ker Mc
    ker = multMatrixK(d, n).nullspace()
    
    # tests for whether this one is composition of others
    def containsMinor(Comps, I):
        F = [I - set((i,)) for i in I]
        for Iprime in F:
            if (Iprime in Comps):
                return True
        return False
        
    # tests the intersection of the I-plane and ker
    # return true if non-trivial
    def intersect(ker, I):
        card = len(simps)
        n = zeros(card, 1)
        for i in I:
            tmp = [0 for x in range(card)]
            tmp[i] = 1
            n = n.row_join(Matrix(tmp))
        n.col_del(0)
        for b in ker:
            n = n.row_join(b)
        return (n.nullspace() != [])
    
    # Iterate over standard subspaces
    MinComps = set([]) # set of minimally balanceable complexes
    Comps = set([]) # set of balanceable complexes
    f.write("Simplices:\n")
    f.write(list(powerset(list(range(len(simps))))))
    f.write("\n")
    for I in powerset(list(range(len(simps)))):
        print(I)
        I = frozenset(I)
        if containsMinor(Comps, I):
            Comps.add(I)
        else:
            if intersect(ker, I):
                print("******", I)
                MinComps.add(I)
                Comps.add(I)
                f.write(str(min))
                f.write("\n")

# An optimized linear algebra search
def BetterLASearch(d, n, outname):
    verts = [v for v in range(n)]
    simps = [s for s in combinations_with_replacement(verts, d+1)]
    f = open(outname, "w")
    
    # Generate complete multiplicity
    def multMatrixK(d, n):
        for i in range(n):
            simps.remove(tuple(i for x in range(d+1)))
        return Complex(simps)._multMatrix()
    
    # Find ker Mc
    Mc = multMatrixK(d, n)
    ker = Mc.nullspace()
    
    # tests for whether this one is composition of others
    def containsMinor(Comps, I):
        F = [I - set((i,)) for i in I]
        for Iprime in F: 
            if (Iprime in Comps):
                return True
        return False
        
    # tests the intersection of the I-plane and ker
    # return true if non-trivial
    def intersect(ker, I):
        n = zeros(Mc.rows, 1)
        for i in I:
            n = n.row_join(Mc[:, i])
        n.col_del(0)
        return (n.nullspace() != [])
    
    # Iterate over standard subspaces
    MinComps = [] # set of minimally balanceable complexes
    Comps = set([]) # set of balanceable complexes
    f.write("Simplices:\n")
    f.write(str(simps) + "\n")
    f.write("\n")
    for I in powerset(list(range(len(simps)))):
        I = frozenset(I)
        if containsMinor(Comps, I):
            Comps.add(I)
        else:
            if intersect(ker, I):
                c = Complex([ simps[i] for i in I ]) 
                print("******", c)
                MinComps.append(c)
                Comps.add(I)
                f.write(c.toStr())
    



    
    
    