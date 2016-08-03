# Simplicial Complex Module
# Defines classes and functions for working with simplicial complices
# written for SUMRY 2016
# by Connor Halleck-Dube

from SCM.Helpers import (
    prod,
    nCr,
    nonDegWeight,
    is_zero
)


from itertools import combinations_with_replacement, combinations
from collections_extended import frozenbag, setlist, bag
from sympy import Symbol
from sympy.matrices import *
from sympy.polys import *
from ast import literal_eval
from networkx import Graph, MultiGraph, is_isomorphic, is_connected
from collections import deque
from copy import copy, deepcopy
import operator

## Defines a simplex
class Simplex(frozenbag):
    def __init__(self, l=[]):
        super(Simplex, self).__init__(l)
        
        # degree of the simplex
        self.deg = self._size-1
    
    # Returns true if simplex contains multiset, false otherwise
    def contains(self, multiset):
        return (self >= frozenbag(multiset))
    
    # Maps f onto t
    def reduced(self, f, t):
        return Simplex([i for i in self if (i != f)] + [t for i in self if (i == f)])
    
    # Swaps v1 and v2 in the simplex
    def swapped(self, v1, v2):
        l = ([i for i in self if ((i != v1) and (i != v2))] 
            + [v1 for i in self if (i == v2)])
        if (v1 != v2):
            l += [v2 for i in self if (i == v1)]
        return Simplex(l)
    
    # Performs permutation of vertices
    def vertexPerm(self, perm):
        return Simplex([perm[v] for v in self])
            

## Defines a simplicial complex
class Complex(setlist):
    def __init__(self, l=[]):
        super(Complex, self).__init__([Simplex(i) for i in sorted([list(Simplex(i)) for i in l])])
        try:
            self.d = len(self._list[0])-1
        except IndexError:
            self.d = 0
        self._vertbag = bag()
        for simp in self:
            self._vertbag += simp
        self.vertices = set(self._vertbag)
    
    def n(self):
        return len(self.vertices)
      
    def toStr(self):
        return (str(sorted([sorted(list(simp)) for simp in self])) + "\n")
    
    def display(self):
        print(self.toStr())
        
    def output(self, f):
        f.write(self.toStr())
    
    def write(self, f):
        f.write(str([list(x) for x in list(self)]))
        f.write("\n")
    
    # adds a simplex to the complex
    # Had to do this one manually for some reason
    # Library may have off-by-one error
    def add(self, value):
        simp = Simplex(value)
        if simp in self:
            pass
        else:
            # Do this first in case value isn't Hashable
            self._dict[simp] = len(self)
            self._list.append(simp)
            self._vertbag += simp
            self.vertices |= set(simp)
        
        if (self.d == 0):
            self.d = simp.deg
    
    # removes a simplex from the complex
    def remove(self, simp):
        try:
            super(Complex, self).remove(Simplex(simp))
            self._vertbag -= simp
            self.vertices = set(self._vertbag)
        except Exception:
            pass
    
    # generalized star construction
    def star(self, vlist):
        scomp = Complex([])
        slist = Simplex(vlist)
        for simp in self:
            if (slist <= simp):  # multisubset
                scomp.append(simp-vlist) # multiset minus
        scomp.d = self.d - len(vlist)
        return scomp
    
    # generates the multiplicity matrix
    def _multMatrix(self):
        # generate multiplicities for degree d
        m = zeros(1, len(self))
        i = 1
        for subset in combinations_with_replacement(self.vertices, self.d):
            row = []
            for simp in self:
                e = prod([nCr(simp.count(i), subset.count(i)) for i in set(subset)])
                row.append(e)
            if not is_zero(row):
                m = m.row_insert(i, Matrix([row]))
            i += 1
        m.row_del(0)
        return m
    
    # Returns whether the complex contains a given face
    def containsFace(self, tup):
        for simp in self:
            if simp.contains(tup):
                return True
        return False
    
    # Returns the count of a face within self
    def simpsCount(self, tup):
        i = 0
        for simp in self:
            if simp.contains(tup):
                i += 1
        return i
    
    # returns 1 if i-th vertex and j-th vertex are adjacent
    # 0 otherwise
    def isAdjacent(self, i, j):
        for simp in self:
            if simp.contains([self.vertices, self.vertices]):
                return 1
        return 0
    
    # Tests if complex is balanceable, returns weighting vector if so, 0 vector if not
    ## MAY BE OBSOLETE
    def canBalance(self):
        Null = self._multMatrix().nullspace()
        if (Null == []):
            return zeros(1, 1)
        n = Null[0] 
        for e in Null[1:]:
            n = n.row_join(e)
        
        return nonDegWeight(n)
    
    # Tests if the complex or any subcomplex is balanceable
    def canGenBalance(self):
        m = self._multMatrix()
        return (m.rank() != m.cols)
        
    # tests simple balancing (based on dim of nullspace)
    def canSimBalance(self):
        Null = self._multMatrix().nullspace()
        if (len(Null) > 1):
            return False
        if (Null == []):
            return False
            
        for weight in Null[0]:
            if (weight <= 0.01) and (weight >= -0.01):
                return False
            
        return True
        
    # Convenience
    def Balance(self):
        return self._multMatrix().nullspace()
    
    # Check the isomorphism of complexes
    def is_isomorphic(self, other, certify=False):
        g1 = MultiGraph()
        g2 = MultiGraph()
        g1.add_edges_from(((v,), tuple(f)) for f in self for v in f)
        g2.add_edges_from(((v,), tuple(f)) for f in other for v in f)
        g1.add_edges_from(((-1,), (v,))
                        for v in self.vertices)
        g2.add_edges_from(((-1,), (v,))
                        for v in other.vertices)
        return is_isomorphic(g1, g2)
    
    # Checks complex isomorphism against list l
    # also rejects if self is not connected
    def isoCheck(self, l):
        g = MultiGraph()
        g.add_edges_from(((v,), tuple(f)) for f in self for v in f)
        
        # Guarantee connected complex
        if (not is_connected(g)):
            return False
        
        # Determine left from right side
        g.add_edges_from(((-1,), (v,))
                    for v in self.vertices)
        
        # Check against list
        for comp in l:
            g1 = MultiGraph()
            g1.add_edges_from(((v,), tuple(f)) for f in comp for v in f)
            g1.add_edges_from(((-1,), (v,))
                        for v in comp.vertices)
            if is_isomorphic(g, g1):
                return False
        return True
    
    # Individually checks if complex is connected
    def is_connected(self):
        g = MultiGraph()
        g.add_edges_from(((v,), tuple(f)) for f in self for v in f)
        return is_connected(g)
        
        # Guarantee connected complex
        if (not is_connected(g)):
            return False
    
    # Maps f vertex onto t vertex
    # if t already used, is a pinch
    # otherwise performs a change of variable
    def vertexIdentify(self, f, t):
        l = []
        for simp in self:
            l.append(simp.reduced(f, t))
        return Complex(l)
    
    # Swaps two vertices
    def vertexSwap(self, v1, v2):
        l = []
        for simp in self:
            l.append(simp.swapped(v1, v2))
        return Complex(l)
    
    # Permutes two sets of n vertices
    def vertexPerm(self, flist, tlist):
        if (len(flist) != self.n()) or (len(tlist) != self.n()):
            raise IndexError('Wrong number of vertices for permutation!')
        else:
            # Generates permutation as dictionary
            perm = dict()
            for i in range(self.n()):
                perm[flist[i]] = tlist[i]
                
            # perform permutation
            l = []
            for simp in self:
                l.append(simp.vertexPerm(perm))
            return Complex(l)
            
    
    # Converts to polynomial
    # Returns 0 if not simply balanceable
    def toPoly(self):
        # Test simple balanceability
        #TODO
        
        # Generate the variables for the polynomial
        varstrs = []
        for i in range(1, self.n()+1):
            varstrs.append("x_" + str(i))
        varstrs.sort()
        vars = []
        for varstr in varstrs:
            vars.append(Symbol(varstr))
        
        # generate a monomial for each simplex
        poly = 0
        for simp in self:
            m = 1
            for v in simp:
                m *= vars[v-1]
            poly += m
        print(poly)
        print(factor(poly))
        return poly
    
    # Returns the f-vector (think like degree list) for the complex
    def fvect(self):
        l = [1]
        for d in range(1, (self.d+1)+1):
            c = 0
            for face in combinations_with_replacement(self.vertices, d):
                if self.containsFace(face):
                    c += 1
            l.append(c)
        return tuple(l)
        
                
# Takes the product of two complexes
def ComplexProd(c1, c2):
    prod = Complex([])
    for s1 in c1:
        for s2 in c2:
            prod.add(s1 | s2)
    prod.d = c1.d + c2.d + 1
    return prod

# Generates the complete d-complex on (n-1)-vertices in Mbar(0, n)
# no singularity
def Kcomp(n, d, verts = None):
    if (verts == None):
        verts = list(range(1, n+1))
    simps = [s for s in combinations(verts, d+1)]
    return Complex(simps)

# Singular version
def KcompSing(n, d, verts = None):
    if (verts == None):
        verts = list(range(1, n+1))
    simps = [s for s in combinations_with_replacement(verts, d+1)]
    return Complex(simps)

# Returns the union of two complexes
def Union(c1, c2):
    c = Complex(c1._list + c2._list)
    c._vertbag = c1._vertbag + c2._vertbag
    c.vertices = set(c._vertbag)
    return c
    
    
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
    comp = comp.vertexPerm([1,2,3,d+2,2*d+2,2*d+3], [u,v,smallcol,mid,w,bigcol])
    return comp
    
# Generates the exodia being
# m2 contains mostly middle points
def Exodia(d, u, v, l, m1, m2, r):
    comp = MirCon(2*d+1, d, d)
    for i in range(3, (d+1)+1):
        comp = comp.vertexIdentify(i, 3)
    for i in range(d+3, (2*d+2)+1):
        comp = comp.vertexIdentify(i, d+3)
    for i in range(2*d+3, (3*d+1)+1):
        comp = comp.vertexIdentify(i, d+2)
    comp = comp.vertexPerm([1,2,3,d+2,d+3,3*d+2], [u,v,l,m1,m2,r])
    return comp
    
    
    
    
    
    