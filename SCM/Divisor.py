# Simplicial Complex Divisor Module
# Implements divisor class and associated functions
# written for SUMRY 2016
# by Connor Halleck-Dube

from SCM.Complex import (
    Complex, 
    Stevelike, 
    Kcomp, 
    KcompSing,
    MirCon,
    Exodia
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

BigList = []

## Divisors on Mbar(0, n) and conversions to/from simplicial complexes

# Defines a divisor
class Divisor(dict):
    def __init__(self, n, coeffs = None, model = None):
        if (coeffs == None):
            coeffs = [0 for i in range(exceptionalCount(n))]
        self.coeffs = coeffs
        self.n = n
        super(Divisor, self).__init__(self)
        if (model == None):
            verts = set(range(1, n))
            self.model = n
        else:
            print("IT HAPPENED!")
            print(coeffs, model)
            exit(1)
            verts = list(range(1, n+1))
            verts.remove(model)
            self.model = model
        if (len(verts) != n-1):
            raise ValueError('Incorrect length of vertex names!')
    
        i = 0
        for index in powerset(list(verts)):
            if (len(index) < n-3):
                try:
                    self[index] = coeffs[i]
                except IndexError:
                    self[index] = 0
                except TypeError:
                    self[index] = 0
                
                i += 1
        self.deg = self.coef(tuple()) - 1 
        self.vertices = set(verts)
        
    # Coefficient of element associated with i (a tuple)
    def coef(self, i):
        return self[i]
    
    # Properly sorts the indices of a divisor
    def indexSort(self):
        return sorted(self, key=lambda key: (len(key), key))
    
    # Displays a divisor
    def display(self, newLines = False, compat=False):
        if compat:
            l = [i for i in [(list(i), self[i]) for i in self.indexSort()] if (i[1] != 0)]
        else:
            l = [i for i in [[i, self[i]] for i in self.indexSort()] if (i[1] != 0)]
        if newLines:
            for EI in l:
                print(EI)
        else:
            print(l)
    
    # outputs a divisor
    def output(self, file, newLines = False):
        l = [i for i in [[i, self[i]] for i in self.indexSort()] if (i[1] != 0)]
        if newLines:
            for EI in l:
                file.write(str(EI))
                file.write("\n")
        else:
            file.write(str(l))
            file.write("\n")
    
    # adds two divisors
    def __add__(self, other):
        if (self.n != other.n):
            raise ValueError('Added two divisors from different spaces!')
        elif (self.model != other.model):
            raise ValueError('Added two divisors from different models!')
        else:
            return Divisor(self.n, [self[i] + other[i] for i in self.indexSort()])
    
    # subtract other from self
    def __sub__(self, other):
        if (self.n != other.n):
            raise ValueError('Added two divisors from different spaces!')
        else: 
            return Divisor(self.n, [self[i] - other[i] for i in self.indexSort()])
    
    # Change model from current to ith model
    def ChangeModel(self, i):
        return ChangeCoords(self, i)
            
    
    # Tests to see if self is a sum of other and exceptional divisors
    def constructedBy(self, other):
        for i in self:
            if (self.coef(i) < other.coef(i)):
                return False
        return True
            
    # Checks against list of known divisors for minimality
    def isMinimal(self):
        # Only check the non-positive divisors
        for i in self:
            if (self[i] > 0) and (i != ()):
                return False
            
        for elem in BigList[self.deg]:
            if (elem[()] != self.deg+1):
                print("AHHHHH!")
                exit(1)
            if self.constructedBy(elem) and self != elem:
                elem.display()
                return False
        return True
            
    
    # Reduces a divisor in Mbar(0,n+1) to Mbar(0, n) by forgetting the i-th vertex
    # Auto-converts to 7th model
    def forgetful(self, i):
        l = []
        for index in self.indexSort():
            if i in index:
                if (len(index) == 1):
                    l.append((0-1)*self[index])
                else:
                    l.append(self[index])
        return Divisor(self.n-1, l, list(range(1, self.n)).remove(i))
    
    # converts to complex if possible, empty complex otherwise
    def toComp(self):
        return divToMaxComp(self)
            
    # Tests whether the divisor corresponds to a complex
    def divCorComp(div):
        return (compToDiv(div.n, divToMaxComp(div)) == div)
        
    # Extension operation
    def ext(self):
        return Divisor(self.n, [i-1 for i in self.coeffs], self.model) 
    
    
    # Define a /consistent divisor/ as a divisor such that kI <= kJ for all I containing J
    
    # Tests a divisor for consistency
    def isConsistent(self):
        raise NotImplementedError("I didn't think this would actually be needed!\n")
        
    
    # Increments a divisor 'index', returning another divisor which is consistent if index is 
    def incr(self, index):
        for key in self.indexSort():
            try:
                self[tuple(sorted(index + key))] += 1
            except KeyError:
                pass

    # Decrements a divisor 'index', undoing the action of incr
    def decr(self, index):
        for key in self.indexSort():
            try:
                self[tuple(sorted(index + key))] -= 1
            except KeyError:
                pass
    
    def minimality(self):
        return sum([self[i] for i in self])
    
    # Returns true if the divisor corresponds to a simplicial complex
    def isComplex(self):
        return (compToDiv(self.n, divToMaxComp(self)) == self)
        
    # Returns true if the divisor is effective
    def isEffective(self):
        return (self.toComp().Balance() != [])

    # Returns true if the divisor is strictly effective
    def isStrictEff(self):
        if not self.isEffective():
            return False
        else:
            for key in self.indexSort():
                self[key] -= 1
                if self.isEffective():
                    self[key] += 1
                    return False
                self[key] += 1
            return True

    
    # Permutes the vertices of a divisor in place
    # FROM and TO should be an unzipped dictionary for the permutation
    # E.g. if FROM = [1,2,3] and TO = [2,1,3], then the permutation will take
    #       1 ==> 2
    #       2 ==> 1
    #       3 ==> 3
    # Does error checking to prevent undefined permutations
    # TODO IMPLEMENT
    """
    def vertPerm(self, FROM, TO):
        if (len(FROM) != len(TO)):
            
        
    # Same as vertPerm, but no error checking (for SPEED)
    def _vertPerm(self, FROM, TO):
        for I in self:
    """
            
        
## Conversion Functions

# Generates the complete complex on a divisor
def divToMaxComp(div):
    comp = KcompSing(div.n-1, div.deg)
    
    # removes impossible terms
    for index in div.indexSort():
        for simp in reversed(comp):
            if (sumCount(simp, index) > div.deg+1 + div[index]):
                comp.remove(simp)
    return comp
    
# Generates a divisor from a complex
def compToDiv(n, comp):
    if (n <= len(comp.vertices)):
        raise IndexError('Too many vertices!')
    verts = set(range(1, n))
    div = Divisor(n)
    card = comp.d + 1
    l = [card]
    for i in powerset(verts):
        if (i != tuple()):
            l.append(max([sum([simp.count(v) for v in i]) for simp in comp]) - card)
    return Divisor(n, l)
    



## Searches on the Divisor Space

# Freezes a divisor by converting it to a tuple
def freeze(div):
    l = []
    for I in div.indexSort():
        l.append(div[I])
    return tuple(l)


# Searches for those degree d+1 divisors with n vertices (Mbar(n+1)) 
def DivisorSearch(n, d, outfile):
    # Starting place
    # To be improved
    # starts at all -d's for e_is and e_ijs, -d+1 for all e_ijks
    dstart = Divisor(n+1, [d+1] + [-1 for i in range(exceptionalCount(n+1)-1)])
    dstart.display()
    stricteffs = []
    out = open(outfile, "w")
    checked = dict()
    Q = deque()
    Q.append(dstart)
            
    """
    # Add a divisor and all permutations of vertices to the checklist
    def addChecked(div):
        for perm in 
    """
    
    # recursive step
    def Step(div):
        f = freeze(div)
        
        # check if seen before
        if (f in checked):
            #print("- Seen before")
            return
        #print("- First time seen")
        
        # check a_Is
        for I in div:
            if (div[I] >= 0) and (I != ()):
                return
        
        # test current divisor
        if (div.isEffective()):
            div.display()
            print("- Effective")
            # record div
            checked[f] = True 
            
            # Check the predecessors of div
            for I in div:
                div.decr(I)
                fd = freeze(div)
                if (fd in checked):
                    if (checked[fd]==True):
                        return
                else: 
                    if (div.isEffective()):   
                        checked[fd] = True  
                        return
                    else:
                        Q.append(div)
                div.incr(I)
            
            # if made it this far, div is strictly effective
            div.display()
            stricteffs.append(div)
        else:
            #print("- Not Effective")
            # record div
            checked[f] = False
            
            # continue recursion
            #print("- Adding descendents")
            for I in div:
                div.incr(I)
                Q.append(deepcopy(div))
                div.decr(I)
        
            
    max = 0
    while (len(Q) > 0):
        div = Q.pop()
        #div.display()
        Step(div)
        
        """
        if (max < len(Q)):
            max = len(Q)
            if (max%10 == 0):
                print(max)
        """
        
        #print(len(Q))
    
    print("\n", len(stricteffs))

