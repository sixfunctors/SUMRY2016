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

# Converts an unofficial divisor to an official divisor
def listToDiv(divlist):
    verts = sum([set(i[0]) for i in divlist], set())
    return Divisor(len(verts)+1, [i[1] for i in divlist])

  
# Generates the largest possible complex from a divisor
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
    


    
## Special Generating functions

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



## Premade things

# Checks that the permutations occur in a reasonable order 
# (makes addition of same degree minimals cognizant of commutivity)
def orderCheck(comb, mlist):
    for d in range(max(mlist)+1):
        if not (is_sorted([comb[i] for i in range(len(mlist)) if (mlist[i] == d)])):
            return False
    return True
        
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

# The list of Minimal divisors
# index i contains those with i+1 as coefficient of H (the i complexes)
Minimals = list()

# Generate boundary divisors
B = dict()
for pair in combinations(range(1, 7), 2):
    B[pair] = Boundary(7, pair[0], pair[1])
Minimals.append(list(B.values()))
print("Boundary done!")


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



## Tests if there is any more minimal version of a divisor which is effective
## returns the most minimal version, which is the divisor itself if 
def Reduced(div):
    for i in div:
        div[i] -= 1
        if isComplex(div):
            return Reduced(div)
        div[i] += 1
    return div

"""
# Parses a divisor string, returns divisor
def fromStr(string):
    d = Divisor(7)
    d[tuple()] = literal_eval(string(0))
    while (string != ""):
"""


"""
# Checks that incr and decr are functioning as inverses.
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
"""

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

