# Graph Isomorphism Hash Table Module
# Implements a hash table for using with complex isomorphism
# written for SUMRY 2016
# by Connor Halleck-Dube


from SCM.Helpers import degList, Laplacian, deg2List
from SCM.Complex import Simplex, Complex

from sympy.matrices import *
from ast import literal_eval
from networkx import Graph, MultiGraph, is_isomorphic
from collections import deque
from copy import copy, deepcopy
import operator
import time

# Hash table for testing complex isomorphism
class CompHT(dict):
    def __init__(self):
        self.time = 0
    
    # adds a complex to self
    def add(self, comp):
        dl = degList(comp)
        try:
            self[dl].append(comp)
        except KeyError:
            self[dl] = [comp]
    
    # Tests if comp is in the hash table up to isomorphism
    # if not, add it and return true
    # if so return false
    def checkAdd(self, comp):
        tb = time.clock()
        key = degList(comp)
        try:
            l = self[key]
            if comp.isoCheck(l):
                self[key].append(comp)
                
                ta = time.clock()
                self.time += ta - tb
                
                return True
            else:
                
                ta = time.clock()
                self.time += ta - tb
                
                return False
        except KeyError:
            if comp.is_connected():
                self[key] = [comp]
                
                ta = time.clock()
                self.time += ta - tb
                
                return True
            else:
                
                ta = time.clock()
                self.time += ta - tb
                
                return False
            
    # return list of bucket lengths
    def Info(self):
        l = []
        for key in self:
            l.append(len(self[key]))
        print(max(l), "size of biggest bucket")
        print(len(l), "buckets")
        print(sum(l), "entries")
        print(round(self.time, 3), "sec spent in HT")