#! python3

# GenSearch.py
# A script to search for balanceable non-isomorphic complexes
# written for SUMRY 2016
# by Connor Halleck-Dube

# Usage: GenSearch.py n d outfile [debugfile [cur]]

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
from queue import Empty
import sys

def GenSearch(n, d, outfile, debugfile = None):
    verts = set([v for v in range(1, n+1)])
    simps = set([Simplex(s) for s in combinations_with_replacement(verts, d+1)])
    for i in range(1, n+1):
        simps.remove(Simplex([i for x in range(d+1)]))
    
    # Output files
    out = open(outfile, "w")
    try:
        debug = open(debugfile, "w")
    except TypeError:
        pass
    
    bq = deque()
    fq = deque()
    bq.append(Complex())
    HT = CompHT()
    
    # generates a list of non-isomorphic complexes one step away from cur
    def Step(cur, q):
        cands = simps - cur
        for cand in cands:
            cur.add(cand)
            tmp = deepcopy(cur)
            cur.remove(cand)
            if HT.checkAdd(tmp):
                q.append(tmp)
    
    t0 = time.clock()
    qlen = 0
    time_in_step = 0
    
    while (len(bq) > 0):
        if (len(bq) > qlen):
            qlen = len(bq)
        cur = bq.popleft()
        if (len(cur.vertices) == n):
            fq.append(cur)
        else:
            tb = time.clock()
            Step(cur, bq)
            ta = time.clock()
            time_in_step += ta - tb
    
    t1 = time.clock()
    
    print(len(fq))
    print(qlen)
    print(round(t1 - t0, 3))
    print(round(time_in_step, 3))
    HT.Info()
    
    for comp in fq:
        out.write(comp.toStr())
