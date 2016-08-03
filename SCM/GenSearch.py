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

"""


def loop1(n, HT, bq, fq):
    while (True):
        try:
            cur = bq.get(True, 10)
            if (len(cur.vertices) == n):
                fq.put(cur) # Move to the next step of the process
            else:
                Step(cur, HT, bq)
        except Empty:
            break
    
    
    
    print("Terminating a process")

## MAIN





    
# Check dimension of cur (if present)
try:
    if (sys.argv[5].d != d):
        raise IndexError('Wrong degree cur!')
except IndexError:
    pass
except AttributeError:
    pass

# Generates those non-isomorphic classes of d complexes on exactly n vertices
# which it encounters _first_, until it gets to n vertices 
if __name__ == '__main__':
    t0 = time.clock()
   
    workers = Pool()
    while True:
        try:
            bq.get(True, 10)
            
    
    t1 = time.clock()
    print("Starts generated in", round(t1-t0, 3), "seconds")
    exit(0)
    while True:
        try:
            out.write(fq.get(False).toStr())
        except Empty:
            exit(0)
    
# Checkpoint everything
debug.write(str(HT) + "\n")
for i in fq:
    debug.write(i.toStr())
debug.write("\n")

# Second step: add vertices by a secondary step process 
# until a simply balanceable complex is reached
# add these to a final list
final = []

def eval2(value):
    if (not value):
        job_server.submit(Step, (cur, HT, fq))

def eval1(value):
    if (value):
        final.append(cur)
        cur.output(out)
    else:
        job_server.submit(cur.canGenBalance, (), callback=eval2)

def can_empty(fq):
    if (len(fq) == 0):
        return False
        
    while (len(fq) > 0):
        # Checkpoint the data every 1000
        if (count%1000 == 0):
            debug.write(str(HT) + "\n")
            for i in fq:
                debug.write(i.toStr())
            debug.write("\n")
        
        
        cur = fq.popleft()
        job_server.submit(cur.canSimBalance, (), callback=eval1)
    return True

while(can_empty(fq)):
    print("This happened.")
    job_server.wait()
    
print(HT.buckLen())
t2 = time.clock()

print(round(t1 - t0, 3))
print(round(t2 - t1, 3))
"""