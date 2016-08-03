# generates a list of non-isomorphic complexes one step away from cur
def Step(cur, HT, q):
    cands = simps - cur
    for cand in cands:
        cur.add(cand)
        tmp = deepcopy(cur)
        cur.remove(cand)
        if HT.checkAdd(tmp):
            q.put(tmp)

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
if (len(sys.argv) < 4):
    print("Usage: GenSearch.py n d outfile [debugfile [cur]]")
    exit(1)

n = int(sys.argv[1])
d = int(sys.argv[2])
verts = set([v for v in range(1, n+1)])
simps = set([Simplex(s) for s in combinations_with_replacement(verts, d+1)])
for i in range(1, n+1):
    simps.remove(Simplex([i for x in range(d+1)]))


# Output files
out = open(sys.argv[3], "w")
try:
    debug = open(sys.argv[4], "w")
except IndexError:
    pass

    
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
    bq = mp.Queue()
    fq = mp.Queue()
    bq.put(Complex())
    HT = CompHT()
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
    
"""
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