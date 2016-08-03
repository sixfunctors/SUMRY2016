
### NO LONGER NEEDED ###

## Defines a weighting on a complex
class WeightedComplex(dict):
    def __init__(self, comp, wlist):
        self.deg = comp.deg
        i = 0
        for simp in comp:
            if (len(wlist) <= i):
                print("Invalid Weighting!")
                break
            self[simp] = wlist[i]
            i += 1
    
    # Tests whether a weighting on a complex is balanced
    def isBalanced(self):
        # generate set of vertices
        verts = set()
        for simp in self.keys():
            verts = verts | simp
            if (self[simp] == 0):
                print("Degenerate weighting at ", simp, "!")
                return False
        verts = set(verts)
        
        # generate sets of size d
        for subset in combinations_with_replacement(verts, self.deg):
            # sum over simplices
            sum = 0
            for simp in self.keys():
                value = self[simp]
                for i in set(subset):
                    value *= nCr(simp.count(i), subset.count(i))
                sum += value
            print(subset, sum)
            if (sum != 0):
                return False
        return True
    
    # Generates the compatible weighting for the generalized star construction
    def compatibleStar(self, vlist):
        # generate compatible weighting
        clist = list(self.values())
        ## TODO actually get compatible weighting
        return WeightedComplex(Complex(self.keys()).star(vlist), clist)
            
    def display(self):
        for simp in self.keys():
            print(list(simp), " :", self[simp])

## Polycycle investigation
def PolyCycle(n):
    s = range(n)
    d = [[x, x] for x in s]
    s = []
    i = 0
    for c in cycle([n-1] + list(range(n-1))):
        s.append(d[i]+[c])
        i += 1
        if (i >= n):
            break
    i = 0
    for c in cycle(list(range(1, n)) + [0]):
        s.append(d[i] + [c])
        i += 1
        if (i >= n):
            break
    return s
    
def PolyCycleCheck(maxn):
    for i in range(3,  maxn+1):
        if (Complex(PolyCycle(i)).canSimBalance() == Matrix([[0]])):
            print("Fails for n=", i)
            return False
        else:
            print("Succeeds for n=", i)
    return True
            
    
## Parses a text file to remove duplicates
def RemoveDup(rname, wname):
    # accounts for different vertex numbers
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
    
    f = open(rname, "r")
    newlist = setlist([])
    for line in f:
        newlist.add(str(mini(literal_eval(line))))
    
    g = open(wname, "w")
    for comp in newlist:
        g.write(str(comp))
        g.write("\n")
    
## returns linear combination of the columns of nullspace with no zeroes
# if none exists, returns zero 1x1 matrix
def nonDegWeight(nullspace):
    Z = zeros(1, 1)
    
    # returns row index+1 if there is a zero row, 0 otherwise
    def zerotest(vec):
        for i in range(vec.rows):
            if (vec[i] == 0):
                return i+1
        return 0
    
    # tests whether possible
    for i in range(n.rows):
        if (n[i, :] == zeros(1, n.cols)):
            return Z
    
    # Try sum
    v = zeros(n.rows, 1)
    for i in range(n.rows):
        s = sum([n[i, k] for k in range(n.cols)])
        v[i] = s
    
    # loop until one is found
    row = zerotest(v)
    while (row):
        print("Got here with n:", n)
        j = 0
        while (n[i, j] == 0):
            j += 1
        v = v + n[:, j]
        row = zerotest(v)
    return v


## Performs an obsolete implementation of GenSearch
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
        
        
        
        
        
        
        
        
        
        
        
