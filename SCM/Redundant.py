




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
    

        
        
        
        
        
        
        
        
        
        
        
        
        
