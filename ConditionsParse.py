# Parsing of Divisor Conditions on M(0,6)
# by Connor Halleck-Dube

from ast import literal_eval

f = open("Conditions.txt", "r")
g = open("ConditionsParsed.txt", "w")

# Parses a line from MacCauley2 Output format into an array
def parsed(line):
    l = []
    first = True
    for char in line:
        try:
            if first:
                l.append(literal_eval(char))
                first = False
            else:
                l.append(literal_eval(char))
        except Exception:
            pass
    return tuple(l)

# Determines if a difference of non-identical tuples has no negative elements
def nonNegative(t1, t2):
    if (t1 == t2):
        return False
    if (t1[0] != t2[0]):
        return False
    for i in range(1, len(t1)):
        if (t1[i] - t2[i] < 0):
            return False
    return True

s = set([])
for line in f:
    s.add(parsed(line))

copy = set(s)
for elem in s:
    for comp in s:
        if nonNegative(elem, comp):
            try:
                copy.remove(comp)
            except KeyError:
                pass

for elem in sorted(copy, key = lambda tup: list(reversed(list(tup)))):
    g.write(str(elem))
    g.write("\n")
                
