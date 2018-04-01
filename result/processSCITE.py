#!/usr/bin/python
import sys

def parseGV(in_stream, m, n):
    f.readline()
    f.readline()
    pi = [-2 for i in range(n+1)]
    for _ in range(n):
        line = f.readline()
        s = line.rstrip(";\n").split()
        source = int(s[0]) - 1
        target = int(s[2]) - 1
        assert 0 <= target < n
        parent = -2
        if source == n:
            parent = -1
        else:
            parent = source
        pi[target] = parent

    samples = [-2 for i in range(m)]
    f.readline()

    for _ in range(m):
        line = f.readline()
        s = line.rstrip(";\n").split()
        source = int(s[0]) - 1
        target = int(s[2][1:])
        parent = -2
        if source == n:
            parent = -1
        else:
            parent = source
        samples[target] = source

    return pi, samples

def printMatrix(pi, samples):
    m = len(samples)
    n = len(pi)
    print m, "# taxa"
    print n-1, "# characters"
    for p in range(m):
        states = [ 0 for c in range(n) ]
        parent = samples[p]
        while parent != -1:
            #print p, parent
            states[parent] = 1
            parent = pi[parent]
        print " ".join(map(str, states[:n-1]))

if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.stderr.write("Usage: %s <SCITE.gv> <NUMBER_OF_TAXA> <NUMBER_OF_CHARACTERS>\n" % sys.argv[0])
        sys.exit(1)

    filename = sys.argv[1]
    m = int(sys.argv[2])
    n = int(sys.argv[3])
    with open(filename) as f:
        pi, samples = parseGV(f, m, n)

    #print pi

    #try:
    printMatrix(pi, samples)
    #except IndexError:
    #  sys.stderr.write("%s\n" % str(sys.argv))
