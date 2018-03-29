#!/usr/bin/python
import sys

if len(sys.argv) != 2:
    sys.stderr.write("Usage: %s <INPUT>\n" % sys.argv[0])
    sys.exit(1)

with open(sys.argv[1]) as f:
    m = int(f.readline().rstrip("\n").split()[0])
    n = int(f.readline().rstrip("\n").split()[0])
    D = []
    for i in range(m):
        D.append(f.readline().split())

    for j in range(n):
        for i in range(m):
            if i != 0:
                sys.stdout.write(" ")
            sys.stdout.write(D[i][j])
        sys.stdout.write("\n")

