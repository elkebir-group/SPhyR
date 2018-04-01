#!/usr/bin/python
import sys

m = int(sys.argv[1])
n = int(sys.argv[2])
k = int(sys.argv[3])

submatrices = m * (m-1) * (m-2) * (n * (n - 1) / 2)

print submatrices * ((k + 1)**4 + 2*k*k*((k+1)**2) + k**4)

