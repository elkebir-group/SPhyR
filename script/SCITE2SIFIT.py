#!/usr/bin/python
import sys

with open(sys.argv[1]) as f:
    i = 0
    for line in f:
        print i, line,
        i += 1
