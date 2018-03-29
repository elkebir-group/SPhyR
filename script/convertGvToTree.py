#!/usr/bin/python
import sys

if len(sys.argv) != 2:
    sys.stderr.write("Usage: %s <SCITE_GRAPHVIZ>\n" % sys.argv[0])
    sys.exit(1)

edges = []
samples = {}
with open(sys.argv[1]) as f:
    f.readline()
    f.readline()

    for line in f:
        if line.startswith("node"):
            break
        s = line.rstrip(";\n").split(" -> ")

        u = int(s[0]) - 1
        v = int(s[1]) - 1

        edges.append((u,v))

    for line in f:
        if line.startswith("}"):
            break

        s = line.rstrip(";\n").split(" -> ")

        u = int(s[0]) - 1
        sample = int(s[1][1:])

        samples[sample] = u

# remap n to -1 (root node)
n = len(edges)
for idx,(u,v) in enumerate(edges):
    if u == n:
        edges[idx] = (-1, v)

for idx,u in enumerate(samples):
    if u == n:
        samples[idx] = -1

print len(edges), "# n"
for (u,v) in edges:
    print u, v

print len(samples), "# m"
for sample in samples:
    print samples[sample]
