#!/usr/bin/python
import sys
from Bio import Phylo
from cStringIO import StringIO

def parseNewickTree(string):
    tree = Phylo.read(StringIO(string), "newick")
    return tree

def label(pi, vertex, edgeLabels, vertexStates, index2Vertex):
    for i in range(len(pi)):
        if vertex == pi[i]:
            for c in range(len(vertexStates[vertex])):
                vertexStates[i][c] = vertexStates[vertex][c]
            for mutation in edgeLabels[i]:
                vertexStates[i][mutation] = 1 - vertexStates[vertex][mutation]
                #print "Yay", "vertex", i, index2Vertex[i], "mutation", mutation
            label(pi, i, edgeLabels, vertexStates, index2Vertex)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.stderr.write("Usage: %s <TREE> <LEAVES> <#CHARACTERS>\n" % sys.argv[0])
        sys.exit(1)

    n = int(sys.argv[3])

    edgeLabels = {}
    with open(sys.argv[1]) as f:
        line = f.readline().rstrip("\n")
        tree = parseNewickTree(line)

        while True:
            line = f.readline()
            if line == "":
                break
            if line.startswith("in"):
                child = line.rstrip("\n")
                f.readline()
                assert child not in edgeLabels
                edgeLabels[child] = []
                line = f.readline().rstrip("\n")
                while len(line.split()) == 1:
                    mutation = int(line)
                    edgeLabels[child].append(mutation)
                    line = f.readline().rstrip("\n")

    leafStates = {}
    with open(sys.argv[2]) as f:
        for line in f:
            s = line.rstrip("\n").split("\t")

            leafStates[s[0]] = map(int, s[1:])

    #print edgeLabels
    #print leafStates

    tree.rooted = True
    network = Phylo.to_networkx(tree)

    # map vertices to integers
    vertexIndex = {}
    index2Vertex = []
    for edge in network.edges():
        if str(edge[0]) not in vertexIndex:
            vertexIndex[str(edge[0])] = len(vertexIndex)
            index2Vertex.append(str(edge[0]))
        if str(edge[1]) not in vertexIndex:
            vertexIndex[str(edge[1])] = len(vertexIndex)
            index2Vertex.append(str(edge[1]))

    pi = [-1 for i in range(len(vertexIndex))]
    for edge in network.edges():
        pi[vertexIndex[str(edge[1])]] = vertexIndex[str(edge[0])]

    #print vertexIndex
    #print edgeLabels
    #print pi

    edgeLabelsArray = [{} for i in range(len(vertexIndex))]
    for edge in network.edges():
        if str(edge[1]) in edgeLabels:
            edgeLabelsArray[vertexIndex[str(edge[1])]] = edgeLabels[str(edge[1])]
        else:
            edgeLabelsArray[vertexIndex[str(edge[1])]] = []

    vertexLabels = [ [0 for c in range(n)] for i in range(len(vertexIndex)) ]
    root = pi.index(-1)
    label(pi, root, edgeLabelsArray, vertexLabels, index2Vertex)

    for leaf in leafStates:
        assert vertexLabels[vertexIndex[leaf]] == leafStates[leaf]
        #print leaf, leafStates[leaf]
        #print index2Vertex[vertexIndex[leaf]], vertexLabels[vertexIndex[leaf]]
        #print index2Vertex[pi[vertexIndex[leaf]]], vertexLabels[pi[vertexIndex[leaf]]]
        #print

    print len(vertexIndex), "#vertices"
    print " ".join(map(str, pi))
    for i in range(len(index2Vertex)):
        label = index2Vertex[i]
        if label in leafStates:
            print int(label[2:]) - 1,
        else:
            print label,
        print " ".join(map(str, vertexLabels[i]))
