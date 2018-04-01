#!/usr/bin/python
import sys

if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.stderr.write("Usage: %s <k>\n" % sys.argv[0])
        sys.exit(1)

    k = int(sys.argv[1])
    forbidden = set()

    # condition 1
    print "Condition 1"
    for i in range(1, k + 2):
        for i_prime in range(1, k + 2):
            for j in range(1, k + 2):
                for j_prime in range(1, k + 2):
                    print i, 0
                    print 0, j
                    print i_prime, j_prime
                    print
                    M = (i,0,0,j,i_prime,j_prime)
                    assert M not in forbidden
                    forbidden.add(M)

    # condition 2
    print "Condition 2"
    for i in range(1, k + 2):
        for i_prime in range(1, k + 2):
            for j in range(2, k + 2):
                for j_prime in range(1, k + 2):
                    if j_prime == j:
                        continue
                    print i, j_prime
                    print 0, j
                    print i_prime, j
                    print
                    M = (i,j_prime,0,j,i_prime,j)
                    assert M not in forbidden
                    forbidden.add(M)

    # condition 3
    print "Condition 3"
    for i in range(2, k + 2):
        for i_prime in range(1, k + 2):
            if i_prime == i:
                continue
            for j in range(1, k + 2):
                for j_prime in range(1, k + 2):
                    print i, 0
                    print i_prime, j
                    print i, j_prime
                    print
                    M = (i,0,i_prime,j,i,j_prime)
                    assert M not in forbidden
                    forbidden.add(M)

    # condition 4
    print "Condition 4"
    for i in range(2, k + 2):
        for i_prime in range(1, k + 2):
            if i_prime == i:
                continue
            for j in range(2, k + 2):
                for j_prime in range(1, k + 2):
                    if j_prime == j:
                        continue
                    print i, j_prime
                    print i_prime, j
                    print i, j
                    print
                    M = (i,j_prime,i_prime,j,i,j)
                    assert M not in forbidden
                    forbidden.add(M)

    print "Number of forbidden submatrices:", len(forbidden)
