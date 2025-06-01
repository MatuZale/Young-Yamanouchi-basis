
import numpy as np
import toqito
import sys
from sage.all import *
from qutip import *
import itertools
from toqito.perms import permutation_operator

d = 3
N = 4 #Add 1
o=[]
for i in range(1,N):
    o.append(i)
two = list(itertools.permutations(o))
two = str(two)
q = eval(two)
l = [list(i) for i in q]
SSTccc = 0
def SplitInLParts(n, l):  #Make partitions
    a = [0 for i in range(n + 1)]
    k = 1
    a[0] = 0
    a[1] = n
    while k != 0:
        x = a[k - 1] + 1
        y = a[k] - 1
        k -= 1
        while x <= y and k < l - 1:
            a[k] = x
            y -= x
            k += 1
        a[k] = x + y
        yield a[:k + 1]

pt = SplitInLParts(N-1,d)
pt1 = [list(x)[::-1] for x in pt] #Generating set of all posibles young diagrams 
pt1.sort()
qqq =0
E = {}
Mu = {}
Dim = {}
for p in pt1: 
    E_ij = {}
    print('====================================================')
    print(p)
    SSTc = SemistandardTableaux(p, max_entry=d).cardinality()
    STc = StandardTableaux(p, max_entry=d).cardinality()
    print ('dim mu ' + str(STc))
    print ('multiplicity mu ' + str(SSTc))
    v = list(p)
    vp = len(p)
    spc = SymmetricGroupRepresentation(v, 'orthogonal')
    perm = [spc.representation_matrix(Permutation(one)) for one in l]
    op_perm = [permutation_operator(d, [x - 1 for x in one]) for one in l]
    for pp in range(0,STc):
            for qq in range(0,STc):
                    E_ij[pp,qq] = np.zeros((d**(N-1),d**(N-1)))
    if STc == 1:
        E_ij = np.zeros((d**(N-1),d**(N-1)))
        for one in range(0,len(l)):
            a = perm[one]
            b = op_perm[one]
            E_ij = E_ij + (STc/factorial(N-1))*a*b
    else:
        for one in range(0,len(l)):
            for j in range(0,STc):
                for i in range(0,STc):
                    a = perm[one][j,i] 
                    b = op_perm[one]
                    E_ij[i,j] = E_ij[i,j] + (STc/factorial(N-1))*a*b
    E[qqq] = E_ij
    Mu[qqq] = SSTc
    Dim[qqq] = STc
    qqq += 1
print(E)
