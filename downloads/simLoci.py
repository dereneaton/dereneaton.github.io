###############################################################
##  Simulate sequence data on an input tree w/ introgression ##
##  output in .loci readable format                          ##
##  author:  Deren Eaton                                     ##
##  contact: deren.eaton@yale.edu                            ##
##  date:    5/20/14                                         ##
##  verstion: 1.0                                            ##
###############################################################

## load modules                                        
import egglib
import sys
import numpy

## read in arguments
Nloci = int(float(sys.argv[1]))
locuslength = int(float(sys.argv[2]))
N = int(float(sys.argv[3]))
tree = str(sys.argv[4])
migration = str(sys.argv[5])

## split multiple migrations
migscenarios = [m.lstrip("[").rstrip("]") for m in migration.split('/')]

## fixed Bio data
mu = 1e-9  ## per site mutation rate/gen
lu = (mu)*locuslength ## per loc mutation rate/gen
theta = 4*N*lu

## if tree, parse it.
if tree == '0': tree = 0
if tree:
    tre = egglib.Tree(string=tree)
    tiptax = tre.all_leaves()
else:
    tiptax = list("ABCDEFGHI")    
    "clade 1 = ((A,B),C)"
    "clade 2 = ((D,E),F)"
    "clade 3 = ((G,H),I)"
    nodeABCDEFGHI = 3.0 
    nodeABCDEF = 1.5  
    nodeGHI   =  0.75 
    nodeDEF   =  0.75 
    nodeABC   =  0.75 
    nodeGH    =  0.5 
    nodeDE    =  0.5 
    nodeAB    =  0.5 


## function to return node height
def lengthtotip(node):
    dec = 0.
    while node.descendants():
        dec += node.branch_to(node.descendants()[0])
        node = node.descendants()[0]
    return dec


## draw tree in ascii
#    from Bio import Phylo
#    draw = Phylo.read(StringIO.StringIO(str(tree)),'newick')
#    draw.ladderize()
#    Phylo.draw_ascii(draw)


## sets the two parameter classes
paramSet = egglib.simul.CoalesceParamSet(singleSamples=None,
                                         doubleSamples=[1]*len(tiptax),
                                         M=0.0)

## traverse tree fusing populations at each node
if tree:
    for node in tre:
        if node.descendants():
            date = lengthtotip(node)
            s1 = min([tiptax.index(l) for l in node.descendants()[0].leaves_down()])
            s2 = min([tiptax.index(l) for l in node.descendants()[1].leaves_down()])
            paramSet.populationFusion(date,s1,s2)
else:
    "clade 1"
    paramSet.populationFusion(nodeABC, 0,2)     
    paramSet.populationFusion(nodeAB, 0,1)      
    "clade 2"
    paramSet.populationFusion(nodeDEF, 3,5)     
    paramSet.populationFusion(nodeDE, 3,4)      
    "clade 3"
    paramSet.populationFusion(nodeGHI, 6,8)   
    paramSet.populationFusion(nodeGH, 6,7)    
    "together and outgroup"
    paramSet.populationFusion(nodeABCDEFGHI, 0,6)  
    paramSet.populationFusion(nodeABCDEF, 0,3)     


## sets migration patterns 
for mig in migscenarios:
    p2,p3,s,e,m = mig.split(',')
    p2 = tiptax.index(p2)
    p3 = tiptax.index(p3)
    M = 4.*N*float(m)
    paramSet.changePairwiseMigrationRate(0.0, int(p2), int(p3), 0.)
    paramSet.changePairwiseMigrationRate(float(s), int(p2), int(p3), float(M))
    paramSet.changePairwiseMigrationRate(float(e), int(p2), int(p3), 0.)

mutator = egglib.simul.CoalesceFiniteAlleleMutator(theta=theta,
                                                   alleles= 4,
                                                   randomAncestralState=True)
mutator.setSites(locuslength)
aligns = egglib.simul.coalesce(paramSet, mutator, Nloci)



def hetero(n1,n2):
    """
    returns IUPAC symbol for ambiguity bases,
    used for polymorphic sites.
    """
    D = {('G','A'):"R",
         ('G','T'):"K",
         ('G','C'):"S",
         ('T','C'):"Y",
         ('T','A'):"W",
         ('C','A'):"M"}
    a = D.get((n1,n2))
    b = D.get((n2,n1))
    if a:
        return a
    else:
        return b


def unstruct(amb):
    amb = amb.upper()
    " returns bases from ambiguity code"
    D = {"R":["G","A"],
         "K":["G","T"],
         "S":["G","C"],
         "Y":["T","C"],
         "W":["T","A"],
         "M":["C","A"],
         "A":["A","A"],
         "T":["T","T"],
         "G":["G","G"],
         "C":["C","C"],
         "N":["N","N"],
         "-":["-","-"]}
    return D.get(amb)


def makehetero(seq1,seq2):
    seq = ""
    for base in zip(seq1,seq2):
        if base[0] != base[1]:
            seq += hetero(base[0],base[1])
        else:
            seq += base[0]
    return seq


def snpstring(array):
    snpsite = " "*9
    for site in array.transpose():
        reals = site.tolist()
        if len(set(reals)) > 1:
            " convert ambigs to reals"
            for i in xrange(len(reals)):
                if reals[i] in list("RWMSYK"):
                    for j in unstruct(reals[i]):
                        reals.append(j)
            reals = [i for i in reals if i not in "RWMSYK"]
            if sorted([reals.count(i) for i in set(reals)], reverse=True)[1] > 1:
                snpsite += "*"
            else:
                snpsite += "-"
        else:
            snpsite += " "
    snpsite += "|"
    return snpsite


" append names"
names = []
for name in tiptax:
    for i in range(1):
        names.append(name) #+str(i))
        names.append(name) #+str(i))

" appends names for allele 1 vs. allele 2 for each diploid sample "
for i in aligns:
    for s,j in zip(i,names):
        s.name = j
    
" make dictionary with list of loci for each sample "
D = {}
for i in set(names):
    D[i] = []
for samp in D:
    for i in aligns:
        l = []
        for j in i:
            if j.name == samp:
                l.append(j.sequence)
        D[samp].append(l) 


" print in .loci readable format used by pyRAD " 
nn = D.keys()
nn.sort()

if max([len(i) for i in nn]) > 9:
    print 'error: name lengths too long'; sys.exit()

    
for i in range(Nloci):
    l = []
    for n in nn:
        l.append(list(makehetero(D[n][i][0],D[n][i][1])))
        print ">"+n + " "*(10-len(n))+ makehetero(D[n][i][0],D[n][i][1])
    N = numpy.array(l)
    print "//"+snpstring(N)


