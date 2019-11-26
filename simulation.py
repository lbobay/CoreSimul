
# Parameters to implement: 

########################################################################
import math

 
def mean( echantillon ) :
    size = len( echantillon )
    moyenne = float(sum( echantillon )) / float(size)
    return moyenne


def stat_variance( echantillon ) :
    n = float(len( echantillon )) # taille
    mq = mean( echantillon )**2
    s = sum( [ x**2 for x in echantillon ] )
    variance = s / n - mq
    return variance


def stat_ecart_type( echantillon ) :
    variance = stat_variance( echantillon )
    ecart_type = math.sqrt( variance )
    return ecart_type

def median( echantillon) :
        echantillon.sort()
        size = len( echantillon )
        if len( echantillon ) % 2 == 0:
                M= float(echantillon[size / 2 - 1] + echantillon[size / 2]) / 2
        else:
                M= echantillon[size / 2]
        return M




########################################################################



print("USAGE:")

# Length: provided by user
L = 22000
# GC%
GC=50
if GC < 1:
	print("WARNING: the GC% option must be provided in % (expected value ranges betwenn 0 and 100)")
# Branch re-scaling 
coeff = 1
# Kappa
kappa = 1

################# Recombination parameters #################
# Size of recombination events follow an exponential distribution of size delta
# Recombination events occur at rate R (COEFF) which is a function of the mutation rate
DELTA = 100
COEFF = 1

import numpy as np
import numpy
import sys
import random
import os

print(sys.argv)

path = sys.argv[-1]
coeff = float(sys.argv[-2])
DELTA = int(sys.argv[-3])
COEFF = float(sys.argv[-4])
L = int(sys.argv[-5])
GC = float(sys.argv[-6])
kappa = float(sys.argv[-7])
sub = sys.argv[-8].split(",")



un,deux,trois= sub[0].strip(" "),sub[1].strip(" "),sub[2].strip(" ")
un,deux,trois = float(un),float(deux),float(trois)
S = un + deux + trois
if S>= 0.99 and S <= 1.01:
	pass
elif S==100:
	un,deux,trois = int(un/100),int(deux/100),int(trois/100)
else:
	print("PROBLEM: sum of codons frequency is not 1 or 100. Resetting to default: 0.2,0.1,0.7")
	un,deux,trois=0.2,0.1,0.7


un,deux,trois = int(round(100*float(un),0)) , int(round(100*float(deux),0)) , int(round(100*float(trois),0))




comp={}
alpha=["A","T","C","G"]
transition={}
transition["A"]="G"
transition["G"]="A"
transition["C"]="T"
transition["T"]="C"
other1,other2={},{}
other1["A"],other2["A"]="C","T"
other1["T"],other2["T"]="G","A"
other1["C"],other2["C"]="G","A"
other1["G"],other2["G"]="C","T"


# Creates the substitution matrix with Kappa
probability={}
if kappa ==1:
	for N in alpha:
		toto = list(alpha)
		toto.remove(N)
		probability[N] = toto
else:
	for N in alpha:
		probability[N]=[]
		i=0
		while i < round(10 * kappa,0):
			probability[N].append(transition[N])
			i+=1
		i=0
		while i < 5:
			probability[N].append(other1[N])
			probability[N].append(other2[N])
			i+=1



for N in alpha:
	print("Probability:",N," ",probability[N]," ",len(probability[N]))




CODONS = [un,un + deux,100]

ALL_POS = range(100)

# MODIFY THIS FUNCTION TO INCLUDE CODON BIAS
def single(seq):
	i = random.choice(range(len(seq)))
	frame = random.choice(ALL_POS) + 1
	if frame <= CODONS[0]:
		F=1
	elif frame <= CODONS[1]:
		F=2
	elif frame > CODONS[1]:
		F=0
	if i%3==F:
		pass
	elif (i+1)%3 == F:
		i=i+1
	elif (i+2)%3 == F:
		i=i+2
	try:
		N = seq[i]
	except IndexError:
		i=i-3
		N = seq[i]
	new = random.choice(probability[N])
	out = seq[:i] + new + seq[i+1:]	
	return out


def mutate(seq,val):
	L = len(seq)
	mutations = numpy.random.poisson(L * val)
	out=seq
	selection = range(len(out))
	j=0
	while j <= mutations:
		i = random.choice(selection)
		N = out[i]
		new = random.choice(probability[N])
		out = out[:i] + new + out[i+1:]	
		j+=1
	return out



type = {}
branch = {}
f=open(path + "renamed.txt","r")
for l in f:
	a=l.strip("\n").split("\t")
	node = a[0]
	length = float(a[2])
	branch[node] = length * coeff
	type[node] = a[3]

f.close()


parent={}
dicho={}
f=open(path + "dichotomies.txt","r")
for l in f:
	a=l.strip("\n").split("\t")
	dicho[a[0]] = [a[1],a[2]]
	parent[a[1]] = a[0]
	parent[a[2]] = a[0]

f.close()



print("Initialization")
roots=[]
f=open(path + "roots.txt","r")
for l in f:
	roots.append(l.strip("\n"))

f.close()

GC = round(GC * 100,0)
AT = 10000-GC
AT = round(AT,0)

print(GC," ",AT)
ALPHA=[]
i=0
while i < GC:
	ALPHA.append("G")
	ALPHA.append("C")
	i+=1

i=0
while i < AT:
	ALPHA.append("A")
	ALPHA.append("T")
	i+=1


print("ALPHA= ", len(ALPHA))



tmp=[]
while len(tmp) < L:
	tmp.append(random.choice(ALPHA))

seq = "".join(tmp)


sequence={}
for node in roots:
	sequence[node] = seq 

print("Traceback")



tmp = list(roots)


########## Simulate evolution without recombination

cumul={}
for node in roots:
	cumul[node] = branch[node]


if 1==1:
	#h=open("tips.fa","w")
	while len(tmp)>0:
		node = tmp[0]
		if type[node] != "tip":
			node1,node2 = dicho[node][0],dicho[node][1]
			sequence[node1],sequence[node2] = mutate(sequence[node] , branch[node1]), mutate(sequence[node] , branch[node2])
			tmp.append(node1)
			tmp.append(node2)
			if node in cumul:
				cumul[node1] = cumul[node] + branch[node1]
				cumul[node2] = cumul[node] + branch[node2]
			else:
				cumul[node1] =  branch[node1]
				cumul[node2] =  branch[node2]	
			tmp.remove(node)
		else:
			#print node," is a tip"
			new = mutate(sequence[node] , branch[node])
			tmp.remove(node)
			#h.write(">" + node + "\n" + new + "\n")
	#h.close()



#print "Exiting"

#exit()

#for node in cumul:
#	print node ," ",cumul[node]," ",branch[node]

rev={}
for node in cumul:
	if cumul[node] in rev:
		rev[cumul[node]].append(node)
	else:
		rev[cumul[node]]=[node]


values = list(rev.keys())

values.sort()

#print values

# Build the dictionary interval to define which branches overlapped in time
interval={}
for node in cumul:
	interval[node] = [cumul[node] - branch[node] , cumul[node] ]
	#print node ," ",cumul[node]," to ",cumul[node] - branch[node]

# Creat new internal nodes

MIN,MAX=10,0
for node in cumul:
	if cumul[node]>MAX:
		MAX = cumul[node]
	if branch[node] != 0:
		if branch[node] < MIN:
			MIN= branch[node]

print("MIN= ",MIN)
print("MAX= ",MAX)
scan = {}
nb=0
i = 0
while i < MAX:
	tmp=[]
	for node in interval:
		deb,fin = interval[node][0],interval[node][1]
		if i >= deb and i <= fin:
			tmp.append(node)
	tmp.sort()
	vector = "-".join(tmp)
	if vector in scan:
		pass
	else:
		nb +=1
		scan[vector] = nb
	i+=MIN/2

#for vector in scan:
#	print vector," ",scan[vector]


process={}
for vector in scan:
	nb = scan[vector]
	tmp = vector.split("-")
	DEBUT,FIN=[],[]
	for node in tmp:
		deb,fin = interval[node][0],interval[node][1]
		DEBUT.append(deb)
		FIN.append(fin)
	process[nb]=[tmp,max(DEBUT),min(FIN)]


sequence={}
sequence["root"] = seq 

for node in roots:
	parent[node] = "root"

NB = max(process.keys())



########  Simulate evolution with recombination

print("Simulate evolution with recombination AND SELECTION")

total_m,total_r = 0,0
NU=[]
LISTE_DELTA=[]
nb=1
while nb <= NB:
	LIST= process[nb][0]
	for node in LIST:
		if node in sequence:
			pass
		else:
			sequence[node] = sequence[parent[node]]
	deb,fin = process[nb][1],process[nb][2]
	segment = fin - deb
	MEGA=[]
	dico={}
	for node in LIST:
		mutations = numpy.random.poisson(L * segment)
		i=1
		while i <= mutations:
			resu = node + "_m"
			MEGA.append(resu)
			i+=1
		recomb = numpy.random.poisson(L * segment * COEFF)
		i=1
		while i <= recomb:
			resu = node + "_r"
			MEGA.append(resu)
			i+=1
		dico[node] = [mutations,recomb]
	random.shuffle(MEGA)
	for resu in MEGA:
		a = resu.split("_")
		node = a[0]
		event = a[1]
		if event == "m":
			sequence[node] = single(sequence[node])
			total_m+=1
		elif event == "r":
			if len(LIST)>1:
				copy = list(LIST)
				copy.remove(node)
				donor = random.choice(copy)
				delta = numpy.random.exponential(DELTA)
				start = random.choice(range(L))
				end = int(round(start + delta,0))
				former,newer = sequence[node][start:end], sequence[donor][start:end]
				if len(former)>0:
					total_r +=1
					sequence[node] =   sequence[node][:start] + sequence[donor][start:end] + sequence[node][end:]
					iteration = 0
					compteur = 0
					while iteration < len(former):
						Nalpha,Nbeta = former[iteration] , newer[iteration]
						if Nalpha != Nbeta:
							compteur +=1
						iteration+=1
					nu = float(compteur) / len(former)
					NU.append(nu)
					LISTE_DELTA.append(delta)
					
	#print nb," ",process[nb]," ",mutations," ",recomb
	nb+=1


print("total_m=",total_m,"total_r=",total_r)



h=open(path + "nu.txt","w")
if len(NU)> 0:
	h.write( str(sum(NU) / len(NU)) + "\n")
else:
	h.write("NA\n")
h.close()


rename={}
f=open(path + "names.txt","r")
for l in f:
	a=l.strip("\n").split("\t")
	rename[a[0]]=a[1]

f.close()


h=open(path + "genomes.fa","w")
for node in sequence:
	if node != "root" and type[node] == "tip":
		h.write(">" + rename[node] + "\n" + sequence[node] + "\n")

h.close()




h=open(path + "rm.txt","w")
if COEFF == 0:
	h.write("r/m= 0.0\n")
	print("Effective recombination rate r/m= 0.0")
else:
	h.write("r/m= " + str(COEFF * DELTA *  sum(NU) / len(NU)) + "\n")
	print("Effective recombination rate r/m= ", COEFF * DELTA *  sum(NU) / len(NU))
h.close()


h=open(path + "detail.txt","w")
if COEFF == 0:
	rm="0.0"
	nu = "NA"
else:
	rm=str(COEFF * DELTA *  sum(NU) / len(NU))
	nu = str(sum(NU) / len(NU)) 
h.write("rm= " + rm + "\n" + "Rho= " + str(COEFF) + "\nDelta= " + str(DELTA) + "\nNu= " + str(nu) + "\nTotal_mutations= " + str(total_m) + "\nTotal_recombination_events= " + str(total_r) + "\n" )
h.close()

print("DONE")











