
import sys


TREE = sys.argv[-1]


path = sys.argv[-2]


#TREE = "clonal.tree"
#path=""


f=open(TREE,"r")
l=f.readline()

symbols=[",","(",")",":",";"]



link={}
bag=[]
tag=0
l=l.strip("\n")
i=1
while i < len(l):
	j=i-1
	if l[i] in symbols:
		if tag == 1:
			bag.append(taxa)
			link[taxa] = souvenir + taxa
			tag=0
	else:
		if l[j] == "(" or l[j] == ",":
			souvenir=l[j]
			taxa = l[i]
			tag=1
		else:
			taxa+=l[i]
	i+=1


print(bag)


memo=l


h=open(path + "names.txt","w")
nb=0
for stuff in bag:
	if 1==1:
		if 2==2:
			taxa = stuff
			#print(taxa)
			nb+=1
			if nb < 10:
				new = "BOZO000" + str(nb) + ":"
			elif nb < 100:
				new = "BOZO00" + str(nb)  + ":"
			elif nb < 1000:
				new = "BOZO0" + str(nb)  + ":"
			elif nb < 10000:
				new = "BOZO" + str(nb)  + ":"
			h.write(new.strip(":") + "\t" + taxa + "\n")
			taxa = taxa + ":"
			new = link[taxa.strip(":")][0] + new 
			thingy = link[taxa.strip(":")] + ":"
			memo = memo.replace(thingy,new)

h.close()

print("Found ",nb,"taxa in the tree")


new=""
tmp=""
tag=0
for L in memo:
	if L == ")":
		tag=1
		tmp+=L
	if tag==0:
		new+=L
	elif tag==1:
		if L==":":
			tmp+=L
			new+=tmp
			tag=0
			tmp=""

new += tmp 
new+=";"

h=open(path + "renamed.tree","w")
h.write(new)
h.close()


