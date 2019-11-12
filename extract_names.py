
import sys


TREE = sys.argv[-1]


path = sys.argv[-2]



f=open(TREE,"r")
l=f.readline()
memo=l

l = l.replace("(",",")
l = l.replace(")",",")
l = l.replace(":",",")
l = l.replace(";",",")
a=l.strip("\n").split(",")

#print(a)

h=open(path + "names.txt","w")
nb=0
for stuff in a:
	if stuff != "":
		try:
			taxa = float(stuff)
		except ValueError:
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
			memo = memo.replace(taxa,new)

h.close()

print("Found ",nb,"taxa in the tree")


h=open(path + "renamed.tree","w")
h.write(memo)
h.close()



