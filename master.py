


import os
import sys


version = sys.version


snake= str(version[0])

if snake=="2":
	print("Launching pipeline with Python 2")
elif snake=="3":
	print("Launching pipeline with Python 3")
else:
	print("Don't know this version of Python. Let's try with 3...")

print("Snake=",snake)


control_file = sys.argv[-1]

for stuff in sys.argv:
	if "master.py" in stuff:
		loc = stuff.split("master.py")[0]

print("Usage:  python master.py   control.txt")

print(loc)

# Default settings:
parameters={}
parameters["KAPPA"]=1
parameters["OUTPUT"]="output"
parameters["GC"]=50
parameters["LENGTH"]=10000
parameters["DELTA"]=100
parameters["RHO"]=0
parameters["CODONS"]="0.33,0.33,0.33"
parameters["RESCALE"]=1



f=open(control_file,"r")
for l in f:
	l=l.strip("\n").strip("\r").strip(" ")
	if "=" in l:
		a=l.split("=")
		attribute=a[0]
		if attribute[0] != "#":
			value = a[1].split("#")[0].strip(" ").strip("\t").strip(" ")
			try:
				parameters[attribute]=value
			except KeyError:
				print("Attribute ",attribute ," is not a valid option, ignoring...")

f.close()


#print(parameters)

print("\n### Welcome! ###\n")



path = parameters["OUTPUT"]
if path[-1] != "/":
	path += "/"
L = int(parameters["LENGTH"])
GC = float(parameters["GC"])
coeff = float(parameters["RESCALE"])
kappa = float(parameters["KAPPA"])
DELTA = int(parameters["DELTA"])
COEFF = float(parameters["RHO"])
sub = parameters["CODONS"]
while " " in sub:
	sub=sub.replace(" ","")

if "TREE" in parameters:
	TREE = parameters["TREE"]
else:
	print("A tree is needed for the simulations. Please come back with a tree. Exiting...")
	exit()



if snake=="2":
	try:
		os.mkdir(path)
	except OSError:
		print("Output folder already exists. Previous files will be lost.")
else:
	try:
		os.mkdir(path)
	except FileExistsError:
		print("Output folder already exists. Previous files will be lost.")



print("1. Reading the tree")
os.system("python " +  loc + "extract_names.py " + path + " " + TREE)



print("2. Extracting branch lengths and topology")

os.system("python " + loc + "branch_length.py " + path  + " renamed.tree"  )




print("3. Simulating")

print("\n################")
print("PARAMETERS: ")
print("OUTPUT= ",path)
print("Alignment Length= ",L,"bp")
print("GC%= ",GC,"% (default=50%)")
print("Branch rescaling= ",coeff," (default=1, no rescaling)")
print("Transition/Transversion bias, Kappa= ",kappa," (default=1, no bias)")
print("CODONS Frequency= ",sub,"")
print("################\n")


os.system("python " + loc + "simulation.py " + sub + " " + str(kappa) + " " + str(GC) + " " + str(L) +  " " + str(COEFF) + " " + str(DELTA) + " "  + str(coeff) + " " + path)















