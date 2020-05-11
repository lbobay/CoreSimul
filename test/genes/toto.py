
seq={}
f=open("gene6.fa","r")
for l in f:
	if l[0]==">":
		id = l.strip("\n").strip(">")
		seq[id]=""
	else:
		seq[id] += l.strip("\n")
	

un,deux,trois=0,0,0
i=0
while i < len(seq[id]):
	tmp1,tmp2,tmp3=[],[],[]
	for id in seq:
		codon = seq[id][i:i+3]
		if codon!= "---":
			if codon[0] not in tmp1:
				tmp1.append(codon[0])
			if codon[1] not in tmp2:
				tmp2.append(codon[1])
			if codon[2] not in tmp3:
				tmp3.append(codon[2])
			un = un + len(tmp1) - 1
			deux = deux + len(tmp2) - 1
			trois = trois + len(tmp3) - 1
	i+=3
total = un+deux+trois
print(un/total,deux/total,trois/total)