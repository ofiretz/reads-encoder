import sys
import csv


def reverse(seq):
	tmp = list(seq)
	for i in range(100):
		if seq[i] == 'a':
			tmp[99-i] = 't'
		elif seq[i] == 'c':
			tmp[99-i] = 'g'
		elif seq[i] == 'g':
			tmp[99-i] = 'c'
		elif seq[i] == 't':
			tmp[99-i] = 'a'

	return tmp

def reverseChar(c):
	if c == 'a':
		return 't'
	elif c == 'c':
		return 'g'
	elif c == 'g':
		return 'c'
	elif c == 't':
		return 'a'



FN = open(sys.argv[2], "rb")
fout = open(sys.argv[3], "w")
i=0;

line = FN.readline()

while line:
	fout.write(line)
	line = FN.readline()
	i=i+1

FN.close()
csvfile = open(sys.argv[1], "rb")
reader = csv.reader(csvfile, delimiter="\t")
for readpos, strand, readseq, varstr in reader:
	l = list(readseq)

	var = eval(varstr)
	if strand == '1':
		l = reverse(readseq)
		
	if var:
		for varpos, vartype, varmisc, quals in var:
			if strand == '1':
				l[99-int(varpos)] = varmisc[0]
			else: 
				l[int(varpos)] = varmisc[0]



	fout.write(''.join(l).upper())
	fout.write("\n")



csvfile.close()
fout.close()