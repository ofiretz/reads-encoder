import sys

fin1 = open(sys.argv[1], "r")
line1 = fin1.readline()
list1 = []
list1_cnt = {}
list2_cnt = {}

i=1
j=1

while line1:
	list1.append(line1)
	line1 = fin1.readline()
	#list1_cnt[line1] =+ 1
	i=i+1

fin1.close()

fin2 = open(sys.argv[2], "r")

line2 = fin2.readline()
list2 = []


while line2:
	list2.append(line2)
	line2 = fin2.readline()
	#list2_cnt[line2] =+ 1
	j=j+1

fin2.close()

#sorted(list1, key=str.lower)
#sorted(list2, key=str.lower)
list1.sort()
list2.sort()
#print(list1[0:20])
#print("\n\n\n")
#print(list2[0:20])

if(i != j):
	print("Not equal number of reads!\n")

print(i)
print(j)



t=0
k=0
f=0
while k < i-1:
	#print(s)
	#print(s.lower())
	#if not list1[k] in list2:
	if list1[k] != list2[t]:
		print(len(list1[k]))
		print(len(list2[t]))
		print(list1[k])
		print(list2[t])
		print(list1[k].lower())
		#if f == 0:
		#	k=k-1
		#	f=1
		#else:
		#t=t-1
		#	f=0

	t=t+1
	k=k+1

		
	#if s not in list1:
		#t+=1
		#print("false\n")
		#print(s)
		#print(s.lower())



print(t)


