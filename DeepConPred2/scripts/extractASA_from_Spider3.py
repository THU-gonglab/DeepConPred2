import sys

MaxACC= {'A':115, 'R':225, 'N':160, 'D':150, 'C':135, 'E':190, 'Q':180, 'G':75, 'H':195, 'I':175, 'L':170, 'K':200, 'M':185, 'F':210, 'P':145, 'S':115, 'T':140, 'W':255, 'Y':230, 'V':155}

pro_name=sys.argv[1]



fastafile=open('../data/'+pro_name,'r')
protein=fastafile.readlines()
fastafile.close()

length=len(protein[1].split()[0])

SSfile=open('../data/'+pro_name+'.spd33','r')
residues=SSfile.readlines()[1:]
SSfile.close()

outfile=open('../data/'+pro_name+'.ACC','w')
outfile.write(protein[0])
# outfile.write(protein[1])




if length==len(residues):
	
	for res in range(length):
		
		if residues[res].split()[1]!=protein[1][res]:
			# outfile.write(residues[res].split()[2])
		# else:
			sys.exit(1)
			print 'sequence wrong for ACC extraction!'
			

	# outfile.write('\n')

	for res in range(length):
		ASA=float(residues[res].split()[3])
		rASA=str(int(int(10*ASA/MaxACC[residues[res].split()[1]])))
		outfile.write(rASA)
	outfile.write('\n')
	outfile.close()




