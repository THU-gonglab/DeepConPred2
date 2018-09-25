import sys

pro_name=sys.argv[1]



fastafile=open('../data/'+pro_name,'r')
protein=fastafile.readlines()
fastafile.close()

length=len(protein[1].split()[0])

SSfile=open('../data/'+pro_name+'.spd33','r')
residues=SSfile.readlines()[1:]
SSfile.close()

outfile=open('../data/'+pro_name+'.SSE','w')
outfile.write(protein[0])
outfile.write(protein[1])

SS=[]

if length==len(residues):
	
	for res in range(length):
		
		if residues[res].split()[1]==protein[1][res]:
			outfile.write(residues[res].split()[2])
			SS.append(residues[res].split()[2])

	outfile.write('\n')

	Global_count=0
	current=SS[0]
	count=1
	fit_resid=[protein[1][0]]
	fit_index=[0]
	for ele in range(1,length):
		if current==SS[ele]:
			count=count+1
			fit_index.append(ele)
			fit_resid.append(protein[1][ele])
			current=SS[ele]
		else:
			if current == 'E' and count>3:
				Global_count=Global_count+1
				outfile.write('Strand{}: seq_ID('.format(str(Global_count)))
				for index in fit_index:
					outfile.write(str(index)+',')
				outfile.write(') ')
				for resid in fit_resid:
					outfile.write(resid)
				outfile.write('\n')
			elif current == 'H' and count>6:
				Global_count=Global_count+1
				outfile.write('Helix{}: seq_ID('.format(str(Global_count)))
				for index in fit_index:
					outfile.write(str(index)+',')
				outfile.write(') ')
				for resid in fit_resid:
					outfile.write(resid)
				outfile.write('\n')
			count=1
			fit_resid=[protein[1][ele]]
			fit_index=[ele]
			current=SS[ele]
	outfile.close()




