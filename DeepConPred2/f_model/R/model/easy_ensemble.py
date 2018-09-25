import sys
protein_fa_file = open('../../'+sys.argv[1] , 'r')

contents = protein_fa_file.read()[1:].split('\n>')
protein_fa_file.close()
proteins = {}
for ele in contents:
	ele_list = ele.split('\n')
	name=ele_list[0].split()[0]
	proteins[name]=ele_list[1]

protein_name=proteins.keys()

for identifier in protein_name:
	file40=open('../results/long_40/'+identifier+'.DeepRCon','r')
	file50=open('../results/long_50/'+identifier+'.DeepRCon','r')
	file60=open('../results/long_60/'+identifier+'.DeepRCon','r')
	contents40=file40.readlines()
	contents50=file50.readlines()[1:]
	contents60=file60.readlines()[1:]

	if len(contents40)-1 == len(contents60) and len(contents40)-1 == len(contents50):

		outfile=open('../results/long/'+identifier+'.DeepRCon','w')
		outfile.write(contents40[0])
		contents40=contents40[1:]
		for i in range(len(contents50)):
			ele40=contents40[i].split()
			ele50=contents50[i].split()
			ele60=contents60[i].split()
			contact=(float(ele40[2])+float(ele50[2])+float(ele60[2]))/3
			noncontact=1-contact
			outfile.write(' '.join((ele40[0],ele40[1],str(contact),str(noncontact)))+'\n')
		outfile.close()
	else:
		sys.exit(1)





