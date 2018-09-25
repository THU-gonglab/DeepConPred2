import sys,os

from features import *

infile_name=sys.argv[1]
if not os.path.exists('../features'):
	os.mkdir('../features')
if not os.path.exists('../inter_data'):
	os.mkdir('../inter_data')

	
protein_fa_file = open('../../../data/'+infile_name , 'r')
#sequence information 
contents = protein_fa_file.read()[1:].split('\n>')
protein_fa_file.close()
proteins = {}
for ele in contents:
	ele_list = ele.split('\n')
	name=ele_list[0].split()[0]
	proteins[name]=ele_list[1]

protein_name=proteins.keys()
for identifier in protein_name:
	prot=proteins[identifier]
	N=len(prot)
	features(identifier)

	# infile_true_label = open('../../SS_ture_contact/'+identifier+'_SS_true_contct', 'r')
	# contents = infile_true_label.readlines()
	# infile_true_label.close()
	# ind=0
	# for i in range(len(contents)):
	# 	if contents[i][0] == '#':
	# 		ind = i   # "#..." is a boundary
	# 		break
	# contents = contents[(ind+2):]  # removing forward sequence infos and SSE infos


	# label = {}
	# for ele in contents:
	# 	ele_list = ele.split()
	# 	label[ele_list[0]] = ele_list[1]

	infile_res_num = open('../inter_data/'+identifier+'.res_num', 'r') #6 #6
	infile_flags = open('../inter_data/'+identifier+'.flags', 'r') #6 #12
	infile_even_odd_comp = open('../inter_data/'+identifier+'.even_odd_comp', 'r') #80 #92
	infile_SSE_num = open('../inter_data/'+identifier+'.SSE_num', 'r') #8 #100
	infile_coevo_info = open('../inter_data/'+identifier+'.coevo_info', 'r') #25 #125
	infile_intervening_len = open('../inter_data/'+identifier+'.intervening_len', 'r') #8 #133

	outfile = open('../features/'+identifier+'_features', 'w')
	outfile_obj = open('../features/'+identifier+'_features_obj', 'w')


	while True:
		line1 = infile_res_num.readline()
		if not line1:
			break   # read row by row, line1==NULL is stop condition

		line1 = line1.strip('\n')
		pos1 = line1.index(' ')
		SSE1, line1 = line1[:pos1], line1[(pos1+1):]
		pos2 = line1.index(' ')
		SSE2, feature_values1 = line1[:pos2], line1[(pos2+1):]
		pos = pos1+pos2+2

		line2 = infile_flags.readline().strip('\n')
		line3 = infile_even_odd_comp.readline().strip('\n')
		line4 = infile_SSE_num.readline().strip('\n')
		line5 = infile_coevo_info.readline().strip('\n')
		line6 = infile_intervening_len.readline().strip('\n')   ## !!!!All features extracted before have same format!!!!

		feature_values2 = line2[pos:]
		feature_values3 = line3[pos:]
		feature_values4 = line4[pos:]
		feature_values5 = line5[pos:]
		feature_values6 = line6[pos:]

		outfile.write(feature_values1+feature_values2+feature_values3+feature_values4+feature_values5+feature_values6+'\n')
		outfile_obj.write(SSE1+' '+SSE2+'\n')

	infile_res_num.close()
	infile_flags.close()
	infile_even_odd_comp.close()
	infile_SSE_num.close()
	infile_coevo_info.close()
	infile_intervening_len.close()
	outfile.close()
	outfile_obj.close()

