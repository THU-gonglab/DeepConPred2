import sys
import random

from features_short import *
import global_list_short

infile_name=sys.argv[1]
if not os.path.exists('../features/short'):
	os.mkdir('../features/short')

protein_fa_file = open('../../../data/'+infile_name , 'r')
contents = protein_fa_file.read()[1:].split('\n>')
protein_fa_file.close()
proteins = {}
for ele in contents:
	ele_list = ele.split('\n')
	name=ele_list[0].split()[0]
	proteins[name]=ele_list[1]

protein_name=proteins.keys()

count=0
total=float(len(protein_name))

for identifier in protein_name:
	count+=1
	
	
	features(identifier)

	# print ('end time: {} . \n'.format(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time()))))

	check=[len(global_list_short.outfile_res_pairs_SS),
	len(global_list_short.outfile_res_pairs_ACC),
	len(global_list_short.outfile_res_pairs_smoothed_PSSM),
	len(global_list_short.outfile_res_pairs_coevo_info),
	len(global_list_short.outfile_res_pairs_SSCP ),
	len(global_list_short.outfile_res_pairs_SSConProb ),
	len(global_list_short.outfile_intervening_seq_SS),
	len(global_list_short.outfile_intervening_seq_len ),
	len(global_list_short.outfile_intervening_seq_NV),
	len(global_list_short.outfile_entire_prot_ACC),
	len(global_list_short.outfile_entire_prot_len)]
	# print check
	# print max(check),' ',min(check)
	if max(check)!=min(check):
		print 'length not equal'
		sys.exit(1)


	prot=proteins[identifier]
	N=len(prot)

	# infile_true_label = open('../../true_contact/'+identifier+'_true_contact', 'r')
	# contents = infile_true_label.readlines()
	# infile_true_label.close()

	# label = {}
	# for ele in contents:
	# 	ele_list = ele.split()
	# 	label[ele_list[0]] = ele_list[1]



	
	outfile = open('../features/short/'+identifier+'_features', 'w')
	outfile_obj = open('../features/short/'+identifier+'_features_obj', 'w')

	# for i in range(2):
	for i in range(len(global_list_short.outfile_res_pairs_ACC)):

		line1 = global_list_short.outfile_res_pairs_SS[i]
		line2 = global_list_short.outfile_res_pairs_ACC[i]
		line3 = global_list_short.outfile_res_pairs_smoothed_PSSM[i]
		line4 = global_list_short.outfile_res_pairs_coevo_info[i]
		line5 = global_list_short.outfile_res_pairs_SSCP[i]
		line6 = global_list_short.outfile_res_pairs_SSConProb[i]
		line7 = global_list_short.outfile_intervening_seq_SS[i]
		line8 = global_list_short.outfile_intervening_seq_len[i]
		line9 = global_list_short.outfile_intervening_seq_NV[i]
		line10 = global_list_short.outfile_entire_prot_ACC[i]
		line11 = global_list_short.outfile_entire_prot_len[i]

		pos1 = line1.index(' ')
		res1 = line1[:pos1]
		line1 = line1[pos1+1:]
		pos2 = line1.index(' ')
		res2 = line1[:pos2]
		feature1 = line1[pos2+1:]
		pos = pos1+pos2+2

		feature2 = line2[pos:]
		feature3 = line3[pos:]
		feature4 = line4[pos:]
		feature5 = line5[pos:]
		feature6 = line6[pos:]
		feature7 = line7[pos:]
		feature8 = line8[pos:]
		feature9 = line9[pos:]
		feature10 = line10[pos:]
		feature11 = line11[pos:]

		# print len((feature1+feature2+feature3+feature4+feature5+feature6+feature7+feature8+feature9+feature10+feature11+label[res1+':'+res2]).split())
		


		# print line1
		# print feature1+feature2+feature3+feature4+feature5+feature6+feature7+feature8+feature9+feature10+feature11+label[res1+':'+res2]


		outfile.write(feature1+feature2+feature3+feature4+feature5+feature6+feature7+feature8+feature9+feature10+feature11+'\n')
		outfile_obj.write(res1+' '+res2+'\n')

	
	outfile.close()
	outfile_obj.close()
	print 'processing %.2f' %(count/total*100)
