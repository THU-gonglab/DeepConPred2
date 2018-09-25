#!/usr/bin/env python

import sys, os
import string
import numpy as np

def AA_index(seq_len): # classify relative position of each aa in its sequence 
	ref_loc = []
	for k in range(seq_len):
		loc = (1.0*(k+1)/seq_len)
		if (loc > 0) and (loc <= 0.2):
			ref_loc.append(0)
		elif (loc > 0.2) and (loc <= 0.4):
			ref_loc.append(1)
		elif (loc > 0.4) and (loc <= 0.6):
			ref_loc.append(2)
		elif (loc > 0.6) and (loc <= 0.8):
			ref_loc.append(3)
		elif (loc > 0.8) and (loc <= 1):
			ref_loc.append(4)

	return ref_loc


def calc_avg_AA_dist(res_element): # res_element is a sequence string   ???? average distance is for ????  Probability ????

	AA_LIST = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
	AA_sum_arr = np.zeros((20,), dtype=np.int64)

	for res in res_element:
		AA_sum_arr[AA_LIST.index(res)] += 1  # counting occurence number of each aa

	sum_value = len(res_element) # length of sequence 

	return 1.0*AA_sum_arr/sum_value






def features(identifier):

    #read SSE of examples
	SSE_info_list = []
	SSE_infile = open('../data/'+identifier+'.SSE', 'r')
	contents = SSE_infile.readlines()
	SSE_infile.close()
	N = len(contents[1].strip('\n')) # length of peptide
	for ele in contents[3:]:
		SSE_info_list.append(ele.strip('\n'))  #SSE: second structure entries, restore infos


 #    #read sequance of the same example
	# ifp1 = open('../inter_data/'+identifier+'.MSA_hmmer_fa', 'r')
	# line = ifp1.readline()
	# MSA_line = ifp1.readline().strip()  # strip() without parameters, remove space character
 # 	ifp1.close()

 #    #organize MSA sequence infos
	# tmp_line = MSA_line.replace('.', '')
	# tmp_line = tmp_line.replace('-', '')
	# index_dict = {} # holding sequence info of uppers
	# _NA_list = [] # holding sequence info of "-"s
	# NA_list = [] # holding sequence info of lowers
	# index = 0
	# _index = 0
	# for C in MSA_line:
	# 	if C.isupper():
	# 		_index += 1
	# 		index += 1
	# 		index_dict[_index] = index
	# 	elif C == '-':
	# 		_index += 1
	# 		_NA_list.append(_index)
	# 	elif C.islower():
	# 		index += 1
	# 		NA_list.append(index)
	# 	elif C == '.':
	# 		pass # removing "."s


 #    #restore plmDCA of the same example into a matrix
	# ifp2 = open('../inter_data/'+identifier+'.plmDCA', 'r')
	# plmDCA_arr = np.zeros((N, N))  # N=length of peptide ,create a zero matrix
	# for line in ifp2:
	# 	line_list = line.strip().split(',')
	# 	i1 = string.atoi(line_list[0])
	# 	i2 = string.atoi(line_list[1])
	# 	if (i1 not in _NA_list) and (i2 not in _NA_list):  #i1 and i2 both not "-"
	# 		plmDCA_arr[index_dict[i1]-1, index_dict[i2]-1] = string.atof(line_list[2])  
	# ifp2.close()
	ccmPred = np.loadtxt('../data/'+identifier+'.ccmpred')
	ccmPred = (ccmPred+ccmPred.T)/2.0 


	length_list = []
	flag_arr = np.array([], dtype=np.int64) #feature 9
	odd_even_AA_avg_dist_arr = np.array([]) #featrue 10


	for i in range(len(SSE_info_list)):  # for each sse in this example
		SSE_info = SSE_info_list[i]
		res_element = SSE_info.split()[2]  # Notion: strip and split
		
		length_list.append(len(res_element)) #res_element :  seq_ID(2,3,4,...)

		if i == 0:
			flag_arr = np.array([1, 0, 0, 0], dtype=np.int64)  #first entry of SS
		elif i == 1:
			flag_arr = np.vstack( (flag_arr, np.array([0, 1, 0, 0], dtype=np.int64)) )  #v for vertical, hstack: h for horizontal ,second entry of SS
		elif i == len(SSE_info_list)-2:
			flag_arr = np.vstack( (flag_arr, np.array([0, 0, 1, 0], dtype=np.int64)) )  # the last but one entry of SS
		elif i == len(SSE_info_list)-1:
			flag_arr = np.vstack( (flag_arr, np.array([0, 0, 0, 1], dtype=np.int64)) )  # the last entry of SS
		else:
			flag_arr = np.vstack( (flag_arr, np.array([0, 0, 0, 0], dtype=np.int64)) )  # any othors

		odd_res_element = ''
		even_res_element = ''

		for j in range(len(res_element)):  
			if j%2 == 0:
				odd_res_element += res_element[j]
			else:
				even_res_element += res_element[j]

		tmp_arr = np.hstack( (calc_avg_AA_dist(odd_res_element), calc_avg_AA_dist(even_res_element)) )   # 40 entries; ???? aim ????
		odd_even_AA_avg_dist_arr = np.vstack( (odd_even_AA_avg_dist_arr, tmp_arr) ) if odd_even_AA_avg_dist_arr.size else tmp_arr
        #odd_even_AA_avg_dist_arr contains all infos of all SSEs
		

      #feature recording
	outfile_res_num = open('../inter_data/'+identifier+'.res_num', 'w')
	outfile_flags = open('../inter_data/'+identifier+'.flags', 'w')
	outfile_even_odd_comp = open('../inter_data/'+identifier+'.even_odd_comp', 'w')
	outfile_SSE_num = open('../inter_data/'+identifier+'.SSE_num', 'w')
	outfile_coevo_info = open('../inter_data/'+identifier+'.coevo_info', 'w')
	outfile_intervening_len = open('../inter_data/'+identifier+'.intervening_len', 'w')


      # from the very first neighbouring SSEs to the last pair
	for i in range(len(SSE_info_list)-1):
		for j in range(i+1, len(SSE_info_list)): 
			SSE1_info = SSE_info_list[i]
			SSE2_info = SSE_info_list[j]

			tmp_list1 = SSE1_info.split()
			tmp_list2 = SSE2_info.split()
            

            #change char to int
			seq_ID1 = tmp_list1[1]
			seq_ID1_list = seq_ID1[7:len(seq_ID1)-2].split(',')
			for m in range(len(seq_ID1_list)):
				seq_ID1_list[m] = string.atoi(seq_ID1_list[m])

			seq_ID2 = tmp_list2[1]
			seq_ID2_list = seq_ID2[7:len(seq_ID2)-2].split(',')
			for m in range(len(seq_ID2_list)):
				seq_ID2_list[m] = string.atoi(seq_ID2_list[m])

			
			pos1 = seq_ID1_list[len(seq_ID1_list)-1]
			pos2 = seq_ID2_list[0]
			sep = pos2-pos1-1  # sep -> seperating area, how many aas in this area
			if sep <= 1:  # =1 means there is only one aa between two neighbouring SSEs
				continue


			# write names of SSEs (removing ":")
			SSE1 = tmp_list1[0]
			SSE2 = tmp_list2[0]
			outfile_res_num.write(SSE1[:-1]+' '+SSE2[:-1]+' ')
			outfile_flags.write(SSE1[:-1]+' '+SSE2[:-1]+' ')
			outfile_even_odd_comp.write(SSE1[:-1]+' '+SSE2[:-1]+' ')
			outfile_SSE_num.write(SSE1[:-1]+' '+SSE2[:-1]+' ')
			outfile_coevo_info.write(SSE1[:-1]+' '+SSE2[:-1]+' ')
			outfile_intervening_len.write(SSE1[:-1]+' '+SSE2[:-1]+' ')







			#feature 8: res_num : 6 elements: the first 'for': length of SSEs having relative position of -1, 0, 1 of current SSE i,
			#  the second 'for' , the same as the first but point to the orther SSE in the pair 
			for m in range(-1, 2):
				if i+m >= 0:
					outfile_res_num.write('%d ' % length_list[i+m]) # length of SSEs of example in length_list
				else:
					outfile_res_num.write('%d ' % 0)
			for m in range(-1, 2):
				if j+m <= len(SSE_info_list)-1:
					outfile_res_num.write('%d ' % length_list[j+m])
				else:
					outfile_res_num.write('%d ' % 0)




			#feature9 : flags: 6 elements 
			for n in range(flag_arr.shape[1]-1): # traverse columns of flag matrix except the last column (for the first in the neighbourhood, this column would never shows one)
				outfile_flags.write('%d ' % flag_arr[i, n])
			for n in range(1, flag_arr.shape[1]): # removing the first one, for the latter in the neighbourhood, this column would never shows one
				outfile_flags.write('%d ' % flag_arr[j, n])


		

			#features10 : even_odd_comp , 40 entries for each row of odd_even_AA_avg_dist_arr (20aa for odd and 20aa for even),2 SSEs, 80 entries totally
			for n in range(odd_even_AA_avg_dist_arr.shape[1]): # traverse columns
				outfile_even_odd_comp.write('%f ' % odd_even_AA_avg_dist_arr[i, n])
			for n in range(odd_even_AA_avg_dist_arr.shape[1]):
				outfile_even_odd_comp.write('%f ' % odd_even_AA_avg_dist_arr[j, n])


		

			#feature 6: SSE_num
			if SSE1[0] == 'H':
				SSE1_id = string.atoi(SSE1[5:len(SSE1)-1]) # extract 10 from "Helix10:"
			elif SSE1[0] == 'S':
				SSE1_id = string.atoi(SSE1[6:len(SSE1)-1])
			if SSE2[0] == 'H':
				SSE2_id = string.atoi(SSE2[5:len(SSE2)-1])
			if SSE2[0] == 'S':  #  ???? elif ????
				SSE2_id = string.atoi(SSE2[6:len(SSE2)-1])
			
			num = SSE2_id-SSE1_id-1   # number of interval SSEs
            
            #feature 6: uisng intervels to discretize number of interval SSEs
			if num == 0:
				outfile_SSE_num.write('1 0 0 0 0 0 0 0 ')
			elif num == 1:
				outfile_SSE_num.write('0 1 0 0 0 0 0 0 ')
			elif (num >= 2) and (num <= 5):
				outfile_SSE_num.write('0 0 1 0 0 0 0 0 ')
			elif (num >= 6) and (num <= 13):
				outfile_SSE_num.write('0 0 0 1 0 0 0 0 ')
			elif (num >= 14) and (num <= 22):
				outfile_SSE_num.write('0 0 0 0 1 0 0 0 ')
			elif (num >= 23) and (num <= 31):
				outfile_SSE_num.write('0 0 0 0 0 1 0 0 ')
			elif (num >= 32) and (num <= 40):
				outfile_SSE_num.write('0 0 0 0 0 0 1 0 ')
			else:
				outfile_SSE_num.write('0 0 0 0 0 0 0 1 ')




			#feature 1: coevo_info #plmDCA  score
			ref_loc1 = AA_index(len(seq_ID1_list))  # seq_ID1_list contains aa indexs of SSE sequence 
			ref_loc2 = AA_index(len(seq_ID2_list))

			loc_dict = {}
			for k in range(25):
				loc_dict[k] = []
			for m in range(len(seq_ID1_list)):
				for n in range(len(seq_ID2_list)):
					loc_dict[ref_loc1[m]*5+ref_loc2[n]].append(ccmPred[seq_ID1_list[m], seq_ID2_list[n]]) # plmDCA_arr is a upper trangular matrix

			for k in range(25):
				if len(loc_dict[k]) != 0:
					outfile_coevo_info.write('%f ' % np.amax(np.array(loc_dict[k])))
				else:
					outfile_coevo_info.write('0.0000 ')


			#feature4 : intervening_len
			if (sep >= 2) and (sep <= 5):
				outfile_intervening_len.write('1 0 0 0 0 0 0 0 ')
			elif (sep >= 6) and (sep <= 9):
				outfile_intervening_len.write('0 1 0 0 0 0 0 0 ')
			elif (sep >= 10) and (sep <= 12):
				outfile_intervening_len.write('0 0 1 0 0 0 0 0 ')
			elif (sep >= 13) and (sep <= 15):
				outfile_intervening_len.write('0 0 0 1 0 0 0 0 ')
			elif (sep >= 16) and (sep <= 23):
				outfile_intervening_len.write('0 0 0 0 1 0 0 0 ')
			elif (sep >= 24) and (sep <= 27):
				outfile_intervening_len.write('0 0 0 0 0 1 0 0 ')
			elif (sep >= 28) and (sep <= 65):
				outfile_intervening_len.write('0 0 0 0 0 0 1 0 ')
			else:
				outfile_intervening_len.write('0 0 0 0 0 0 0 1 ')


			outfile_res_num.write('\n')
			outfile_flags.write('\n')
			outfile_even_odd_comp.write('\n')
			outfile_SSE_num.write('\n')
			outfile_coevo_info.write('\n')
			outfile_intervening_len.write('\n')


	outfile_res_num.close()
	outfile_flags.close()
	outfile_even_odd_comp.close()
	outfile_SSE_num.close()
	outfile_coevo_info.close()
	outfile_intervening_len.close()
