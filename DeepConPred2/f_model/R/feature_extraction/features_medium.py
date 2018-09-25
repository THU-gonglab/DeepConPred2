#!/usr/bin/env python

import sys, os
import string
import numpy as np
from funs_in_features_medium import *
# import time
import global_list_medium




def features(identifier):
	global_list_medium.init()
	# print global_list.outfile_res_pairs_SS,'sadasfsdadfadsfasdf'

	w, cw, sw = 9, 3, 5
	# sw: sliding window for smoothing PSSM
    # w: window for res pairs
    # cw: window for quantile of intervening sequence between res pairs  
	AA2id = {'A':0, 'R':1, 'N':2, 'D':3, 'C':4, 'E':5, 'Q':6, 'G':7, 'H':8, 'I':9, 'L':10, 'K':11, 'M':12, 'F':13, 'P':14, 'S':15, 'T':16, 'W':17, 'Y':18, 'V':19}
	

	SS2id = {'H2H':0, 'H2E':1, 'H2C':2, 'E2H':3, 'E2E':4, 'E2C':5, 'C2H':6, 'C2E':7, 'C2C':8}


	infile_SS = open('../data/'+identifier +'.SSE', 'r')
	line = infile_SS.readline()
	prot = infile_SS.readline().strip('\n')  # first struct (sequence)
	SS_line = infile_SS.readline().strip('\n') # which SS aa belongs
	infile_SS.close()

	infile_ACC = open('../data/'+identifier +'.ACC', 'r')
	
	line = infile_ACC.readline()
	ACC_line = infile_ACC.readline().strip() 
	infile_ACC.close()

	N = len(prot)

	smoothed_PSSM_arr = calc_smoothed_PSSM_arr(identifier, prot, sw)
	ccmPred = calc_ccmPred(identifier)
	AA_SSCP_arr = calc_AA_SSCP_arr(prot, AA2id, SS2id)
	index_seqID, index_prob = calc_index_info(identifier)

	ACC2id = {'0':1, '1':1, '2':0, '3':0, '4':0, '5':0, '6':0, '7':0, '8':0, '9':0}
	ACCNum = np.zeros((2,), dtype=np.int64)
	for ele in ACC_line:
		ACCNum[ACC2id[ele]] += 1
	ACCPct = 1.0*ACCNum/len(ACC_line) 

	
	for id1 in range(N-13):   
		tail= N if id1+25>N else id1+25
		for id2 in range(id1+13, tail):

			# print ('time: {} . \n'.format(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time()))))
			temp=''
			AA1, AA2 = prot[id1], prot[id2]
			res_pair = AA1+str(id1)+' '+AA2+str(id2)


			# print ('feature 6 time: {} . \n'.format(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time()))))
			# feature 6: residue pairs: secondary structure
			# feature dimension : 2 * 9 * 3 = 54
			temp=temp+res_pair+' '
			for i in range(id1-(w-1)/2, id1+(w-1)/2+1):
				if (i >= 0) and (i < N):
					if SS_line[i] == 'H':
						temp=temp+'1 0 0 '
					elif SS_line[i] == 'E':
						temp=temp+'0 1 0 '
					elif SS_line[i] == 'C':
						temp=temp+'0 0 1 '
				else:
					temp=temp+'0 0 0 '

			for i in range(id2-(w-1)/2, id2+(w-1)/2+1):
				if (i >= 0) and (i < N):
					if SS_line[i] == 'H':
						temp=temp+'1 0 0 '
					elif SS_line[i] == 'E':
						temp=temp+'0 1 0 '
					elif SS_line[i] == 'C':
						temp=temp+'0 0 1 '
				else:
					temp=temp+'0 0 0 '

				# print len(temp.split())
			global_list_medium.outfile_res_pairs_SS.extend([temp])

			# print ('feature 7 time: {} . \n'.format(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time()))))

			# feature 7 : residue pairs: solvent accessibility
			# feature dimension : 2 * 9 * 2 = 36
			temp=''
			temp=temp+res_pair+' '
			for i in range(id1-(w-1)/2, id1+(w-1)/2+1):
				if (i >= 0) and (i < N):
					if ACC_line[i] in '23456789':
						temp=temp+'1 0 '
					elif ACC_line[i] in '01':
						temp=temp+'0 1 '
				else:
					temp=temp+'0 0 '

			for i in range(id2-(w-1)/2, id2+(w-1)/2+1):
				if (i >= 0) and (i < N):
					if ACC_line[i] in '23456789':
						temp=temp+'1 0 '
					elif ACC_line[i] in '01':
						temp=temp+'0 1 '
				else:
					temp=temp+'0 0 '
			global_list_medium.outfile_res_pairs_ACC.extend([temp])

			# print ('feature 2 time: {} . \n'.format(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time()))))

			# feature 2 : residue pairs: smoothed PSSM
			# feature dimension : 2 * 5 * 20 = 200
			temp=''
			temp=temp+res_pair+' '
			for i in range(id1-(sw-1)/2, id1+(sw-1)/2+1):
				if (i >= 0) and (i < N):
					for j in range(20):
						temp=temp+('%d ' % smoothed_PSSM_arr[i, j])
				else:
					for j in range(20):
						temp=temp+'0 '

			for i in range(id2-(sw-1)/2, id2+(sw-1)/2+1):
				if (i >= 0) and (i < N):
					for j in range(20):
						temp=temp+('%d ' % smoothed_PSSM_arr[i, j])
				else:
					for j in range(20):
						temp=temp+'0 '
			global_list_medium.outfile_res_pairs_smoothed_PSSM.extend([temp])

			# print ('feature 5 time: {} . \n'.format(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time()))))

			# feature 5: residue pairs: coevolutionary information
			# feature dimension : 81
			temp=''
			temp=temp+res_pair+' '
			for i in range(id1-(9-1)/2, id1+(9-1)/2+1):
				for j in range(id2-(9-1)/2, id2+(9-1)/2+1):
					if (i < 0) or (i >= N):
						temp=temp+'0.0 '
					elif (j < 0) or (j >= N):
						temp=temp+'0.0 '
					else:
						temp=temp+('%f ' % ccmPred[id1, id2])
			global_list_medium.outfile_res_pairs_coevo_info.extend([temp])

			# print ('feature 4 time: {} . \n'.format(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time()))))

			# feature 4: residue pairs: contact propensity
			# feature dimension : 1
			temp=''
			temp=temp+res_pair+' '
			SS_info = SS_line[id1]+'2'+SS_line[id2]
			temp=temp+('%f ' % AA_SSCP_arr[SS2id[SS_info], AA2id[AA1], AA2id[AA2]])
			global_list_medium.outfile_res_pairs_SSCP.extend([temp])

			# print ('feature 1 time: {} . \n'.format(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time()))))

	
			#feature 1 :residue pairs: coarse contact information
			# feature dimension : 3
			temp=''
			temp=temp+res_pair+' '
			label = False
			for index in range(len(index_seqID)): #index_seqID stores aa indexes of two SSEs
				
				####Notice : index of aa in my SSE file is from 0, different from Xiong's which starts from 1
				if (id1 in index_seqID[index]) and (id2 in index_seqID[index]):
					temp=temp+('%s %s %s ' % (index_prob[index][0], index_prob[index][1], index_prob[index][2]))
					global_list_medium.outfile_res_pairs_SSConProb.extend([temp])
					label = True
					break
			if not label:
				temp=temp+'0 0 0 '
				global_list_medium.outfile_res_pairs_SSConProb.extend([temp])

			# print ('feature 6 time: {} . \n'.format(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time()))))

			# fearture 6: intervening sequences: secondary structure
			# feature dimension : 3*3*3 = 27
			temp=''
			temp=temp+(res_pair+' ')
			sep = id2-id1-1 # number of aa in intervening sequence
			sep_ = sep/2
			
			        # compute indexes of quantiles 
			if sep%2 == 0:
				mid = id1+sep_
				if sep_%2 == 0:
					lmid = id1+sep_/2
					rmid = mid+sep_/2
				else:
					lmid = id1+sep_/2+1
					rmid = mid+sep_/2+1
			else:
				mid = id1+sep_+1
				if sep_%2 == 0:
					lmid = id1+sep_/2+1
					rmid = mid+sep_/2
				else:
					lmid = id1+sep_/2+1
					rmid = mid+sep_/2+1

			for i in range(lmid-(cw-1)/2, lmid+(cw-1)/2+1):
				if SS_line[i] == 'H':
					temp=temp+('1 0 0 ')
				elif SS_line[i] == 'E':
					temp=temp+('0 1 0 ')
				elif SS_line[i] == 'C':
					temp=temp+('0 0 1 ')
			for i in range(mid-(cw-1)/2, mid+(cw-1)/2+1):
				if SS_line[i] == 'H':
					temp=temp+('1 0 0 ')
				elif SS_line[i] == 'E':
					temp=temp+('0 1 0 ')
				elif SS_line[i] == 'C':
					temp=temp+('0 0 1 ')
			for i in range(rmid-(cw-1)/2, rmid+(cw-1)/2+1):
				if SS_line[i] == 'H':
					temp=temp+('1 0 0 ')
				elif SS_line[i] == 'E':
					temp=temp+('0 1 0 ')
				elif SS_line[i] == 'C':
					temp=temp+('0 0 1 ')
			global_list_medium.outfile_intervening_seq_SS.extend([temp])

			# print ('feature 9 time: {} . \n'.format(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time()))))

			#featrue 9: intervening sequences: length
			# feature dimension : 12
			temp=''
			temp=temp+(res_pair+' ')
			if sep == 12:
				temp=temp+('1 0 0 0 0 0 0 0 0 0 0 0 ')
			elif sep == 13:
				temp=temp+('0 1 0 0 0 0 0 0 0 0 0 0 ')
			elif sep == 14:
				temp=temp+('0 0 1 0 0 0 0 0 0 0 0 0 ')
			elif sep == 15:
				temp=temp+('0 0 0 1 0 0 0 0 0 0 0 0 ')
			elif sep == 16:
				temp=temp+('0 0 0 0 1 0 0 0 0 0 0 0 ')
			elif sep == 17:
				temp=temp+('0 0 0 0 0 1 0 0 0 0 0 0 ')
			elif sep == 18:
				temp=temp+('0 0 0 0 0 0 1 0 0 0 0 0 ')
			elif sep == 19:
				temp=temp+('0 0 0 0 0 0 0 1 0 0 0 0 ')
			elif sep == 20:
				temp=temp+('0 0 0 0 0 0 0 0 1 0 0 0 ')
			elif sep == 21:
				temp=temp+('0 0 0 0 0 0 0 0 0 1 0 0 ')
			elif sep == 22:
				temp=temp+('0 0 0 0 0 0 0 0 0 0 1 0 ')
			elif sep == 23:
				temp=temp+('0 0 0 0 0 0 0 0 0 0 0 1 ')
			global_list_medium.outfile_intervening_seq_len.extend([temp])


			# print ('feature 3 time: {} . \n'.format(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time()))))

			# feature 3: intervening sequences: natural vector
			# feature dimension : 60
			temp=''
			temp=temp+(res_pair+' ')
			sep_frag = prot[id1+1:id2]
			nv_arr = natural_vector(sep_frag)
			for i in range(nv_arr.shape[0]):
				temp=temp+('%f ' % nv_arr[i])
			global_list_medium.outfile_intervening_seq_NV.extend([temp])


			# feature 11 : entire protein: solvent accessibility
			# feature dimension : 2
			temp=''
			temp=temp+(res_pair+' ')+('%f %f ' % (ACCPct[0], ACCPct[1]))
			global_list_medium.outfile_entire_prot_ACC.extend([temp])

			# print ('feature 13 time: {} . \n'.format(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time()))))

			#feature 13: entire protein: length
			# feature dimension : 4
			temp=''
			temp=temp+(res_pair+' ')
			if N <= 80:
				temp=temp+('1 0 0 0 ')
			elif (N > 80) and (N <= 160):
				temp=temp+('0 1 0 0 ')
			elif (N > 160) and (N <= 240):
				temp=temp+('0 0 1 0 ')
			elif N > 240:
				temp=temp+('0 0 0 1 ')
			global_list_medium.outfile_entire_prot_len.extend([temp])

	
