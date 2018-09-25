import sys, os
import string
import numpy as np

def calc_smoothed_PSSM_arr(identifier, prot,sw):
	N = len(prot)
	PSSM_arr = np.zeros((N, 20), dtype=np.int64)
	smoothed_PSSM_arr = np.zeros((N, 20), dtype=np.int64)

	infile_PSSM = open('../data/'+identifier+'.PSSM', 'r')
	for i in range(3):
		line = infile_PSSM.readline()
	if line[11:69] != 'A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V':
		print 'PSSM file format error'
		sys.exit(1)
	seq = ''
	i = 0
	while True:
		line = infile_PSSM.readline()
		if line == '\n':
			break

		begin_id = 9
		for j in range(20):
			if (line[begin_id] != ' ') and (line[begin_id] != '-'):
				print 'PSSM file format exception'
				sys.exit(1)
			PSSM_arr[i, j] = string.atoi(line[begin_id:begin_id+3]) # ignore the last two columns
			begin_id += 3
		i += 1
		seq += line[6] # seq -> sequence 
	if (seq != prot) or (i != N):
		print 'The different protein sequences!\nprotein:\n%s\nPSSM_sequence:\n%s' % (prot, seq)
		print '%d != %d' % (i, N)
		sys.exit(1)
	infile_PSSM.close()

	for i in range(N):
		for j in range(i-(sw-1)/2, i+(sw-1)/2+1):
			if (j >= 0) and (j <= N-1):
				smoothed_PSSM_arr[i, :] += PSSM_arr[j, :] # sw: sliding window for smoothing PSSM

	return smoothed_PSSM_arr # single PSSM value turns into sum of 5 vicinities.		



def calc_ccmPred(identifier): 
	ccmPred = np.loadtxt('../data/'+identifier+'.ccmpred')
	ccmPred = (ccmPred+ccmPred.T)/2.0 
	return ccmPred


def calc_AA_SSCP_arr(prot,AA2id, SS2id): 
	AA_SSCP_arr = np.zeros((9, 20, 20)) 
	infile_AA_SSCP = open('../data/AA_SSCP_db', 'r')  # 9 (20*20) probability matrices
	line = infile_AA_SSCP.readline()
	while True:
		line = infile_AA_SSCP.readline()
		if not line:
			break

		if line[:20] == 'Secondary structure:':
			SS_info = line[21:24]
			line = infile_AA_SSCP.readline()
			col_AA = line.split()
			for i in range(20):
				line = infile_AA_SSCP.readline()
				row_AA = line[0]
				if (AA2id[row_AA] != i) or (row_AA != col_AA[i]):
					print 'AA column orders error in SSCP database file'
					sys.exit(1)
				line = line[1:]
				line_list = line.split()
				for j in range(len(line_list)):
					AA_SSCP_arr[SS2id[SS_info], AA2id[row_AA], AA2id[col_AA[j]]] = string.atof(line_list[j])
	infile_AA_SSCP.close()

	return AA_SSCP_arr

# SSE pair contact probability, return two dictionary 
def calc_index_info(identifier):
	index_seqID, index_prob = {}, {}
	infile_result_SS = open('../data/'+identifier+'.SSE', 'r')
	contents = infile_result_SS.readlines()[3:]
	infile_result_SS.close()

	if contents:
		index = 0

		for i in range(len(contents)-1):
			for j in range(i+1, len(contents)):

				SSE1_info = contents[i]
				SSE2_info = contents[j]

				seq_ID1 = SSE1_info.split()[1]
				seq_ID1_list = seq_ID1[7:len(seq_ID1)-2].split(',')
				for k in range(len(seq_ID1_list)):
					seq_ID1_list[k] = string.atoi(seq_ID1_list[k])
				seq_ID2 = SSE2_info.split()[1]
				seq_ID2_list = seq_ID2[7:len(seq_ID2)-2].split(',')
				for k in range(len(seq_ID2_list)):
					seq_ID2_list[k] = string.atoi(seq_ID2_list[k])
				

				#### This has been checked in feature extraction in DeepCCon, not necessary here and may cause conflict with below procedure
				pos1 = seq_ID1_list[len(seq_ID1_list)-1]
				pos2 = seq_ID2_list[0]
				# ???? these 2 variables does not used ????

				sep = pos2-pos1-1
				if sep <= 1:
					continue

				index_seqID[index] = seq_ID1_list+seq_ID2_list # sum of array, stack together

				# print index
				# print index_seqID[index]

				index += 1

	infile_result_DeepCCon = open('../data/'+identifier+'.DeepCCon', 'r')
	contents = infile_result_DeepCCon.readlines()[1:]
	infile_result_DeepCCon.close()

	if contents:

		index = 0
		for line in contents:
			# print index
			index_prob[index] = line.split()[2:]
			
			# print index_prob[index]

			index += 1

	return index_seqID, index_prob


def natural_vector(prot_seq):
	AA_LIST = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
	AA_dict = {}
	for i in range(len(AA_LIST)):
		AA_dict[i] = []
	for i in xrange(len(prot_seq)):
		AA_dict[AA_LIST.index(prot_seq[i])].append(i) # record appearing position of each aa

	n_arr = np.zeros((20,), dtype=np.int64)
	u_arr = np.zeros((20,))
	D_arr = np.zeros((20,))

	for i in range(len(AA_LIST)):
		n_arr[i] = len(AA_dict[i]) # occurance number of aa
		if n_arr[i] > 0:
			u_arr[i] = 1.0*sum(np.array(AA_dict[i]))/n_arr[i]  # expectation of aa occuring position 
		else:
			u_arr[i] = 0

	n = len(prot_seq)
	for i in range(len(AA_LIST)):
		t = 0
		for j in range(n_arr[i]):
			t += (AA_dict[i][j] - u_arr[i])**2

		if n_arr[i] > 0:
			D_arr[i] = 1.0*t/(n_arr[i]*n)  
		else:
			D_arr[i] = 0

	nv = np.zeros((60,))
	for i in range(20):
		nv[3*i] = n_arr[i]
		nv[3*i+1] = u_arr[i]
		nv[3*i+2] = D_arr[i]

	return nv