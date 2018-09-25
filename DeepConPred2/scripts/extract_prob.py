import os 
import sys

pro_name=sys.argv[1]


with open('./spider_cannot_list','a') as logfile:
	
	if os.path.exists('../data/'+pro_name+'.spd33'):
		spd33=open('../data/'+pro_name+'.spd33','r')
		contents=spd33.readlines()
		spd33.close()
		if contents:
			contents=contents[1:]
			probs=[]
			for line in contents:
				line=line.split()
				probs.append([i for i in line[10:]])
			with open('../f_model/New_Third/data/spd_prob/'+pro_name+'.spd_prob','w') as outfile:
				for line in probs:
					for a in line:
						outfile.write(a+' ') 
					outfile.write('\n')
		else:
			print 'emptyfile:{}'.format(pro_name)
			
			logfile.write(pro_name+'\n')
	else:
		print 'nofile:{}'.format(pro_name)
		
		logfile.write(pro_name+'\n')



