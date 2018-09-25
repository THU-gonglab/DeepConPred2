import tensorflow as tf
import numpy as np
import scipy.io
import sys
import time
import os
import string

class DeepCCon(object):
	"""docstring for DeepCCon"""
	def __init__(self,model_path):
		
		args=scipy.io.loadmat(model_path+'BPEpoch96Par.mat')
		self.w1=args['w1']
		self.w2=args['w2']
		self.w3=args['w3']
		self.w_class=args['w_class']
		args=scipy.io.loadmat(model_path+'max_min_vals.mat')
		self.max_val1=args['max_val1'][0]
		self.max_val2=args['max_val2'][0]
		self.min_val1=args['min_val1'][0]
		self.min_val2=args['min_val2'][0]
	def prediction(self,protein,Dtest_path,prediction_path):
		Dtest = np.loadtxt(Dtest_path+protein+'_features')
		if len(Dtest.shape)==1 and Dtest.shape[0]==0:
			outputfile=open(prediction_path+protein+'.DeepCCon','w')
			outputfile.close()

		else:
			if Dtest.ndim==1:
				Dtest=Dtest[np.newaxis,:]
			# Dtest=self.Dtest
			Dtest[:,:6] = 1.0*(Dtest[:,:6]-self.min_val1)/(self.max_val1-self.min_val1)
			for c in range(100,125):
				Dtest[:,c] = 1.0*(Dtest[:,c]-self.min_val2[c-100])/(self.max_val2[c-100]-self.min_val2[c-100])
			N=Dtest.shape[0]
			# np.savetxt("./f_model/data.txt",Dtest,fmt='%.4f')
			Dtest=np.hstack((Dtest,np.ones((N,1))))
			X=tf.constant(Dtest,dtype=tf.float64)
			w1=tf.constant(self.w1,dtype=tf.float64)
			w2=tf.constant(self.w2,dtype=tf.float64)
			w3=tf.constant(self.w3,dtype=tf.float64)
			w_class=tf.constant(self.w_class,dtype=tf.float64)
			a1= tf.concat( [tf.nn.sigmoid(tf.matmul(X,w1)), tf.ones((N,1),dtype=tf.float64)],1)
			a2= tf.concat( [tf.nn.sigmoid(tf.matmul(a1,w2)), tf.ones((N,1),dtype=tf.float64)],1)
			a3= tf.concat( [tf.nn.sigmoid(tf.matmul(a2,w3)), tf.ones((N,1),dtype=tf.float64)],1)
			out=tf.nn.softmax(tf.matmul(a3,w_class))
			# out_up=tf.exp(tf.matmul(a3,w_class))
			# out_down=tf.reduce_sum(out_up,axis=1,keepdims=True)
			# out=out_up/tf.concat([out_down,out_down,out_down],1)

			with tf.Session() as sess:
				predicion = sess.run(out)
				inputfile=open(Dtest_path+protein+'_features_obj','r')
				contents=inputfile.readlines()
				inputfile.close()
				outputfile=open(prediction_path+protein+'.DeepCCon','w')
				outputfile.write('#SSE1      SSE2      parallel_prob      anti-parallel_prob      no-contact_prob\n')
				for i in range(len(contents)):
					content=contents[i].split()
					outputfile.write('{} {}      {:.6f}      {:.6f}      {:.6f}\n'.format(content[0],content[1],predicion[i][0],predicion[i][1],predicion[i][2]))
				outputfile.close()




class DeepRCon(object):
	"""docstring for DeepCCon"""
	def __init__(self,model_path):
		
		self.w1={}
		self.w2={}
		self.w3={}
		self.w_class={}
		self.max_val1={}
		self.max_val2={}
		self.min_val1={}
		self.min_val2={}
		for i in ['long_40','long_50','long_60','medium','short']:
			args=scipy.io.loadmat(model_path+i+'_model.mat')
			self.w1[i]=args['w1']
			self.w2[i]=args['w2']
			self.w3[i]=args['w3']
			self.w_class[i]=args['w_class']
			args=scipy.io.loadmat(model_path+'max_min_vals_'+i+'.mat')
			self.max_val1[i]=args['max_val1'][0]
			self.max_val2[i]=args['max_val2'][0]
			self.min_val1[i]=args['min_val1'][0]
			self.min_val2[i]=args['min_val2'][0]
	def prediction(self,protein,Dtest_path,prediction_path):
		n=tf.placeholder(dtype=tf.int32)
		X=tf.placeholder(dtype=tf.float64)
		w1=tf.placeholder(dtype=tf.float64)
		w2=tf.placeholder(dtype=tf.float64)
		w3=tf.placeholder(dtype=tf.float64)
		w_class=tf.placeholder(dtype=tf.float64)
		a1= tf.concat( [tf.nn.sigmoid(tf.matmul(X,w1)), tf.ones((n,1),dtype=tf.float64)],1 )
		a2= tf.concat( [tf.nn.sigmoid(tf.matmul(a1,w2)), tf.ones((n,1),dtype=tf.float64)],1 )
		a3= tf.concat( [tf.nn.sigmoid(tf.matmul(a2,w3)), tf.ones((n,1),dtype=tf.float64)],1 )
		out=tf.nn.softmax(tf.matmul(a3,w_class))
		for i in ['long','medium','short']:	
			Dtest = np.loadtxt(Dtest_path+i+'/'+protein+'_features')
			if Dtest.shape[0]!=0 and Dtest.shape[1]!=0:
				predicion={}
				if i=='long':
					for j in ['40','50','60']:
						Dtest = np.loadtxt(Dtest_path+i+'/'+protein+'_features')
						item=i+'_'+j
						for c in range(90,290):
							Dtest[:,c] = 1.0*(Dtest[:,c]-self.min_val1[item])/(self.max_val1[item]-self.min_val1[item])		
						for c in range(412,472):
							Dtest[:,c] = 1.0*(Dtest[:,c]-self.min_val2[item][c-412])/(self.max_val2[item][c-412]-self.min_val2[item][c-412])
						N=Dtest.shape[0]
						Dtest=np.hstack((Dtest,np.ones((N,1))))
						with tf.Session() as sess:
							predicion[item]= sess.run(out,feed_dict={X:Dtest,w1:self.w1[item],w2:self.w2[item],w3:self.w3[item],w_class:self.w_class[item],n:N})
					predicion[i]=(predicion['long_40']+predicion['long_50']+predicion['long_60'])/3

				else:
					item=i
					if item=='medium':
						for c in range(90,290):
							Dtest[:,c] = 1.0*(Dtest[:,c]-self.min_val1[item])/(self.max_val1[item]-self.min_val1[item])		
						for c in range(414,474):
							Dtest[:,c] = 1.0*(Dtest[:,c]-self.min_val2[item][c-414])/(self.max_val2[item][c-414]-self.min_val2[item][c-414])
					else:
						for c in range(50,250):
							Dtest[:,c] = 1.0*(Dtest[:,c]-self.min_val1[item])/(self.max_val1[item]-self.min_val1[item])		
						for c in range(350,410):
							Dtest[:,c] = 1.0*(Dtest[:,c]-self.min_val2[item][c-350])/(self.max_val2[item][c-350]-self.min_val2[item][c-350])
					N=Dtest.shape[0]
					Dtest=np.hstack((Dtest,np.ones((N,1))))
					with tf.Session() as sess:
						predicion[item]= sess.run(out,feed_dict={X:Dtest,w1:self.w1[item],w2:self.w2[item],w3:self.w3[item],w_class:self.w_class[item],n:N})

				inputfile=open(Dtest_path+i+'/'+protein+'_features_obj','r')
				contents=inputfile.readlines()
				inputfile.close()
				outputfile=open(prediction_path+i+'/'+protein+'.DeepRCon','w')
				outputfile.write('#res1 res2      contact_prob      no-contact_prob\n')
				for k in range(len(contents)):
					content=contents[k].split()
					outputfile.write('{} {}      {:.6f}      {:.6f}\n'.format(content[0],content[1],predicion[i][k][0],predicion[i][k][1]))





class new_Third(object):
	def __init__(self,model_path):
		self.model_path=model_path
	def data_process(self,protein,data_path):
		ccmPred = np.loadtxt(data_path+'ccmpred/'+protein+'.ccmpred')
		ccmPred = (ccmPred+ccmPred.T)/2.0 
		self.N=ccmPred.shape[1]

		deepRcon = np.zeros((self.N, self.N))
		short_file = open(data_path+'short/'+protein+'.DeepRCon', 'r')
		medium_file = open(data_path+'medium/'+protein+'.DeepRCon', 'r')
		long_file = open(data_path+'long/'+protein+'.DeepRCon', 'r')

		contents_short = short_file.readlines()[1:]
		contents_medium = medium_file.readlines()[1:]
		contents_long = long_file.readlines()[1:]

		short_file.close()
		medium_file.close()
		long_file.close()


		for line in contents_short:
			line_list = line.split()
			id1 = string.atoi(line_list[0][1:])
			id2 = string.atoi(line_list[1][1:])
			deepRcon[id1, id2] = string.atof(line_list[2])
			deepRcon[id2, id1] = deepRcon[id1, id2]   

		for line in contents_medium:
			line_list = line.split()
			id1 = string.atoi(line_list[0][1:])
			id2 = string.atoi(line_list[1][1:])
			deepRcon[id1, id2] = string.atof(line_list[2])
			deepRcon[id2, id1] = deepRcon[id1, id2]   

		for line in contents_long:
			line_list = line.split()
			id1 = string.atoi(line_list[0][1:])
			id2 = string.atoi(line_list[1][1:])
			deepRcon[id1, id2] = string.atof(line_list[2])
			deepRcon[id2, id1] = deepRcon[id1, id2]
	
		spider=np.loadtxt(data_path+'spd_prob/{}.spd_prob'.format(protein))
	
		length=self.N
		spider2d=np.concatenate((spider[np.newaxis,:,:]*np.ones((length,1,1)),spider[:,np.newaxis,:]*np.ones((1,length,1))),axis=2)
		pos=abs(np.arange(length)[np.newaxis]-np.arange(length)[:,np.newaxis])
		feature2d=np.concatenate((ccmPred[:,:,np.newaxis],deepRcon[:,:,np.newaxis],spider2d,pos[:,:,np.newaxis]),axis=2)
		return feature2d[np.newaxis]
	
	def prediction(self,protein,result_path,feature):

		ensemble=[]

		for model in range(5):
			model=str(model+1)

			tf.reset_default_graph()

			with tf.Session() as sess:
				
				saver = tf.train.import_meta_graph(self.model_path+'Para.ckpt-'+model+'.meta')
				saver.restore(sess,self.model_path+'Para.ckpt-'+model)
				graph = tf.get_default_graph()
				x=graph.get_operation_by_name('x').outputs[0]
				mask = graph.get_operation_by_name('mask').outputs[0]
				istrain=graph.get_operation_by_name('Placeholder').outputs[0]
				realpred=graph.get_operation_by_name('Sigmoid').outputs[0]
				predresult=sess.run(realpred,feed_dict={x: feature,istrain:False})
				predresult=predresult[0,:,:,0]
				predresult=(predresult+predresult.T)/2

				ensemble.append(predresult)
		final=(ensemble[0]+ensemble[1]+ensemble[2]+ensemble[3]+ensemble[4])/5.
		np.savetxt(result_path+protein+'.contactP',final)





if __name__=='__main__':
	pro_name=sys.argv[1]

	if not os.path.exists('./f_model/C/data'):
		os.mkdir('./f_model/C/data')
	os.chdir('scripts')
	os.system('python extractSSE_from_Spider3.py '+pro_name)
	os.system('python extractASA_from_Spider3.py '+pro_name)
	os.chdir('../')
	os.system('cp data/'+pro_name+'.SSE ./f_model/C/data')
	os.system('cp data/'+pro_name+'.ccmpred ./f_model/C/data')

	os.chdir('./f_model/C/feature_extraction/')
	os.system('python extract_features.py '+pro_name)
	os.chdir('../../../')
	deepccon=DeepCCon('./f_model/C/model/')
	deepccon.prediction(pro_name,'./f_model/C/features/','./f_model/C/results/')


	if not os.path.exists('./f_model/R/data'):
		os.mkdir('./f_model/R/data')
	os.system('cp data/'+pro_name+'.SSE ./f_model/R/data')
	os.system('cp data/'+pro_name+'.ACC ./f_model/R/data')
	os.system('cp data/'+pro_name+'.ccmpred ./f_model/R/data')
	os.system('cp data/'+pro_name+'.PSSM ./f_model/R/data')
# os.system('cp inter_process/AA_SSCP_db ./f_model/R/data')
	os.system('cp ./f_model/C/results/'+pro_name+'.DeepCCon ./f_model/R/data')
	os.chdir('./f_model/R/feature_extraction/')
	os.system('python extract_features_long.py '+pro_name)
	os.system('python extract_features_medium.py '+pro_name)
	os.system('python extract_features_short.py '+pro_name)
	os.chdir('../../../')
	deeprcon=DeepRCon('./f_model/R/model/model/')
	deeprcon.prediction(pro_name,'./f_model/R/features/','./f_model/R/results/')
	

	os.system('cp data/'+pro_name+'.ccmpred ./f_model/New_Third/data/ccmpred')
	os.system('cp ./f_model/R/results/long/'+pro_name+'.DeepRCon ./f_model/New_Third/data/long')
	os.system('cp ./f_model/R/results/medium/'+pro_name+'.DeepRCon ./f_model/New_Third/data/medium')
	os.system('cp ./f_model/R/results/short/'+pro_name+'.DeepRCon ./f_model/New_Third/data/short')
	os.chdir('scripts')	
	os.system('python ./extract_prob.py '+pro_name)
	os.chdir('../')
	new_third=new_Third('./f_model/New_Third/model/')
	test_data=new_third.data_process(pro_name,'./f_model/New_Third/data/')
	new_third.prediction(pro_name,'./f_model/New_Third/results/',test_data)
	os.system('cp f_model/New_Third/results/'+pro_name+'.contactP ./result/' )



		
