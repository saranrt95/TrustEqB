import pandas as pd
import numpy as np
#from sklearn.metrics import log_loss
# load data with test set (of size m) and LLM results (class probabilities)
def load_data(filename):
	if filename[-4:]=='.csv':
		data=pd.read_csv(filename)
		#print(data)
		#y_data=data[class_label]
	else:
		if filename[-5:]=='.xlsx':
			data=pd.read_excel(filename)
			#y_data=data[class_label]
		else:
			if filename[-4:]=='.txt':
				data=pd.read_csv(filename,delimiter="\t")
				#print(data)
				#y_data=data[class_label]
	return data

def get_class_probabilities(test_data,nclasses):
	probabs_test=test_data[test_data.columns[-nclasses:]]# get probabilities predicted for each train data for each class;
	probabs_testlist=[]
	for i in range(nclasses):
		probabs_testlist.append(list(probabs_test[probabs_test.columns[i]]))
	probs=list(zip(*probabs_testlist))
	#print(probs)
	return probs
	'''
	class_probs=[]
	for p in probs:
		class_probs.append(max(p))
	return class_probs	
	'''

# cross-entropy loss (non la usiamo)
def cross_entropy_loss(y, probs):
	losstot=[]
	for c,p in zip(tuple(y),probs):
		#print("target label (one-hot encoded), class probabilities")
		#print((c,p))
		loss=0
		for cl in range(len(c)):
			'''
			print("label:")
			print(cl)
			print("predicted prob for label "+str(c[cl]))
			print(p[cl])
			'''
			if c[cl]==1:
				#print("loss value")
				#print(-np.log(p[cl]))
				ls=-np.log(p[cl])
				#loss.append(-np.log(p[cl]))
			else:

				ls=-np.log(1-p[cl])
				#print("loss value")
				#print(ls)
				#loss.append(-np.log(1-p[cl]))
			loss+=ls
		losstot.append(loss)
	return losstot


# 0-1 loss
def zero_one_loss(y_true,y_pred):
	loss=[]
	for i in range(len(y_true)):
		#print((y_true[i],y_pred[i]))
		if y_true[i]==y_pred[i]:
			loss.append(0)
		else:
			loss.append(1)
	#print(loss)
	return loss

def get_empirical_risk(loss,test_size):
	return np.sum(loss)/test_size

def clopper_pearson_bound(emp_risk, test_size,epsilon):
	return emp_risk+np.sqrt(np.log(1/epsilon)/(2*test_size))





