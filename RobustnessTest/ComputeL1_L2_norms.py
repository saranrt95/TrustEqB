import numpy as np

from numpy.linalg import norm

import pandas as pd
#import time # for debug
import itertools as itr


def GetRuleHitRate(nrules,hitx,hity):
	hitxnorm = hitx/nrules
	hitynorm = hity/nrules
	return hitxnorm, hitynorm


def GetHitRate(hitx,hity):
	hitx = hitx/np.sum(hitx)
	hity = hity/np.sum(hity)
	return hitx,hity

def ComputeNorms(hitdata,x,y):
	hitx=np.array(hitdata[x])
	hity=np.array(hitdata[y])
	hitx_v,hity_v=GetHitRate(hitx,hity)
	l1norm=norm(hitx_v-hity_v,1)
	l2norm=norm(hitx_v-hity_v)
	return l1norm, l2norm


def ComputeVectorNormsInTrainAndOps(hitdata,training_idx,operational_idx):
	# input: DataFrame with number of hits for each rule and train/operational indexes (ID 1-ID 5 columns: training; ID 6-ID 10 cols: operational)
	# output: DataFrame with L1 and L2 norms for each considered couple of histograms
	L1list=[]
	L2list=[]
	coupleslist=[]
	for x in training_idx: 
		for y in operational_idx:
			l1norm, l2norm = ComputeNorms(hitdata,x,y)
			coupleslist.append((x,y))
			L1list.append(l1norm)
			L2list.append(l2norm)		
	normsdata=pd.DataFrame(zip(coupleslist,L1list,L2list),columns=["Histogram couples", "L1 norm", "L2 norm"])
	return normsdata 

def ComputeVectorNormsInTraining(hitdata,training_idx):
	L1list=[]
	L2list=[]
	coupleslist=[]
	for tr_couple in list(itr.combinations(training_idx,2)):
		x=tr_couple[0]
		y=tr_couple[1]
		l1norm, l2norm = ComputeNorms(hitdata,x,y)
		coupleslist.append((x,y))
		L1list.append(l1norm)
		L2list.append(l2norm)
	normsdata=pd.DataFrame(zip(coupleslist,L1list,L2list),columns=["Histogram couples", "L1 norm", "L2 norm"])
	return normsdata

# MAIN
# load data with histograms
hitdata=pd.read_excel("NumberOfHits_fi34_fi05.xlsx")
#print(hitdata)


# set training columns
training_idx = ['ID 1','ID 2','ID 3','ID 4','ID 5']
# set operational columns
operational_idx = ['ID 6','ID 7','ID 8','ID 9','ID 10']

normsTrainOps=ComputeVectorNormsInTrainAndOps(hitdata, training_idx, operational_idx)
print("############## TRAINING VS OPERATIONAL ###############")
print(normsTrainOps)

normsTrainOps.to_excel("l1_l2norms_train_ops_fi34_fi05.xlsx",index=False)

normsInTrain = ComputeVectorNormsInTraining(hitdata, training_idx)
print("############## TRAINING vs TRAINING ###############")
print(normsInTrain) 

normsInTrain.to_excel("l1_l2norms_intrain_fi34.xlsx",index=False)
