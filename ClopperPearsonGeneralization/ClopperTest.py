from ClopperPearsonBound import *
from os import listdir
from os.path import isfile, join

# MAIN

#load data
# directory with LLM results
datadir='./dataFiles/'
class_label='out'
predicted_label='pred(out)'
epsilon=[0.1, 0.05, 0.02, 0.01]
datafiles = [f for f in listdir(datadir) if isfile(join(datadir, f)) and f!='.DS_Store']
results=[]
for file in datafiles:
	ID=int(file.split('_')[0][2:])
	if ID >=5 and ID<=10:
		nclasses=2
	else:
		nclasses=3

	data=load_data(datadir+file)
	test_size=len(data)
	# 0-1 loss
	loss=zero_one_loss(data[class_label],data[predicted_label])

	empirical_risk=get_empirical_risk(loss, test_size)
	for eps in epsilon:
		CP_bound=clopper_pearson_bound(empirical_risk,test_size,eps)
		d={'ID':[ID],'TestSize':[test_size],'Epsilon': [eps],'ClopperPearson': [CP_bound], 'EmpiricalRisk': empirical_risk}
		res=pd.DataFrame(data=d)
		results.append(res)
finalRes=pd.concat(results)
finalRes.to_excel('CPresults.xlsx')


