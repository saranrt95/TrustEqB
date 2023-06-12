# TrustEqB
This repository contains the code and data to reproduce the experiments carried out in the paper: _Narteni S., Muselli M., Dabbene F., Mongelli M. (2023). Trustworthy Artificial Intelligence Classification-based Equivalent Bandwidth Control, currently under review for Elsevier-Computer Communications journal_.
It deals with the combined usage of control and rule-based classification for the EqB allocation, in the context of trustworthy AI.
Clopper-Pearson generalization bound is used as an efficient tool to select a rule-based model that performs adequately,
also determining the minimum amount of data required for model training. Robustness, in terms of the model’s ability to recognize out-of-distribution samples is investigated, by comparing the different rates of satisfaction of rules in presence of training or operational data, which is quantified via simple statistics, mutual information, l1 and l2 norms.


# Description
- db building for training.cpp: simulation code to generate the dataset (uses sobol6x10000.txt)
- OutTrain6.txt: the dataset coming from simulation
- train3000_rules.csv: the XAI model for EqB control
- ClopperPearsonGeneralization folder: data and code for model selection via Clopper-Pearson bound
- control application.cpp: code that applies the XAI model to the bandwidth control problem and generates the number of hits for robustness tests
- RobustnessTest: data and code for out-of-distribution detection


