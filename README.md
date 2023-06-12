# TrustEqB
This repository contains the code and data to reproduce the experiments carried out in the paper: 	Narteni S., Muselli M., Dabbene F., Mongelli M. (2023). Trustworthy Artificial Intelligence Classification-based Equivalent Bandwidth Control, currently under review for Elsevier-Computer Communications journal.
It deals with the combined usage of control and rule-based classification for the EqB allocation, in the context of trustworthy AI.
Clopper-Pearson generalization bound is used as an efficient tool to select a rule-based model that performs adequately,
also determining the minimum amount of data required for model training. Robustness, in terms of the model’s ability to recognize out-of-distribution samples is investigated, by comparing the different rates of satisfaction of rules in presence of training or operational data, which is quantified via simple statistics, mutual information, l1 and l2 norms.

The workflow is summarized in the picture below:
![workflow_new](https://github.com/saranrt95/TrustEqB/assets/77918497/7bcc257b-d489-4f73-807b-e1c4456971bc){width = 100 height = 300}

# Description
- OutTrain6.txt: the dataset coming from simulation
- train3000_rules.csv: the XAI model for EqB control
- ClopperPearsonGeneralization folder: data and code for model selection via Clopper-Pearson bound
- RobustnessTest: data and code for out-of-distribution detection

