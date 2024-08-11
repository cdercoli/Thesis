This folder contains three analyses carried out for the dissertation project "How do we think?"
from Chiara D'Ercoli, an MSc Computer Science - Intelligent Systems student at Trinity College, The University of Dublin for the academic year 2023/2024. 

--------------------------------------------------------------------------------


The analyses are organized as follows: 
- folder: The Old Man and The Sea - including the entire code and data necessary to carry out 
the analysis for the first dataset, from embedding extraction to TRF runs. 
- folder: The Adventures of Alice in Wonderland - including the entire code and data for the analysis
of the second dataset. 
- Anaconda notebook: lexical analysis performed on BERT and CLIP embeddings. 

The two folders follow the same structure, explained here below. 
Two folders are presented in each, namely one including the dissimilarity analysis, and the other one for the neurophysiological one. 

-----------------------------
Dissimilarity analysis folder:
- folder "text_data": it includes the original text files from the respective experiment; 
- bert_clip_embeddings_extraction.ipynb: a notebook where the extraction of both BERT and CLIP 
embeddings are extracted from the afore mentioned text files, and re-included in the same files
for the semantic dissimilarity calculation (also performed in this notebook). Here, the dissimilarity scores are added back to the initial text files, saved in the general folder 
dissimilarity_analysis. 

-----------------------------
Neurophysiological analysis folder: 
- folder "code": it includes the code relevant to the study, namely: 
	- addPathDependencies.m : matlab file necessary for running the forwardTRF model from the CNSP workshop: 
	- cndNormaliseChiara.m: normalization function utilized in the forwardTRF code; 
	- CNSP_example1_forwardTRF.m: code for loading, preprocessing EEG data and running the TRF
from the CNSP workshop; 
	- diss2stim.m: code for inserting the previously calculated dissimilarities into the stim file utilized as input of the TRF; 
	- dissimilarity_corr.m: code snippet to calculate the correlation between dissimilarity scores obtained from BERT and CLIP embeddings. 
	- fdr_bh.m, mTRFcrossval_multimetric.m, shuffleFeatures.m: support functions necessary when running the mTRF and the various models; 
	- stats.m: code snippet to calculate the statistics on the models obtained. 

- folder "dataCND": data employed in the study: 
	- folder "subject_data": EEG recordings from the respective experiment; 
	- chanlocs128.mat: channel file including electrode locations for plotting purposes; 
	- dataStim.mat: original file stim from the downloaded dataset; 
	- dataStimBase.mat: file stim with some additional features (i.e.: envelope derivative);
	- dataStimComplete.mat: file stim including features considering dissimilarity scores (obtained after running the diss2stim.m code). 
	
	Note that for the_old_man_and_the_sea folder, in dataCND no subject_data will be found because due to file size, the upload on GitHub was not possible. Therefore, a readme file including the link to the OneDrive download of the respective subject data folder has been included and the download should be performed before running the analysis. Apologies for the inconvenience. 

- folder "results": including the results that have been obtained.
----------------------------------------------------------

The flow for running an entire analysis would be the following:

1. Embeddings extraction: run all cells in bert_clip_embeddings_extraction.ipynb
2. Insert dissimilarity scores into stim file: run diss2stim.m
3. TRF model run: run CNSP_example1_forwardTRF.m
4. Run statistics on the obtained model: run stats.m

-----------------------------------------------------------

