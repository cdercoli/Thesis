dataMainFolder = "./neurophys_analysis/";
dataCNDSubfolder = "dataCND/";
stimFilename = dataMainFolder+dataCNDSubfolder+"dataStim.mat";
load(stimFilename,'stim')

dissFile = "./dissimilarity_analysis/text_data/AliceChapterOne.csv";

for tr = 1:12 
    
    t = readtable(dissFile,'Delimiter', ',');

    tSt = table2array(t(:,3));
    tFin = table2array(t(:,4));
    BERTdiss = table2array(t(:,18));
    CLIPdiss = table2array(t(:,20));

    fs = stim.fs;
    wordOnsetIdx = 2;
    
    x = stim.data{wordOnsetIdx,tr};

    xDissBERT = x;
    xDissCLIP = x;

    xDissBERT(xDissBERT~=0) = BERTdiss(1:nnz(x==1));
    xDissCLIP(xDissCLIP~=0) = CLIPdiss(1:nnz(x==1));

    env = stim.data{1,tr}; 
    x = stim.data{2,tr}; 

    stim.names{3} = 'dissBERT';
    stim.data{3,tr} = xDissBERT;

    stim.names{4} = 'dissCLIP';
    stim.data{4,tr} = xDissCLIP;

    stim.names{5} = 'env+onset+BERT';
    stim.data{5,tr} = [env, x, xDissBERT];
   
    stim.names{6} = 'env+onset+CLIP';
    stim.data{6,tr} = [env, x, xDissCLIP];

    stim.names{7} = 'env+onset+BERT+CLIP';
    stim.data{7,tr} = [env, x, xDissBERT,xDissCLIP];

end

% Save the stim file
save(dataMainFolder+dataCNDSubfolder+'/dataStimComplete.mat', 'stim')