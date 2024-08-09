clear all; 

dataMainFolder = "./neurophys_analysis/dataCND/";
dataDissSubfolder = "./dissimilarity_analysis/text_data/";
stimFilename = dataMainFolder+"dataStimBase.mat";
load(stimFilename,'stim')

for tr = 1:size(stim.data, 2)
    dissFilename = dataDissSubfolder+"file"+tr+"_diff.csv";
    t = readtable(dissFilename,'Delimiter', ',');
    
    %extracting relevant columns from dissimilarity csv files
    tSt = table2array(t(:,2));
    tFin = table2array(t(:,3));

    bertDiss_csv = table2array(t(:,6));
    clipDiss_csv = table2array(t(:,8));
    
    fs = stim.fs;
    wordOnsetIdx = 2;
    
    x = zeros(size(stim.data{wordOnsetIdx,tr})); % word onset vector
    x(1+round(tSt*fs),1) = 1;

    xDiss = x;
    clipDiss = x; 
    xDiss(xDiss~=0) = bertDiss_csv;
    clipDiss(clipDiss~=0) = clipDiss_csv;

    env = stim.data{1,tr}; 
    x = stim.data{2,tr}; 
    envDeriv = stim.data{5,tr};
    envDeriv(envDeriv < 0) = 0;
    BrodDiss = stim.data{6,tr};

    indicesToZero = xor(x > 0, xDiss > 0);
    xDiss(indicesToZero) = 0;
    BrodDiss(indicesToZero) = 0;
    clipDiss(indicesToZero) = 0;
    x(indicesToZero) = 0;

    %removing extra zeros (i.e.: coming from stopwords)
    rowsToRemove = (xDiss == 0 & clipDiss ~= 0) | (xDiss == 0 & BrodDiss ~= 0) | ...
               (clipDiss == 0 & xDiss ~= 0) | (clipDiss == 0 & BrodDiss ~= 0) | ...
               (BrodDiss == 0 & xDiss ~= 0) | (BrodDiss == 0 & clipDiss ~= 0);
    xDiss(rowsToRemove) = [];
    clipDiss(rowsToRemove) = [];
    BrodDiss(rowsToRemove) = [];
    x(rowsToRemove) = [];
    env(rowsToRemove) = [];
    envDeriv(rowsToRemove) = [];

    stim.names{7} = 'dissBERT';
    stim.data{7,tr} = xDiss;
  
    stim.names{8} = 'dissCLIP';
    stim.data{8,tr} = clipDiss;

    stim.names{9} = 'env+onset+BERT';
    stim.data{9,tr} = [env, x, xDiss];
    
    stim.names{10} = 'env+onset+CLIP';
    stim.data{10,tr} = [env, x, clipDiss];

    stim.names{11} = 'env+onset+BERT+CLIP+word2vec';
    stim.data{11,tr} = [env, x, xDiss, clipDiss, BrodDiss];
    
    stim.names{12} = 'env+onset+BERT+CLIP';
    stim.data{12,tr} = [env, x, xDiss,clipDiss];

    stim.names{13} = 'env+onset+BERT+brod';
    stim.data{13,tr} = [env, x, xDiss, BrodDiss];

    stim.names{14} = 'env+onset+BERT+envDerivative';
    stim.data{14,tr} = [env, x, xDiss, envDeriv];

    stim.names{15} = 'env+onset+BERT+brod+CLIP+envDerivative';
    stim.data{15,tr} = [env, x, xDiss, BrodDiss, clipDiss, envDeriv];

    stim.names{16} = 'env+onset+BERT+CLIP+envDerivative';
    stim.data{16,tr} = [env, x, xDiss, clipDiss, envDeriv];

    stim.names{17} = 'env+onset+BERT+brod+envDerivative';
    stim.data{17,tr} = [env, x, xDiss, BrodDiss, envDeriv];
    
end

% Save the stim file
save(dataMainFolder+'/dataStimComplete.mat', 'stim')
