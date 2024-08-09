clear all;

clear all; 
stimFolder = './neurophys_analysis/dataCND/';
load([stimFolder, '/dataStimComplete.mat'])

clear corrTr
for trial=1:size(stim.data, 2)
    bert = stim.data{14, trial}(:, 3);
    clip = stim.data{14, trial}(:, 5);
    r = corrcoef(bert(bert>0), clip(clip>0));
    corrTr(trial) = r(1,2);
end

display(["BERT vs CLIP", mean(corrTr)]);

clear corrTr
for trial=1:size(stim.data, 2)
    bert = stim.data{14, trial}(:, 3);
    brod = stim.data{14, trial}(:, 4);
    indzeroe = xor(brod>0,bert>0); 
    
    % brod(indzeroe)=0
    % bert(indzeroe)=0
    r = corrcoef(bert(bert>0), brod(brod>0));
    corrTr(trial) = r(1,2);
end

display(["BERT vs Brod", mean(corrTr)]);

clear corrTr
for trial=1:size(stim.data, 2)
    clip = stim.data{14, trial}(:, 5);
    brod = stim.data{14, trial}(:, 4);
    r = corrcoef(clip(clip>0), brod(brod>0));
    corrTr(trial) = r(1,2);
end

display(["CLIP vs Brod",mean(corrTr)]);