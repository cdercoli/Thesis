function [cnd] = cndNormaliseChiara(cnd,celldim)
%CNDNORMALISE Rescales stimulus or neural data in CND format. The function
% will implement a number of normalisations and standardisations. This
% first version only implements one type of standardisation, consisting of
% dividing the data by the overall standard deviation i.e., std calculated
% across all dimensions (e.g., EEG channels). This procedure preserves the
% relative magnitude across feature dimensions (e.g., freq bands of a
% speech sgram, EEG channels).
%   CND = CNDNORMALISE(CND,CELLDIM) 
%
%      'cnd'     -- stimulus or neural data in the Continuous-event Neural
%                    Data format (CND)
%      'celldim' -- dimensions to normalise. If stimulus data, dimensions
%                   correspond to different stimulus feature-sets
%                   (e.g., sgram, word dissimilarity)
%
%   Note: cnd.data: featureSets x trials. 'celldim' refers to the first
%   dimension of that cell array.
%
%   Author: 
%   Last update:  
%   Copyright 2022 Di Liberto Lab, Trinity College Dublin

% stimFilename = fullfile("./BERT_analysis/dataCND/dissimilarity_calculation/dataStimClipBertBrod.mat")
% load(stimFilename,'stim')
% stimFeature = stim;
% cnd = stimFeature; 

    if isempty(cnd) || isempty(cnd.data)
        disp('The CND structure is empty or not a cell array')
        return
    elseif ~iscell(cnd.data)
        disp('The CND.data structure is not a cell array')
        return
    end
    
    % If celldim was unspecified or empty
    if ~exist('celldim') || isempty(celldim)
        celldim = 1:size(cnd.data,1);
    end
    
    % For each celldim (e.g., stimulus feature)
    for iCell = celldim
        % Get all values from all trials
        tempCnd = cnd.data{iCell,1};
        for tr = 2:length(cnd.data) 
            tempCnd = cat(1,tempCnd,cnd.data{iCell,tr});
        end
        normFactors = std(tempCnd,1); %if the cell has four columns, here there will be 4 values
        
        % standardisation may not be what we want to apply
        for col = size(tempCnd, 2) 
            curr = tempCnd(:,col); 
            if sum(curr==0)/length(curr) > 0.75
                disp(['Warning: dimension ',num2str(iCell),' has over 3 times more zeros than non-zero values.' ...
                    ' Switching to normalising only non-zero values.'])
                nonZerotmpCnd = curr(curr~=0);
                normFactors(col) = std(nonZerotmpCnd);
            end
            if normFactors(col) == 0 % Presumed onset feature
                continue;
            end
            % Rescaling per feature
            for tr = 1:length(cnd.data)
                cnd.data{iCell,tr}(:,col) = cnd.data{iCell,tr}(:,col)/normFactors(col);
            end
        end   
        clear tempCnd;
    end


