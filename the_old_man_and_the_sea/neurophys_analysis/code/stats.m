%load chanlocs
load('./neurophys_analysis/dataCND/chanlocs128.mat')

%load clip+bert model
load('./neurophys_analysis/results/BERT_CLIP_results.mat')
rOrigCLIPBERT = rAllElec;
rTCCLIPBERT = rAllElecTC;
rGT_TCCLIPBERT = rAllElecGT_TC;
rGTCLIPBERT = rAllElecGT;
modelTrueCLIPBERT = modelAll;
modelTrueAvgCLIPBERT = mTRFmodelAvg(modelAll, 1);

%load shuffled bert in clip+bert  model
load('./neurophys_analysis/results/shBert_results_BERTCLIP.mat')
rOrigshBERT = rAllElec;
rTCShBERT = rAllElecTC;
rGT_TCShBERT = rAllElecGT_TC;
rGTShBERT = rAllElecGT;
modelTrueshBERT = modelAll;
modelTrueAvgshBERT = mTRFmodelAvg(modelAll, 1);

rOrigDiffBERT = rOrigCLIPBERT-rOrigshBERT;
rTCDiffBERT  = rTCCLIPBERT-rTCShBERT;
rGT_TCDiff_BERT = rGT_TCCLIPBERT-rGT_TCShBERT;
rGTDiffBert  = rGTCLIPBERT-rGTShBERT;

%load shuffled clip in clip+bert  model
load('./neurophys_analysis/results/shCLIP_results_BERTCLIP.mat')
rOrigshCLIP_CLIPBERT = rAllElec;
rTCShCLIP_CLIPBERT = rAllElecTC;
rGT_TCshCLIP_CLIPBERT = rAllElecGT_TC;
rGTshCLIP_CLIPBERT = rAllElecGT;
modelTrueshCLIP_CLIPBERT = modelAll;
modelTrueAvgshCLIP_CLIPBERT = mTRFmodelAvg(modelAll, 1);

%load shuffled clip and bert in clip+bert  model
load('./neurophys_analysis/results/shBothBERTCLIP_results_BERTCLIP.mat')
rOrigshBothBERTCLIP = rAllElec;
rTCshBothBERTCLIP = rAllElecTC;
rGT_TCshBothBERTCLIP = rAllElecGT_TC;
rGTshBothBERTCLIP = rAllElecGT;
modelTrueshBothBERTCLIP = modelAll;
modelTrueAvgshBothBERTCLIP = mTRFmodelAvg(modelAll, 1);


%% 
% for the error bars
means = [mean(mean(rTCCLIPBERT)), mean(mean(rTCShCLIP_CLIPBERT)), mean(mean(rTCShBERT)), mean(mean(rTCshBothBERTCLIP))];
errors = [std(mean(rTCCLIPBERT))/sqrt(numel(rTCCLIPBERT)), ...
    std(mean(rTCShCLIP_CLIPBERT))/sqrt(numel(rTCShCLIP_CLIPBERT)), ...
    std(mean(rTCShBERT))/sqrt(numel(rTCShBERT)), ...
    std(mean(rTCshBothBERTCLIP))/sqrt(numel(rTCshBothBERTCLIP))];

bar(1:4,[mean(mean(rTCCLIPBERT)), mean(mean(rTCShCLIP_CLIPBERT)), mean(mean(rTCShBERT)), mean(mean(rTCshBothBERTCLIP))])
hold on
er = errorbar(1:4, means, errors, 'LineStyle', 'none', 'Color', 'k');
xticklabels({'Acoustics Text and Visual', 'Acoustics and Text', 'Acoustics and Visual', 'Acoustics'})
ax = gca; 
ax.FontSize = 16; 
ylim([0.0415 0.0441])
ylabel("Prediction Correlations")
xlabel("Models")

mean_rTCCLIPBERT = mean(rTCCLIPBERT, 1);
mean_rTCShBERT = mean(rTCShBERT, 1);
mean_rTCShCLIP = mean(rTCShCLIP_CLIPBERT, 1);
mean_rTCshBothBERTCLIP = mean(rTCshBothBERTCLIP, 1);

rtcc_vals = {mean_rTCShCLIP, mean_rTCShBERT, mean_rTCshBothBERTCLIP}; 
rtcc_names = {"ShCLIP", "ShBERT", "ShBothBERTCLIP"}; 

for x = 1: length(rtcc_vals)
    [h, p, ci, s] = ttest(mean_rTCCLIPBERT, rtcc_vals{x});
    
    if h == 0
        disp('CLIPBERT and '+ string(rtcc_names(x)) + '  are not significantly different.');
    else
        disp('CLIPBERT and ' + string(rtcc_names(x)) + ' are significantly different.');
    end
    disp(['p-value: ', num2str(p)]);
    disp(['t-statistic: ', num2str(s.tstat)]);
    disp(['Confidence Interval: ', num2str(ci(1)), ' to ', num2str(ci(2))]);
end


%%
bert_samples_distrib = cell(1, 104);
clip_samples_distrib = cell(1, 104);

for t = 1:104
    bert_vals = [];
    clip_vals = [];
    
    for sub = 1:19
        elec = 65; % plotted on the trf
        bert_vals = [bert_vals, modelTrueCLIPBERT(sub).w(3,t,elec)];
        clip_vals = [clip_vals, modelTrueCLIPBERT(sub).w(4,t,elec)];
    end
    bert_samples_distrib{t} = bert_vals;
    clip_samples_distrib{t} = clip_vals;
end

pvals = zeros(1, 104);
for i = 1:104
    [~, pvals(i)] = ttest(bert_samples_distrib{i}, clip_samples_distrib{i});
end

%note that the mafdr function works only if you have the bioinformatics
%toolbox (already installed with the academics licence)
[p_fdr, ~] = mafdr(pvals);

fprintf('pvals list length: %d\n', length(p_fdr));
fprintf('pvals mean: %f\n', mean(p_fdr));
fprintf('min val in pvals: %f\n', min(p_fdr));
fprintf('max val in pvals: %f\n', max(p_fdr));

signif_timepoints = find(p_fdr < 0.05);
disp('Significant time points:');
disp(signif_timepoints);

%% Plot p-values
figure;
y =  [];
for i = 1:104
    y = [y, modelTrueCLIPBERT(1).t(i)];
end
plot(y, pvals);
xlabel('Time Points');
ylabel('P-values');
title('P-values for each time point');

%% plot pvals on trf

%load clip+bert model
load('./neurophys_analysis/results/BERT_CLIP_results.mat')
rOrigCLIPBERT = rAllElec;
rTCCLIPBERT = rAllElecTC;
rGT_TCCLIPBERT = rAllElecGT_TC;
rGTCLIPBERT = rAllElecGT;
modelTrueCLIPBERT = modelAll;
modelTrueAvgCLIPBERT = mTRFmodelAvg(modelAll, 1);

for s=1:19
    new = modelTrueCLIPBERT(s).w;
    for ffea = 1:size(new,1)
        tmp = new(ffea,:,:);
        stdev = std(tmp(:));
        modelTrueCLIPBERT(s).w(ffea,:,:) = modelTrueCLIPBERT(s).w(ffea,:,:)./stdev;
    end
end

new = modelTrueAvgCLIPBERT.w;
for ffea = 1:size(new,1)
    tmp = new(ffea,:,:);
    stdev = std(tmp(:));
    modelTrueAvgCLIPBERT.w(ffea,:,:) = modelTrueAvgCLIPBERT.w(ffea,:,:)./stdev;
end

for s=1:19
    subjallweightsEnv(s,:,:) = modelTrueCLIPBERT(s).w(1,:,:);
    subjallweightsWordonset(s,:,:) = modelTrueCLIPBERT(s).w(2,:,:);
    subjallweightsBert(s,:,:) = modelTrueCLIPBERT(s).w(3,:,:);
    subjallweightsClip(s,:,:) = modelTrueCLIPBERT(s).w(4,:,:);
end

stderrorenv = std(subjallweightsEnv,0,1)/sqrt(19);
stderrorwordonset = std(subjallweightsWordonset,0,1)/sqrt(19);
stderrorClip = std(subjallweightsClip,0,1)/sqrt(19);
stderrorBert = std(subjallweightsBert,0,1)/sqrt(19);

%1 ENV
curve1 = modelTrueAvgCLIPBERT.w(1,:,65) + stderrorenv(1,:,65);
curve2 = modelTrueAvgCLIPBERT.w(1,:,65) - stderrorenv(1,:,65);
x = modelTrueAvgCLIPBERT.t;
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];
h = fill(x2, inBetween, "b");
xlim([0 600])
set(h,'facealpha',.5)
hold on
plot(modelTrueAvgCLIPBERT.t,squeeze(modelTrueAvgCLIPBERT.w(1,:,65)),'LineWidth',2)
hold on
%2 WO
curve1 = modelTrueAvgCLIPBERT.w(2,:,65) + stderrorwordonset(1,:,65);
curve2 = modelTrueAvgCLIPBERT.w(2,:,65) - stderrorwordonset(1,:,65);
inBetween = [curve1, fliplr(curve2)];
h = fill(x2, inBetween, "r");
set(h,'facealpha',.5)
hold on
plot(modelTrueAvgCLIPBERT.t,squeeze(modelTrueAvgCLIPBERT.w(2,:,65)),'LineWidth',2)
hold on
%3  BERT
curve1 = modelTrueAvgCLIPBERT.w(3,:,65) + stderrorBert(1,:,65);
curve2 = modelTrueAvgCLIPBERT.w(3,:,65) - stderrorBert(1,:,65);
inBetween = [curve1, fliplr(curve2)];
h = fill(x2, inBetween, "y");
set(h,'facealpha',.5)
hold on
plot(modelTrueAvgCLIPBERT.t,squeeze(modelTrueAvgCLIPBERT.w(3,:,65)),'LineWidth',2)
%4 CLIP
curve1 = modelTrueAvgCLIPBERT.w(4,:,65) + stderrorClip(1,:,65);
curve2 = modelTrueAvgCLIPBERT.w(4,:,65) - stderrorClip(1,:,65);
inBetween = [curve1, fliplr(curve2)];
h = fill(x2, inBetween, "m");
set(h,'facealpha',.5)
hold on
plot(modelTrueAvgCLIPBERT.t,squeeze(modelTrueAvgCLIPBERT.w(4,:,65)),'LineWidth',2)

%line for significance stats 
signif_timepoints_wrong = find(pvals < 0.05);
signif_timepoints = [];
for i = signif_timepoints_wrong
    signif_timepoints = [signif_timepoints, modelTrueCLIPBERT(1).t(i)];
end

signif_timepoints_filled = NaN * ones(size(modelAll(1).t));
[~, idx] = ismember(modelAll(1).t, signif_timepoints);
signif_timepoints_filled(idx > 0) = modelAll(1).t(idx > 0);
disp(signif_timepoints_filled);

y = -5.* ones(length(modelAll(1).t),1);

line(signif_timepoints_filled, y); 
yline(0, '--');
ylabel("Normalized amplitude")
xlabel("Time (ms)")
legend({'','Env', '', 'Wo', '', 'BERT Diss', '',  'CLIP Diss'})
title('TRF weights for CLIPBERT model for electrode Cz')
