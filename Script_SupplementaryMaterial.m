%% This script contains the codes used to generate and plot the various 
% representations discussed in the Supplementary Material of the article
% "The arrow-of-time in neuroimaging time series identifies causal triggers
% of brain function"



%% Definition of parameters

CM = cbrewer('qual','Set1',8);
CM_RB = flipud(cbrewer('div','RdBu',1000));
CM_RB(CM_RB<0) = 0;
CM_DYN_FC = cbrewer('seq','Reds',1000);
CM_DYN_FC(CM_DYN_FC<0) = 0;
CM_DYN_AOT = flipud(cbrewer('div','PuOr',1000));
CM_DYN_AOT(CM_DYN_AOT < 0) = 0;
CM_DYN_AOT(CM_DYN_AOT > 1) = 1;

CM_TCs = cbrewer('qual','Set3',10);
CM_TCs(CM_TCs<0) = 0;
CM_TCs(CM_TCs>1) = 1;

CM_TCs2 = flipud(cbrewer('div','BrBG',1000));
CM_TCs2(CM_TCs2 < 0) = 0;
CM_TCs2(CM_TCs2 > 1) = 1;

CM_boxplots = 1/255*[182,213,235; 181,221,197; 181,221,197; 231,194,178; 252,198,174; 255,244,220; 241,175,218; 249,182,186; 212,201,229];


% Number of regions in our main atlas
n_regions = 419;

% Number of subjects available in total
n_subjects = 100;

% Order of the MAR model to consider for the autoregressive modeling in the
% main results
Order = 1;

% Number of runs to select and concatenate for the main results
run_ids = 1;

% Number of folds to bootstrap for
n_folds = 100;

% Random subject orderings to consider
ridx = randperm(n_folds,n_subjects);

% Number of nulls
n_null = 35;

% Do we pick random subjects when concatenating data? 
is_random = 1;

% TR of the data
TR = 0.72;

% Numbers of samples available for each paradigm
n_samples_RS = 1200;
n_samples_EMOTION = 176;
n_samples_LANGUAGE = 316;
n_samples_MOTOR = 284;
n_samples_SOCIAL = 274;
n_samples_WM = 405;
n_samples_GAMBLING = 253;
n_samples_RELATIONAL = 232;

% Parameters of the sliding window scheme used for dynamic AoT tracking:
% window size (W, in samples) and shift Delta (in samples)
W = 20;
Delta = 1;

% Do we plot stuff?
is_plot = 0;



%% SF1. Results for resting state with more samples
% For this, we consider the resting state recording, which is by far the
% longest available one. We see how the estimates vary when using as many
% samples as available (i.e., up to 100000).

% Loads the appropriate data
load('SM_RS_MORE_SAMPLES.mat');

DK = squeeze(median(DK_RS,2));
NS = n_samples_tot;

% Mean and SEM across regions
Mu_AOT = mean(DK,2);
SEM_AOT = std(DK,[],2)/sqrt(n_regions);

% Plots the AoT evolution with more samples
figure;
hold on
set(gca,'Box','off');
plot(NS,Mu_AOT,'color',CM_boxplots(1,:),'LineWidth',2);
errorbar(NS,Mu_AOT,SEM_AOT,'LineStyle','None','Color',[0.4,0.4,0.4]);
plot(NS(7),Mu_AOT(7),'rs','MarkerFaceColor','r');

% Plots the heatmaps
figure;
hold on
imagesc(1:n_regions,NS,DK);
set(gca,'Box','off');
colormap(CM_RB);
xlabel('Brain region index');
ylabel('Total number of samples');
xlim([1,n_regions]);
caxis([-max(abs(DK(:))),max(abs(DK(:)))]);
ylim([NS(1),NS(end)]);
plot([1,n_regions],[8000,8000],'k--');

% Plots the correlation
for i = 1:length(NS)-1
    Sim(i) = corr(DK(i,:)',DK(i+1,:)');
end

figure;
hold on
plot(NS(2:end),Sim,'LineWidth',2,'Color',CM_boxplots(1,:));
plot([8000,8000],[0.75,1],'k--');
plot(NS(7),Sim(6),'rs','MarkerFaceColor','r');
set(gca,'Box','off');



%% SF2. Correlation plots and heatmaps for all paradigms

load('DATA_FIGURE2.mat');


% MOTOR
DK = squeeze(median(DK_MOTOR));

% Plots the heatmaps
figure;
hold on
imagesc(1:n_regions,NS,DK);
set(gca,'Box','off');
colormap(CM_RB);
xlabel('Brain region index');
ylabel('Total number of samples');
xlim([1,n_regions]);
caxis([-max(abs(DK(:))),max(abs(DK(:)))]);
ylim([NS(13),NS(end)]);
plot([1,n_regions],[8000,8000],'k--');

% Plots the correlation
for i = 1:length(NS)-1
    Sim(i) = corr(DK(i,:)',DK(i+1,:)');
end

figure;
hold on
plot(NS(2:end),Sim,'LineWidth',2,'Color',CM_boxplots(2,:));
plot([8000,8000],[0,1],'k--');
plot(NS(76),Sim(76),'rs','MarkerFaceColor','r');
set(gca,'Box','off');


% WM
DK = squeeze(median(DK_WM));

% Plots the heatmaps
figure;
hold on
imagesc(1:n_regions,NS,DK);
set(gca,'Box','off');
colormap(CM_RB);
xlabel('Brain region index');
ylabel('Total number of samples');
xlim([1,n_regions]);
caxis([-max(abs(DK(:))),max(abs(DK(:)))]);
ylim([NS(13),NS(end)]);
plot([1,n_regions],[8000,8000],'k--');

% Plots the correlation
for i = 1:length(NS)-1
    Sim(i) = corr(DK(i,:)',DK(i+1,:)');
end

figure;
hold on
plot(NS(2:end),Sim,'LineWidth',2,'Color',CM_boxplots(4,:));
plot([8000,8000],[0,1],'k--');
plot(NS(76),Sim(76),'rs','MarkerFaceColor','r');
set(gca,'Box','off');



% RELATIONAL
DK = squeeze(median(DK_RELATIONAL));

% Plots the heatmaps
figure;
hold on
imagesc(1:n_regions,NS,DK);
set(gca,'Box','off');
colormap(CM_RB);
xlabel('Brain region index');
ylabel('Total number of samples');
xlim([1,n_regions]);
caxis([-max(abs(DK(:))),max(abs(DK(:)))]);
ylim([NS(13),NS(end)]);
plot([1,n_regions],[8000,8000],'k--');

% Plots the correlation
for i = 1:length(NS)-1
    Sim(i) = corr(DK(i,:)',DK(i+1,:)');
end

figure;
hold on
plot(NS(2:end),Sim,'LineWidth',2,'Color',CM_boxplots(5,:));
plot([8000,8000],[0,1],'k--');
plot(NS(76),Sim(76),'rs','MarkerFaceColor','r');
set(gca,'Box','off');



% EMOTION
DK = squeeze(median(DK_EMOTION));

% Plots the heatmaps
figure;
hold on
imagesc(1:n_regions,NS,DK);
set(gca,'Box','off');
colormap(CM_RB);
xlabel('Brain region index');
ylabel('Total number of samples');
xlim([1,n_regions]);
caxis([-max(abs(DK(:))),max(abs(DK(:)))]);
ylim([NS(13),NS(end)]);
plot([1,n_regions],[8000,8000],'k--');

% Plots the correlation
for i = 1:length(NS)-1
    Sim(i) = corr(DK(i,:)',DK(i+1,:)');
end

figure;
hold on
plot(NS(2:end),Sim,'LineWidth',2,'Color',CM_boxplots(6,:));
plot([8000,8000],[0,1],'k--');
plot(NS(76),Sim(76),'rs','MarkerFaceColor','r');
set(gca,'Box','off');




% SOCIAL
DK = squeeze(median(DK_SOCIAL));

% Plots the heatmaps
figure;
hold on
imagesc(1:n_regions,NS,DK);
set(gca,'Box','off');
colormap(CM_RB);
xlabel('Brain region index');
ylabel('Total number of samples');
xlim([1,n_regions]);
caxis([-max(abs(DK(:))),max(abs(DK(:)))]);
ylim([NS(13),NS(end)]);
plot([1,n_regions],[8000,8000],'k--');

% Plots the correlation
for i = 1:length(NS)-1
    Sim(i) = corr(DK(i,:)',DK(i+1,:)');
end

figure;
hold on
plot(NS(2:end),Sim,'LineWidth',2,'Color',CM_boxplots(7,:));
plot([8000,8000],[0,1],'k--');
plot(NS(76),Sim(76),'rs','MarkerFaceColor','r');
set(gca,'Box','off');



% LANGUAGE
DK = squeeze(median(DK_LANGUAGE));

% Plots the heatmaps
figure;
hold on
imagesc(1:n_regions,NS,DK);
set(gca,'Box','off');
colormap(CM_RB);
xlabel('Brain region index');
ylabel('Total number of samples');
xlim([1,n_regions]);
caxis([-max(abs(DK(:))),max(abs(DK(:)))]);
ylim([NS(13),NS(end)]);
plot([1,n_regions],[8000,8000],'k--');

% Plots the correlation
for i = 1:length(NS)-1
    Sim(i) = corr(DK(i,:)',DK(i+1,:)');
end

figure;
hold on
plot(NS(2:end),Sim,'LineWidth',2,'Color',CM_boxplots(8,:));
plot([8000,8000],[0,1],'k--');
plot(NS(76),Sim(76),'rs','MarkerFaceColor','r');
set(gca,'Box','off');


% GAMBLING
DK = squeeze(median(DK_GAMBLING));

% Plots the heatmaps
figure;
hold on
imagesc(1:n_regions,NS,DK);
set(gca,'Box','off');
colormap(CM_RB);
xlabel('Brain region index');
ylabel('Total number of samples');
xlim([1,n_regions]);
caxis([-max(abs(DK(:))),max(abs(DK(:)))]);
ylim([NS(13),NS(end)]);
plot([1,n_regions],[8000,8000],'k--');

% Plots the correlation
for i = 1:length(NS)-1
    Sim(i) = corr(DK(i,:)',DK(i+1,:)');
end

figure;
hold on
plot(NS(2:end),Sim,'LineWidth',2,'Color',CM_boxplots(9,:));
plot([8000,8000],[0,1],'k--');
plot(NS(76),Sim(76),'rs','MarkerFaceColor','r');
set(gca,'Box','off');





%% SF3. Evolution when "no rest" task paradigms are taken

% ADD RELATIONAL

load('SM_NOREST_CURVES2.mat');

%%%%%% MOTOR

% Mean and SEM across regions
Mu_AOT = squeeze(mean(median(DK_MOTOR,2),3));
SEM_AOT = squeeze(std(median(DK_MOTOR,2),[],3))/sqrt(n_regions);

% Plots the AoT evolution with more samples
figure;
hold on
set(gca,'Box','off');
plot(n_samples_tot,Mu_AOT,'color',CM_boxplots(1,:),'LineWidth',2);
errorbar(n_samples_tot,Mu_AOT,SEM_AOT,'LineStyle','None','Color',[0.4,0.4,0.4]);
plot(n_samples_tot(16),Mu_AOT(16),'rs','MarkerFaceColor','r');

% Plots the heatmaps
tmp = squeeze(median(DK_MOTOR,2));

figure;
hold on
imagesc(1:n_regions,n_samples_tot,tmp);
set(gca,'Box','off');
colormap(CM_RB);
xlabel('Brain region index');
ylabel('Total number of samples');
xlim([1,n_regions]);
caxis([-max(tmp(:)),max(tmp(:))]);
ylim([n_samples_tot(1),n_samples_tot(end)]);
plot([1,n_regions],[8000,8000],'k--');


% Plots the correlation
for i = 1:length(n_samples_tot)-1
    Sim(i) = corr(tmp(i,:)',tmp(i+1,:)');
end

figure;
hold on
plot(n_samples_tot(2:end),Sim,'LineWidth',2,'Color',CM_boxplots(1,:));
plot([8000,8000],[0,1],'k--');
plot(n_samples_tot(16),Sim(15),'rs','MarkerFaceColor','r');
set(gca,'Box','off');


%%%%%% EMOTION

% Mean and SEM across regions
Mu_AOT = squeeze(mean(median(DK_EMOTION,2),3));
SEM_AOT = squeeze(std(median(DK_EMOTION,2),[],3))/sqrt(n_regions);

% Plots the AoT evolution with more samples
figure;
hold on
set(gca,'Box','off');
plot(n_samples_tot,Mu_AOT,'color',CM_boxplots(6,:),'LineWidth',2);
errorbar(n_samples_tot,Mu_AOT,SEM_AOT,'LineStyle','None','Color',[0.4,0.4,0.4]);
plot(n_samples_tot(16),Mu_AOT(16),'rs','MarkerFaceColor','r');

% Plots the heatmaps
tmp = squeeze(median(DK_EMOTION,2));

figure;
hold on
imagesc(1:n_regions,n_samples_tot,tmp);
set(gca,'Box','off');
colormap(CM_RB);
xlabel('Brain region index');
ylabel('Total number of samples');
xlim([1,n_regions]);
caxis([-max(tmp(:)),max(tmp(:))]);
ylim([n_samples_tot(1),n_samples_tot(end)]);
plot([1,n_regions],[8000,8000],'k--');


% Plots the correlation
for i = 1:length(n_samples_tot)-1
    Sim(i) = corr(tmp(i,:)',tmp(i+1,:)');
end

figure;
hold on
plot(n_samples_tot(2:end),Sim,'LineWidth',2,'Color',CM_boxplots(6,:));
plot([8000,8000],[0,1],'k--');
plot(n_samples_tot(16),Sim(15),'rs','MarkerFaceColor','r');
set(gca,'Box','off');



%%%%%% SOCIAL

% Mean and SEM across regions
Mu_AOT = squeeze(mean(median(DK_SOCIAL,2),3));
SEM_AOT = squeeze(std(median(DK_SOCIAL,2),[],3))/sqrt(n_regions);

% Plots the AoT evolution with more samples
figure;
hold on
set(gca,'Box','off');
plot(n_samples_tot,Mu_AOT,'color',CM_boxplots(7,:),'LineWidth',2);
errorbar(n_samples_tot,Mu_AOT,SEM_AOT,'LineStyle','None','Color',[0.4,0.4,0.4]);
plot(n_samples_tot(16),Mu_AOT(16),'rs','MarkerFaceColor','r');

% Plots the heatmaps
tmp = squeeze(median(DK_SOCIAL,2));

figure;
hold on
imagesc(1:n_regions,n_samples_tot,tmp);
set(gca,'Box','off');
colormap(CM_RB);
xlabel('Brain region index');
ylabel('Total number of samples');
xlim([1,n_regions]);
caxis([-max(tmp(:)),max(tmp(:))]);
ylim([n_samples_tot(1),n_samples_tot(end)]);
plot([1,n_regions],[8000,8000],'k--');


% Plots the correlation
for i = 1:length(n_samples_tot)-1
    Sim(i) = corr(tmp(i,:)',tmp(i+1,:)');
end

figure;
hold on
plot(n_samples_tot(2:end),Sim,'LineWidth',2,'Color',CM_boxplots(7,:));
plot([8000,8000],[0,1],'k--');
plot(n_samples_tot(16),Sim(15),'rs','MarkerFaceColor','r');
set(gca,'Box','off');



%%%% LANGUAGE


% Mean and SEM across regions
Mu_AOT = squeeze(mean(median(DK_LANGUAGE,2),3));
SEM_AOT = squeeze(std(median(DK_LANGUAGE,2),[],3))/sqrt(n_regions);

% Plots the AoT evolution with more samples
figure;
hold on
set(gca,'Box','off');
plot(n_samples_tot,Mu_AOT,'color',CM_boxplots(8,:),'LineWidth',2);
errorbar(n_samples_tot,Mu_AOT,SEM_AOT,'LineStyle','None','Color',[0.4,0.4,0.4]);
plot(n_samples_tot(16),Mu_AOT(16),'rs','MarkerFaceColor','r');

% Plots the heatmaps
tmp = squeeze(median(DK_LANGUAGE,2));

figure;
hold on
imagesc(1:n_regions,n_samples_tot,tmp);
set(gca,'Box','off');
colormap(CM_RB);
xlabel('Brain region index');
ylabel('Total number of samples');
xlim([1,n_regions]);
caxis([-max(tmp(:)),max(tmp(:))]);
ylim([n_samples_tot(1),n_samples_tot(end)]);
plot([1,n_regions],[8000,8000],'k--');


% Plots the correlation
for i = 1:length(n_samples_tot)-1
    Sim(i) = corr(tmp(i,:)',tmp(i+1,:)');
end

figure;
hold on
plot(n_samples_tot(2:end),Sim,'LineWidth',2,'Color',CM_boxplots(8,:));
plot([8000,8000],[0,1],'k--');
plot(n_samples_tot(16),Sim(15),'rs','MarkerFaceColor','r');
set(gca,'Box','off');



%%%% WM


% Mean and SEM across regions
Mu_AOT = squeeze(mean(median(DK_WM,2),3));
SEM_AOT = squeeze(std(median(DK_WM,2),[],3))/sqrt(n_regions);

% Plots the AoT evolution with more samples
figure;
hold on
set(gca,'Box','off');
plot(n_samples_tot,Mu_AOT,'color',CM_boxplots(4,:),'LineWidth',2);
errorbar(n_samples_tot,Mu_AOT,SEM_AOT,'LineStyle','None','Color',[0.4,0.4,0.4]);
plot(n_samples_tot(16),Mu_AOT(16),'rs','MarkerFaceColor','r');

% Plots the heatmaps
tmp = squeeze(median(DK_WM,2));

figure;
hold on
imagesc(1:n_regions,n_samples_tot,tmp);
set(gca,'Box','off');
colormap(CM_RB);
xlabel('Brain region index');
ylabel('Total number of samples');
xlim([1,n_regions]);
caxis([-max(tmp(:)),max(tmp(:))]);
ylim([n_samples_tot(1),n_samples_tot(end)]);
plot([1,n_regions],[8000,8000],'k--');


% Plots the correlation
for i = 1:length(n_samples_tot)-1
    Sim(i) = corr(tmp(i,:)',tmp(i+1,:)');
end

figure;
hold on
plot(n_samples_tot(2:end),Sim,'LineWidth',2,'Color',CM_boxplots(4,:));
plot([8000,8000],[0,1],'k--');
plot(n_samples_tot(16),Sim(15),'rs','MarkerFaceColor','r');
set(gca,'Box','off');



%% SF4. AoT patterns for all paradigm cases

% COMPLETE WITH RAPHAEL



%% SF5. Null data are very similar regardless of paradigm and 
% inclusion/exclusion of baseline epochs (for task recordings)

load('SM_NULL.mat');

Paradigms = {'RS','MOTOR','MOTOR2','WM','RELATIONAL','EMOTION','SOCIAL','LANGUAGE','GAMBLING'};

for p = 1:length(Paradigms)
    
    tmp_name = ['DKN_',Paradigms{p}];
    
    % Gets the data in tmp
    tmp = eval(tmp_name);
    
    % Gets the median across folds and then concatenates across null
    % realizations
    tmp = squeeze(median(tmp,1));
    tmp = tmp(:);

    NullData(:,p) = tmp;
end

CM_boxplots = 1/255*[182,213,235; 181,221,197; 181,221,197; 231,194,178; 252,198,174; 255,244,220; 241,175,218; 249,182,186; 212,201,229];

figure;
boxplot(NullData,'plotstyle','compact','colors',CM_boxplots,'medianstyle','line','outliersize',3,'symbol','.','whisker',15);
hold on
for p = 1:length(Paradigms)
    plot(p,mean(NullData(:,p)),'sk','MarkerSize',3,'MarkerFaceColor','k');
end


%% SF6. Supplementary analyses: results with coarser or finer scale atlas

% Loads the right data
load('Data_SCH200.mat','TS_MOTOR_GSR','TS_RS_GSR','Scrubbing_RS','Scrubbing_MOTOR');
n_regions = 219;

% How many runs do we wish to concatenate in total (max amount
% available)
n_runs_toreach = length(run_ids)*n_subjects;

for f = 1:n_folds
    
    f
    
    % For each fold, we have a different subject ordering, stored in ridx.
    % We want to consider the same subject ordering for all sample number
    % values (second loop below). We also want to re-consider the same
    % ordering in other sections of the work, so we store the random
    % orderings across folds
    ridx(f,:) = randperm(n_subjects);

    % Resting state data generation (concatenation across subjects)
    [Data_RS,SL_RS,TL_RS,RL_RS] = ...
    AoT_ConcatenateRuns(TS_RS_GSR,run_ids,n_subjects,100,...
    is_random,ridx(f,:),'None',n_regions,Scrubbing_RS,...
    'None',n_runs_toreach,[],[]);

    % Motor task
    [Data_MOTOR,SL_MOTOR,TL_MOTOR,RL_MOTOR] = ...
    AoT_ConcatenateRuns(TS_MOTOR_GSR,run_ids,n_subjects,100,...
    is_random,ridx(f,:),'None',n_regions,Scrubbing_MOTOR,...
    'None',n_runs_toreach,[],[]);
    
    n_TP = n_samples_opt/n_subjects/length(run_ids);

    % The metrics of interest are
    % computed for each paradigm; in addition to Delta kurtosis (DK*),
    % we also compute an alternative measure (DKL*), functional
    % connectivity (FC*) and effective connectivity as inferred through
    % MAR-1 coefficients (EC*)
    [DK_RS(f,:),DKL_RS(f,:),FC_RS(f,:),...
        EC_RS(f,:),NATP_RS{f}] = ...
        AoT_ComputeKurtosis_Full_Shiney(Data_RS,Order,n_TP,...
                SL_RS,TL_RS,RL_RS,'All');

    [DK_MOTOR(f,:),DKL_MOTOR(f,:),FC_MOTOR(f,:),...
        EC_MOTOR(f,:),NATP_MOTOR{f}] = ...
        AoT_ComputeKurtosis_Full_Shiney(Data_MOTOR,Order,n_TP,...
                SL_MOTOR,TL_MOTOR,RL_MOTOR,'All');
end



%% SF7-14. Extended comparison across all processing subchoices

[DK_SM_WM,DKL_SM_WM,NATP_SM_WM,Final_Data_WM,SimMat_WM] = ...
    AOT_ComputeSupplementaryAnalyses(TS_WM_GSR,TS_WM_NoRest_GSR,...
    TS_WM_NoGSR,TS_WM_NoRest_NoGSR,Scrubbing_WM,Scrubbing_WM_NoRest,...
    n_samples_opt,n_subjects,n_folds,n_regions,n_runs_toreach,is_random,Order);

[DK_SM_SOCIAL,DKL_SM_SOCIAL,NATP_SM_SOCIAL,Final_Data_SOCIAL,SimMat_SOCIAL] = ...
    AOT_ComputeSupplementaryAnalyses(TS_SOCIAL_GSR,TS_SOCIAL_NoRest_GSR,...
    TS_SOCIAL_NoGSR,TS_SOCIAL_NoRest_NoGSR,Scrubbing_SOCIAL,Scrubbing_SOCIAL_NoRest,...
    n_samples_opt,n_subjects,n_folds,n_regions,n_runs_toreach,is_random,Order);


[DK_SM_LANGUAGE,DKL_SM_LANGUAGE,NATP_SM_LANGUAGE,Final_Data_LANGUAGE,SimMat_LANGUAGE] = ...
    AOT_ComputeSupplementaryAnalyses(TS_LANGUAGE_GSR,TS_LANGUAGE_NoRest_GSR,...
    TS_LANGUAGE_NoGSR,TS_LANGUAGE_NoRest_NoGSR,Scrubbing_LANGUAGE,Scrubbing_LANGUAGE_NoRest,...
    n_samples_opt,n_subjects,n_folds,n_regions,n_runs_toreach,is_random,Order);


[DK_SM_EMOTION,DKL_SM_EMOTION,NATP_SM_EMOTION,Final_Data_EMOTION,SimMat_EMOTION] = ...
    AOT_ComputeSupplementaryAnalyses(TS_EMOTION_GSR,TS_EMOTION_NoRest_GSR,...
    TS_EMOTION_NoGSR,TS_EMOTION_NoRest_NoGSR,Scrubbing_EMOTION,Scrubbing_EMOTION_NoRest,...
    n_samples_opt,n_subjects,n_folds,n_regions,n_runs_toreach,is_random,Order);

% [DK_SM_RELATIONAL,DKL_SM_RELATIONAL,NATP_SM_RELATIONAL,Final_Data_RELATIONAL,SimMat_RELATIONAL] = ...
%     AOT_ComputeSupplementaryAnalyses(TS_RELATIONAL_GSR,TS_RELATIONAL_NoRest_GSR,...
%     TS_RELATIONAL_NoGSR,TS_RELATIONAL_NoRest_NoGSR,Scrubbing_RELATIONAL,Scrubbing_RELATIONAL_NoRest,...
%     n_samples_opt,n_subjects,n_folds,n_regions,n_runs_toreach,is_random,Order);

[DK_SM_MOTOR,DKL_SM_MOTOR,NATP_SM_MOTOR,Final_Data_MOTOR,SimMat_MOTOR] = ...
    AOT_ComputeSupplementaryAnalyses(TS_MOTOR_GSR,TS_MOTOR_NoRest_GSR,...
    TS_MOTOR_NoGSR,TS_MOTOR_NoRest_NoGSR,Scrubbing_MOTOR,Scrubbing_MOTOR_NoRest,...
    n_samples_opt,n_subjects,n_folds,n_regions,n_runs_toreach,is_random,Order);


[DK_SM_RS,DKL_SM_RS,NATP_SM_RS,Final_Data_RS,SimMat_RS] = ...
    AOT_ComputeSupplementaryAnalyses_RS(TS_RS_GSR,...
    TS_RS_NoGSR,Scrubbing_RS,...
    n_samples_opt,n_subjects,n_folds,n_regions,n_runs_toreach,is_random,Order);





%% SF15-18. Dynamic evolution for all time-locked paradigms

% Emotion paradigm
[DK_dyn_EMOTION,DKL_dyn_EMOTION,FC_dyn_EMOTION,EC_dyn_EMOTION,Act_dyn_EMOTION] = ...
    AoT_Compute_Dynamic_Evolution(TS_EMOTION_GSR,W,Delta,run_ids,Order);

figure;
imagesc(1:size(DK_dyn_EMOTION,2),1:n_regions,DK_dyn_EMOTION);
colormap(CM_RB);
caxis([-2,2]);
xlabel('Time [window index]');
ylabel('Regions');
set(gca,'Box','off');

figure;
imagesc(Paradigm1_EMOTION');
colormap(flipud(cbrewer('seq','Greys',1000)));
caxis([-0.6,1.6]);
set(gca,'Box','off');


% Working memory paradigm
[DK_dyn_WM,DKL_dyn_WM,FC_dyn_WM,EC_dyn_WM,Act_dyn_WM] = ...
    AoT_Compute_Dynamic_Evolution(TS_WM_GSR,W,Delta,run_ids,Order);

figure;
imagesc(1:size(DK_dyn_WM,2),1:n_regions,DK_dyn_WM);
colormap(CM_RB);
caxis([-2,2]);
xlabel('Time [window index]');
ylabel('Regions');
set(gca,'Box','off');

figure;
imagesc(Paradigm1_WM');
colormap(flipud(cbrewer('seq','Greys',1000)));
caxis([-0.6,1.6]);
set(gca,'Box','off');


% Social paradigm
[DK_dyn_SOCIAL,DKL_dyn_SOCIAL,FC_dyn_SOCIAL,EC_dyn_SOCIAL,Act_dyn_SOCIAL] = ...
    AoT_Compute_Dynamic_Evolution(TS_SOCIAL_GSR,W,Delta,run_ids,Order);

figure;
imagesc(1:size(DK_dyn_SOCIAL,2),1:n_regions,DK_dyn_SOCIAL);
colormap(CM_RB);
caxis([-2,2]);
xlabel('Time [window index]');
ylabel('Regions');
set(gca,'Box','off');

figure;
imagesc(Paradigm1_SOCIAL');
colormap(flipud(cbrewer('seq','Greys',1000)));
caxis([-0.6,1.6]);
set(gca,'Box','off');


% Language paradigm
[DK_dyn_LANGUAGE,DKL_dyn_LANGUAGE,FC_dyn_LANGUAGE,EC_dyn_LANGUAGE,Act_dyn_LANGUAGE] = ...
    AoT_Compute_Dynamic_Evolution(TS_LANGUAGE_GSR,W,Delta,run_ids,Order);

figure;
imagesc(1:size(DK_dyn_LANGUAGE,2),1:n_regions,DK_dyn_LANGUAGE);
colormap(CM_RB);
caxis([-2,2]);
xlabel('Time [window index]');
ylabel('Regions');
set(gca,'Box','off');

figure;
imagesc(Paradigm1_LANGUAGE');
colormap(flipud(cbrewer('seq','Greys',1000)));
caxis([-0.6,1.6]);
set(gca,'Box','off');


% Relational paradigm
[DK_dyn_RELATIONAL,DKL_dyn_RELATIONAL,FC_dyn_RELATIONAL,EC_dyn_RELATIONAL,Act_dyn_RELATIONAL] = ...
    AoT_Compute_Dynamic_Evolution(TS_RELATIONAL_GSR,W,Delta,run_ids,Order);

figure;
imagesc(1:size(DK_dyn_RELATIONAL,2),1:n_regions,DK_dyn_RELATIONAL);
colormap(CM_RB);
caxis([-2,2]);
xlabel('Time [window index]');
ylabel('Regions');
set(gca,'Box','off');

figure;
imagesc(Paradigm1_RELATIONAL');
colormap(flipud(cbrewer('seq','Greys',1000)));
caxis([-0.6,1.6]);
set(gca,'Box','off');
