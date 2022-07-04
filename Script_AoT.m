%% This script performs the computations for the "arrow of time" project%



%% 1. Loading of the data
% Here, we read the HCP data preprocessed by MIP post-docs into dedicated
% variables. We want to use non-filtered data, and to compare the outcomes
% with and without GSR. We use z-scored data, and we will compare different
% atlases: Schaefer 200, Schaefer 400, Schaefer 800 and Glasser 360

% Path to subject directories
PathToData = '/Volumes/HCP-Data/HCP_100unrelated_preprocessed_ERG/data';

% Selects the 100 subject folders
a = dir;
a = a(3:end);

% We want to read each subject sequentially...
for s = 1:length(a)
    
    s
    
    % Moving to the directory
    cd(a(s).name);
    
    % Sampling the different types of data: first, resting-state
    cd('rfMRI_REST1_LR');
    
    % Scrubbing data is also loaded
    load('scrubbing.mat');
    Scrubbing_RS{s}(:,1) = scrubbing;
    
    % Z-scored data with global signal regression
    cd('Glasser360');
    load('TS_Glasser360S_gsr_z.mat');
    TS_RS_GSR{s}(:,:,1) = double(TS);
    
    % Same without global signal regression
    load('TS_Glasser360S_z.mat');
    TS_RS_NoGSR{s}(:,:,1) = double(TS);
    
    cd('..');
    cd('..');
    
    % The same is run on all the other types of tasks
    cd('rfMRI_REST1_RL');
    
    load('scrubbing.mat');
    Scrubbing_RS{s}(:,2) = scrubbing;
    
    cd('Glasser360');
    load('TS_Glasser360S_gsr_z.mat');
    TS_RS_GSR{s}(:,:,2) = double(TS);
    
    load('TS_Glasser360S_z.mat');
    TS_RS_NoGSR{s}(:,:,2) = double(TS);
    
    cd('..');
    cd('..');
    
    cd('rfMRI_REST2_LR');
    
    load('scrubbing.mat');
    Scrubbing_RS{s}(:,3) = scrubbing;
    
    cd('Glasser360');
    load('TS_Glasser360S_gsr_z.mat');
    TS_RS_GSR{s}(:,:,3) = double(TS);
    
    load('TS_Glasser360S_z.mat');
    TS_RS_NoGSR{s}(:,:,3) = double(TS);
    
    cd('..');
    cd('..');
    
    cd('rfMRI_REST2_RL');
    
    load('scrubbing.mat');
    Scrubbing_RS{s}(:,4) = scrubbing;
    
    cd('Glasser360');
    load('TS_Glasser360S_gsr_z.mat');
    TS_RS_GSR{s}(:,:,4) = double(TS);
    
    load('TS_Glasser360S_z.mat');
    TS_RS_NoGSR{s}(:,:,4) = double(TS);
    
    
    cd('..');
    cd('..');
    
    
    % EMOTION
    
    cd('tfMRI_EMOTION_LR');
    
    load('scrubbing.mat');
    Scrubbing_EMOTION{s}(:,1) = scrubbing;
    
    cd('Glasser360');
    load('TS_Glasser360S_gsr_z.mat');
    TS_EMOTION_GSR{s}(:,:,1) = double(TS);
    
    load('TS_Glasser360S_z.mat');
    TS_EMOTION_NoGSR{s}(:,:,1) = double(TS);
    
    
    cd('..');
    cd('..');
    
    cd('tfMRI_EMOTION_RL');
    
    load('scrubbing.mat');
    Scrubbing_EMOTION{s}(:,2) = scrubbing;
    
    cd('Glasser360');
    load('TS_Glasser360S_gsr_z.mat');
    TS_EMOTION_GSR{s}(:,:,2) = double(TS);
    
    load('TS_Glasser360S_z.mat');
    TS_EMOTION_NoGSR{s}(:,:,2) = double(TS);
    
    
    cd('..');
    cd('..');
    
    % GAMBLING
    
    cd('tfMRI_GAMBLING_LR');
    
    load('scrubbing.mat');
    Scrubbing_GAMBLING{s}(:,1) = scrubbing;
    
    cd('Glasser360');
    load('TS_Glasser360S_gsr_z.mat');
    TS_GAMBLING_GSR{s}(:,:,1) = double(TS);
    
    load('TS_Glasser360S_z.mat');
    TS_GAMBLING_NoGSR{s}(:,:,1) = double(TS);
    
    
    cd('..');
    cd('..');
    
    cd('tfMRI_GAMBLING_RL');
    
    load('scrubbing.mat');
    Scrubbing_GAMBLING{s}(:,2) = scrubbing;
    
    cd('Glasser360');
    load('TS_Glasser360S_gsr_z.mat');
    TS_GAMBLING_GSR{s}(:,:,2) = double(TS);
    
    load('TS_Glasser360S_z.mat');
    TS_GAMBLING_NoGSR{s}(:,:,2) = double(TS);
    
    
    cd('..');
    cd('..');
    
    % LANGUAGE
    
    cd('tfMRI_LANGUAGE_LR');
    
    load('scrubbing.mat');
    Scrubbing_LANGUAGE{s}(:,1) = scrubbing;
    
    cd('Glasser360');
    load('TS_Glasser360S_gsr_z.mat');
    TS_LANGUAGE_GSR{s}(:,:,1) = double(TS);
    
    load('TS_Glasser360S_z.mat');
    TS_LANGUAGE_NoGSR{s}(:,:,1) = double(TS);
    
    
    cd('..');
    cd('..');
    
    cd('tfMRI_LANGUAGE_RL');
    
    load('scrubbing.mat');
    Scrubbing_LANGUAGE{s}(:,2) = scrubbing;
    
    cd('Glasser360');
    load('TS_Glasser360S_gsr_z.mat');
    TS_LANGUAGE_GSR{s}(:,:,2) = double(TS);
    
    load('TS_Glasser360S_z.mat');
    TS_LANGUAGE_NoGSR{s}(:,:,2) = double(TS);
    
    
    cd('..');
    cd('..');
    
    % MOTOR
    
    cd('tfMRI_MOTOR_LR');
    
    load('scrubbing.mat');
    Scrubbing_MOTOR{s}(:,1) = scrubbing;
    
    cd('Glasser360');
    load('TS_Glasser360S_gsr_z.mat');
    TS_MOTOR_GSR{s}(:,:,1) = double(TS);
    
    load('TS_Glasser360S_z.mat');
    TS_MOTOR_NoGSR{s}(:,:,1) = double(TS);
    
    
    cd('..');
    cd('..');
    
    cd('tfMRI_MOTOR_RL');
    
    load('scrubbing.mat');
    Scrubbing_MOTOR{s}(:,2) = scrubbing;
    
    cd('Glasser360');
    load('TS_Glasser360S_gsr_z.mat');
    TS_MOTOR_GSR{s}(:,:,2) = double(TS);
    
    load('TS_Glasser360S_z.mat');
    TS_MOTOR_NoGSR{s}(:,:,2) = double(TS);
    
    
    cd('..');
    cd('..');
    
    
    % RELATIONAL
    
    cd('tfMRI_RELATIONAL_LR');
    
    load('scrubbing.mat');
    Scrubbing_RELATIONAL{s}(:,1) = scrubbing;
    
    cd('Glasser360');
    load('TS_Glasser360S_gsr_z.mat');
    TS_RELATIONAL_GSR{s}(:,:,1) = double(TS);
    
    load('TS_Glasser360S_z.mat');
    TS_RELATIONAL_NoGSR{s}(:,:,1) = double(TS);
    
    
    cd('..');
    cd('..');
    
    cd('tfMRI_RELATIONAL_RL');
    
    load('scrubbing.mat');
    Scrubbing_RELATIONAL{s}(:,2) = scrubbing;
    
    cd('Glasser360');
    load('TS_Glasser360S_gsr_z.mat');
    TS_RELATIONAL_GSR{s}(:,:,2) = double(TS);
    
    load('TS_Glasser360S_z.mat');
    TS_RELATIONAL_NoGSR{s}(:,:,2) = double(TS);
    
    
    cd('..');
    cd('..');
    
    % SOCIAL
    
    cd('tfMRI_SOCIAL_LR');
    
    load('scrubbing.mat');
    Scrubbing_SOCIAL{s}(:,1) = scrubbing;
    
    cd('Glasser360');
    load('TS_Glasser360S_gsr_z.mat');
    TS_SOCIAL_GSR{s}(:,:,1) = double(TS);
    
    load('TS_Glasser360S_z.mat');
    TS_SOCIAL_NoGSR{s}(:,:,1) = double(TS);
    
    
    cd('..');
    cd('..');
    
    cd('tfMRI_SOCIAL_RL');
    
    load('scrubbing.mat');
    Scrubbing_SOCIAL{s}(:,2) = scrubbing;
    
    cd('Glasser360');
    load('TS_Glasser360S_gsr_z.mat');
    TS_SOCIAL_GSR{s}(:,:,2) = double(TS);
    
    load('TS_Glasser360S_z.mat');
    TS_SOCIAL_NoGSR{s}(:,:,2) = double(TS);
    
    
    cd('..');
    cd('..');
    
    % WM
    
    cd('tfMRI_WM_LR');
    
    load('scrubbing.mat');
    Scrubbing_WM{s}(:,1) = scrubbing;
    
    cd('Glasser360');
    load('TS_Glasser360S_gsr_z.mat');
    TS_WM_GSR{s}(:,:,1) = double(TS);
    
    load('TS_Glasser360S_z.mat');
    TS_WM_NoGSR{s}(:,:,1) = double(TS);
    
    
    cd('..');
    cd('..');
    
    cd('tfMRI_WM_RL');
    
    load('scrubbing.mat');
    Scrubbing_WM{s}(:,2) = scrubbing;
    
    cd('Glasser360');
    load('TS_Glasser360S_gsr_z.mat');
    TS_WM_GSR{s}(:,:,2) = double(TS);
    
    load('TS_Glasser360S_z.mat');
    TS_WM_NoGSR{s}(:,:,2) = double(TS);
    
    
    cd('..');
    cd('..');
    
    
    cd('..');
    
end



%% 2. Estimation as a function of sample size
% We want to assess to what extent the estimation of kurtosis changes when
% adding in more and more data. We also perform bootstrapping to quantify
% how the selection of subjects included impacts the estimates. We run this
% for all different task settings, using the same number of samples every
% time (176, the minimum across all tasks). We only consider the
% region-wise kurtosis for this section

%%%%%% CHANGE ALONG WITH INPUTS

% Loads the data
load('Data');

% Number of regions in our atlas
n_regions = 419;


%%%%%% PARAMETERS

% Number of subjects available in total
n_subjects = 100;

% Order of the MAR model to consider
Order = 1;

% Number of runs to select; we use 2 for matching the different conditions
n_runs = 2;

% Number of folds to bootstrap for
n_folds = 100;

% Max number of subjects selected for this part; we do not use all for
% computational reasons
n_max_subjects = 41;

% Do we pick random subjects when concatenating data? Here, we do because
% we run bootstrapping
is_random = 1;

% Number of time points retained
n_TP = 176;


% Running a loop for computations
idx_subj = 1;

% We start with 3 subjects to avoid issues with badly conditioned matrices
% due to too few data points
for s = 3:2:n_max_subjects
    
    s
    
    % This is our loop for the number of folds
    for f = 1:n_folds
    
        % For each task, we first concatenate the data across subjects, and
        % we then compute the kurtosis. We only consider the region-wise
        % estimates
        [Data_RS] = AoT_ConcatenateRuns(TS_RS_GSR,s,n_runs,n_subjects,n_TP,is_random,'None');
        [deltak_RS(idx_subj,f,:),deltam_RS(idx_subj,f)] = AoT_ComputeKurtosis(Data_RS,Order);
        
        [Data_SOCIAL] = AoT_ConcatenateRuns(TS_SOCIAL_GSR,s,n_runs,n_subjects,n_TP,is_random,'None');
        [deltak_SOCIAL(idx_subj,f,:),deltam_SOCIAL(idx_subj,f)] = AoT_ComputeKurtosis(Data_SOCIAL,Order);
        
        [Data_MOTOR] = AoT_ConcatenateRuns(TS_MOTOR_GSR,s,n_runs,n_subjects,n_TP,is_random,'None');
        [deltak_MOTOR(idx_subj,f,:),deltam_MOTOR(idx_subj,f)] = AoT_ComputeKurtosis(Data_MOTOR,Order);
        
        [Data_GAMBLING] = AoT_ConcatenateRuns(TS_GAMBLING_GSR,s,n_runs,n_subjects,n_TP,is_random,'None');
        [deltak_GAMBLING(idx_subj,f,:),deltam_GAMBLING(idx_subj,f)] = AoT_ComputeKurtosis(Data_GAMBLING,Order);
        
        [Data_WM] = AoT_ConcatenateRuns(TS_WM_GSR,s,n_runs,n_subjects,n_TP,is_random,'None');
        [deltak_WM(idx_subj,f,:),deltam_WM(idx_subj,f)] = AoT_ComputeKurtosis(Data_WM,Order);
        
        [Data_RELATIONAL] = AoT_ConcatenateRuns(TS_RELATIONAL_GSR,s,n_runs,n_subjects,n_TP,is_random,'None');
        [deltak_RELATIONAL(idx_subj,f,:),deltam_RELATIONAL(idx_subj,f)] = AoT_ComputeKurtosis(Data_RELATIONAL,Order);
        
        [Data_LANGUAGE] = AoT_ConcatenateRuns(TS_LANGUAGE_GSR,s,n_runs,n_subjects,n_TP,is_random,'None');
        [deltak_LANGUAGE(idx_subj,f,:),deltam_LANGUAGE(idx_subj,f)] = AoT_ComputeKurtosis(Data_LANGUAGE,Order);
        
        [Data_EMOTION] = AoT_ConcatenateRuns(TS_EMOTION_GSR,s,n_runs,n_subjects,n_TP,is_random,'None');
        [deltak_EMOTION(idx_subj,f,:),deltam_EMOTION(idx_subj,f)] = AoT_ComputeKurtosis(Data_EMOTION,Order);
    end
    
    idx_subj = idx_subj + 1;
        
end

% We want to plot our results
deltak_RS_mean = squeeze(mean(deltak_RS,2));
deltak_RS_std = squeeze(std(deltak_RS,[],2));

deltak_SOCIAL_mean = squeeze(mean(deltak_SOCIAL,2));
deltak_SOCIAL_std = squeeze(std(deltak_SOCIAL,[],2));

deltak_MOTOR_mean = squeeze(mean(deltak_MOTOR,2));
deltak_MOTOR_std = squeeze(std(deltak_MOTOR,[],2));

deltak_GAMBLING_mean = squeeze(mean(deltak_GAMBLING,2));
deltak_GAMBLING_std = squeeze(std(deltak_GAMBLING,[],2));

deltak_WM_mean = squeeze(mean(deltak_WM,2));
deltak_WM_std = squeeze(std(deltak_WM,[],2));

deltak_RELATIONAL_mean = squeeze(mean(deltak_RELATIONAL,2));
deltak_RELATIONAL_std = squeeze(std(deltak_RELATIONAL,[],2));

deltak_LANGUAGE_mean = squeeze(mean(deltak_LANGUAGE,2));
deltak_LANGUAGE_std = squeeze(std(deltak_LANGUAGE,[],2));

deltak_EMOTION_mean = squeeze(mean(deltak_EMOTION,2));
deltak_EMOTION_std = squeeze(std(deltak_EMOTION,[],2));


CM = cbrewer('qual','Set1',8);


% Plots the results: we see a difference between task and rest!
figure;
subplot(2,4,1);
plot(3:2:n_max_subjects,deltak_RS_mean,'color',CM(1,:),'LineWidth',0.5);
title('Resting-state');
xlabel('Number of subjects (2 sessions each)');
ylabel('Regional \Delta kurtosis (F - B)');
xlim([3,41]);

subplot(2,4,2);
plot(3:2:n_max_subjects,deltak_SOCIAL_mean,'color',CM(2,:),'LineWidth',0.5);
title('SOCIAL');
xlabel('Number of subjects (2 sessions each)');
ylabel('Regional \Delta kurtosis (F - B)');
xlim([3,41]);

subplot(2,4,3);
plot(3:2:n_max_subjects,deltak_MOTOR_mean,'color',CM(3,:),'LineWidth',0.5);
title('MOTOR');
xlabel('Number of subjects (2 sessions each)');
ylabel('Regional \Delta kurtosis (F - B)');
xlim([3,41]);

subplot(2,4,4);
plot(3:2:n_max_subjects,deltak_GAMBLING_mean,'color',CM(4,:),'LineWidth',0.5);
title('GAMBLING');
xlabel('Number of subjects (2 sessions each)');
ylabel('Regional \Delta kurtosis (F - B)');
xlim([3,41]);

subplot(2,4,5);
plot(3:2:n_max_subjects,deltak_WM_mean,'color',CM(5,:),'LineWidth',0.5);
title('WM');
xlabel('Number of subjects (2 sessions each)');
ylabel('Regional \Delta kurtosis (F - B)');
xlim([3,41]);

subplot(2,4,6);
plot(3:2:n_max_subjects,deltak_RELATIONAL_mean,'color',CM(6,:),'LineWidth',0.5);
title('RELATIONAL');
xlabel('Number of subjects (2 sessions each)');
ylabel('Regional \Delta kurtosis (F - B)');
xlim([3,41]);

subplot(2,4,7);
plot(3:2:n_max_subjects,deltak_LANGUAGE_mean,'color',CM(7,:),'LineWidth',0.5);
title('LANGUAGE');
xlabel('Number of subjects (2 sessions each)');
ylabel('Regional \Delta kurtosis (F - B)');
xlim([3,41]);

subplot(2,4,8);
plot(3:2:n_max_subjects,deltak_EMOTION_mean,'color',CM(8,:),'LineWidth',0.5);
title('EMOTION');
xlabel('Number of subjects (2 sessions each)');
ylabel('Regional \Delta kurtosis (F - B)');
xlim([3,41]);



%% 3. Comparison, for enough data, across computational schemes
% We want to confirm which paradigms have
% the greatest extent of directionality, and to compare different
% approaches: the regional Kurtosis, the multivariate kurtosis, and the use
% of the Kullback Leibler divergence



n_max_subjects = 50;
n_folds = 30;


for f = 1:n_folds
    
    f
    
    % For each task, we first concatenate the data across subjects, and
    % we then compute the kurtosis. We only consider the region-wise
    % estimates
    [Data_RS] = AoT_ConcatenateRuns(TS_RS_GSR,n_max_subjects,n_runs,n_subjects,n_TP,is_random,'None');
    [DK_RS(f,:),DM_RS(f),DKL_RS(f,:)] = AoT_ComputeKurtosis_Full(Data_RS,Order);

    [Data_SOCIAL] = AoT_ConcatenateRuns(TS_SOCIAL_GSR,n_max_subjects,n_runs,n_subjects,n_TP,is_random,'None');
    [DK_SOCIAL(f,:),DM_SOCIAL(f),DKL_SOCIAL(f,:)] = AoT_ComputeKurtosis_Full(Data_SOCIAL,Order);

    [Data_MOTOR] = AoT_ConcatenateRuns(TS_MOTOR_GSR,n_max_subjects,n_runs,n_subjects,n_TP,is_random,'None');
    [DK_MOTOR(f,:),DM_MOTOR(f),DKL_MOTOR(f,:)] = AoT_ComputeKurtosis_Full(Data_MOTOR,Order);

    [Data_GAMBLING] = AoT_ConcatenateRuns(TS_GAMBLING_GSR,n_max_subjects,n_runs,n_subjects,n_TP,is_random,'None');
    [DK_GAMBLING(f,:),DM_GAMBLING(f),DKL_GAMBLING(f,:)] = AoT_ComputeKurtosis_Full(Data_GAMBLING,Order);

    [Data_WM] = AoT_ConcatenateRuns(TS_WM_GSR,n_max_subjects,n_runs,n_subjects,n_TP,is_random,'None');
    [DK_WM(f,:),DM_WM(f),DKL_WM(f,:)] = AoT_ComputeKurtosis_Full(Data_WM,Order);

    [Data_RELATIONAL] = AoT_ConcatenateRuns(TS_RELATIONAL_GSR,n_max_subjects,n_runs,n_subjects,n_TP,is_random,'None');
    [DK_RELATIONAL(f,:),DM_RELATIONAL(f),DKL_RELATIONAL(f,:)] = AoT_ComputeKurtosis_Full(Data_RELATIONAL,Order);

    [Data_LANGUAGE] = AoT_ConcatenateRuns(TS_LANGUAGE_GSR,n_max_subjects,n_runs,n_subjects,n_TP,is_random,'None');
    [DK_LANGUAGE(f,:),DM_LANGUAGE(f),DKL_LANGUAGE(f,:)] = AoT_ComputeKurtosis_Full(Data_LANGUAGE,Order);
    
    [Data_EMOTION] = AoT_ConcatenateRuns(TS_EMOTION_GSR,n_max_subjects,n_runs,n_subjects,n_TP,is_random,'None');
    [DK_EMOTION(f,:),DM_EMOTION(f),DKL_EMOTION(f,:)] = AoT_ComputeKurtosis_Full(Data_EMOTION,Order);
end


% Multivariate Kurtosis
figure;
boxplot([DM_RS',DM_SOCIAL',DM_MOTOR',DM_GAMBLING',DM_WM',DM_RELATIONAL',DM_LANGUAGE',DM_EMOTION'],'plotstyle','compact','colors',CM,'labels',{'RS','SOCIAL','MOTOR','GAMBLING','WM','RELATIONAL','LANGUAGE','EMOTION'});
set(gca,'Box','off');
ylabel('Multivariate kurtosis');

CM2 = cbrewer('qual','Paired',16);

[CORR_K_KL(1)] = AoT_PlotResults_Measures(DK_RS,DKL_RS,repmat(CM(1,:),2,1),n_regions,'RS');
[CORR_K_KL(2)] = AoT_PlotResults_Measures(DK_SOCIAL,DKL_SOCIAL,repmat(CM(2,:),2,1),n_regions,'SOCIAL');
[CORR_K_KL(3)] = AoT_PlotResults_Measures(DK_MOTOR,DKL_MOTOR,repmat(CM(3,:),2,1),n_regions,'MOTOR');
[CORR_K_KL(4)] = AoT_PlotResults_Measures(DK_GAMBLING,DKL_GAMBLING,repmat(CM(4,:),2,1),n_regions,'GAMBLING');
[CORR_K_KL(5)] = AoT_PlotResults_Measures(DK_WM,DKL_WM,repmat(CM(5,:),2,1),n_regions,'WM');
[CORR_K_KL(6)] = AoT_PlotResults_Measures(DK_RELATIONAL,DKL_RELATIONAL,repmat(CM(6,:),2,1),n_regions,'RELATIONAL');
[CORR_K_KL(7)] = AoT_PlotResults_Measures(DK_LANGUAGE,DKL_LANGUAGE,repmat(CM(7,:),2,1),n_regions,'LANGUAGE');
[CORR_K_KL(8)] = AoT_PlotResults_Measures(DK_EMOTION,DKL_EMOTION,repmat(CM(8,:),2,1),n_regions,'EMOTION');



%% 3 bis. Comparison to the TENET framework results?



%% 4. Selection of the influential regions across paradigms
% We want to determine which are the regions for which the arrow of time is
% "significant"

n_max_subjects = 100;

% Number of null realizations to generate
n_nulls = 500;

% We want to always pick the same set of subjects
is_random = 0;

for n = 1:n_nulls
    
    n
    
    % For each task, we first concatenate the data across subjects, and
    % we then compute the kurtosis. We only consider the region-wise
    % estimates
    [Data_RS] = AoT_ConcatenateRuns(TS_RS_GSR,n_max_subjects,n_runs,n_subjects,n_TP,is_random,'Shuffle');
    [DKNull_Shuffle_RS(n,:)] = AoT_ComputeKurtosis(Data_RS,Order);
    
    [Data_RS] = AoT_ConcatenateRuns(TS_RS_GSR,n_max_subjects,n_runs,n_subjects,n_TP,is_random,'AAFT');
    [DKNull_AAFT_RS(n,:)] = AoT_ComputeKurtosis(Data_RS,Order);
         
end

[Data_RS] = AoT_ConcatenateRuns(TS_RS_GSR,n_max_subjects,n_runs,n_subjects,n_TP,is_random,'None');
[DKActual_RS] = AoT_ComputeKurtosis(Data_RS,Order);


% Selects the regions with significant AoT
is_sign_NULL = AoT_Find_Significant_Regions(DKActual_RS,DKNull_Shuffle_RS);
is_sign_BS = AoT_Find_Significant_Regions(DK_RS,[]);



%% 5. Network-level investigations

% Pools together the data from the same network
for net = 1:8
    tmp = DK_MOTOR(:,Network_Labels(:,1)==net);
    tmp = median(tmp);
    DK_Net{net} = tmp(:);
end

figure;
set(gca,'Box','off');
    
for net = 1:8
    
    h1 = bar(net,mean(DK_Net{net}));
    set(h1,'FaceColor',CM2(net,:),'EdgeColor','None');
    hold on;
    h2 = errorbar(net,mean(DK_Net{net}),abs(mean(DK_Net{net})-prctile((DK_Net{net}),2.5)),abs(mean(DK_Net{net})-prctile((DK_Net{net}),97.5)),'Color',CM2(net,:),'LineWidth',0.1,'LineStyle','None','CapSize',realmin);

    alpha = 0.2;   
    % Set transparency (undocumented)
    set([h2.Bar, h2.Line], 'ColorType', 'truecoloralpha', 'ColorData', [h2.Line.ColorData(1:3); 255*alpha])

end


%% 6. Dynamic analysis
% We want a sliding window-based analysis
% First, we will try on only the data from the WM paradigm


% Parameters of the window
W = 35;
Delta = 1;

t_start = 1;
t_end = W;

t_max = 284;

idx_dyn = 1;

while t_start <= t_max-W+1
    
    t_start
    
    Data = [];
    Data2 = [];
    Data3 = [];
    
    for s = 1:n_max_subjects
        Data = [Data,squeeze(TS_WM_GSR{s}(:,t_start:t_end,1))];
        Data2 = [Data2,squeeze(TS_MOTOR_GSR{s}(:,t_start:t_end,1))];
        Data3 = [Data3,squeeze(TS_RS_GSR{s}(:,t_start:t_end,1))];
    end
    
    DK_DYN_WM(:,idx_dyn) = AoT_ComputeKurtosis(Data',Order);
    DK_DYN_MOTOR(:,idx_dyn) = AoT_ComputeKurtosis(Data2',Order);
    DK_DYN_RS(:,idx_dyn) = AoT_ComputeKurtosis(Data3',Order);
    
    DK_dFC_WM(:,idx_dyn) = jUpperTriMatToVec(AoT_Compute_dFC(Data'));
    DK_dFC_WM_regional(:,idx_dyn) = sum(abs(AoT_Compute_dFC(Data')));
    DK_dFC_MOTOR(:,idx_dyn) = jUpperTriMatToVec(AoT_Compute_dFC(Data2'));
    DK_dFC_MOTOR_regional(:,idx_dyn) = sum(abs(AoT_Compute_dFC(Data2')));
    DK_dFC_RS(:,idx_dyn) = jUpperTriMatToVec(AoT_Compute_dFC(Data3'));
    DK_dFC_RS_regional(:,idx_dyn) = sum(abs(AoT_Compute_dFC(Data3')));
    
    idx_dyn = idx_dyn + 1;
    t_start = t_start + Delta;
    t_end = t_end + Delta;
    
end

% Z-scored plots
figure;
imagesc(zscore(DK_DYN_WM,[],2));
colormap(flipud(cbrewer('div','RdBu',1000)));
caxis([-10,10]);

figure;
imagesc(zscore(DK_DYN_MOTOR,[],2));
colormap(flipud(cbrewer('div','RdBu',1000)));
caxis([-10,10]);

figure;
imagesc(zscore(DK_DYN_RS,[],2));
colormap(flipud(cbrewer('div','RdBu',1000)));
caxis([-10,10]);

% Non z-scored plots
figure;
subplot(2,1,1);
imagesc(DK_DYN_WM);
colormap(flipud(cbrewer('div','RdBu',1000)));
caxis([-10,10]);

subplot(2,1,2);
imagesc(DK_dFC_WM_regional);
colormap(flipud(cbrewer('div','RdBu',1000)));
%caxis([-0.8,0.8]);

figure;
subplot(2,1,1);
imagesc(DK_DYN_MOTOR);
colormap(flipud(cbrewer('div','RdBu',1000)));
caxis([-5,5]);

subplot(2,1,2);
imagesc(DK_dFC_MOTOR_regional);
colormap(flipud(cbrewer('div','RdBu',1000)));
%caxis([-0.8,0.8]);

figure;
subplot(2,1,1);
imagesc(DK_DYN_RS);
colormap(flipud(cbrewer('div','RdBu',1000)));
caxis([-10,10]);

subplot(2,1,2);
imagesc(DK_dFC_RS_regional);
colormap(flipud(cbrewer('div','RdBu',1000)));
%caxis([-0.8,0.8]);
