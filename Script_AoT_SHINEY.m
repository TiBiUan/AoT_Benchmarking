%% This script contains all the steps of analyses pertaining to the
% quantification of the arrow-of-time (AoT), a project done in 
% collaboration with the one and only Raphael Liegeois!
% We consider 8 different paradigms: resting state as well as the 7 tasks
% from the Human Connectome Project dataset
%
% The content includes:
% 1. Assessment of the impact of the number of samples on the results: do
% we observe an arrow of time? Is it significantly different from zero? 
% Does it differ across paradigms?
% 2. Generation of surrogate data (amplitude-adjusted phase randomization) 
% for each paradigm, for which causality has been destroyed
% 3. Assessment of regional AoT patterns in each paradigm, comparison (when
% applicable) to "no rest epoch" task paradigms, and quantification of
% network-wise involvement
% 4. Dynamic tracking of AoT changes in time-locked task paradigms across
% subjects



%% 1. Definition of parameters and loading of the data

% Loads the data that we use (Human Connectome Project, 100 subjects)
load('Data_SCH400');
load('task_paradigm.mat');

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

% Colormaps to use
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



%% 2. Generation of the "no rest epochs" task data
% Creates all the "no rest" paradigms, taking the HRF into account

% Working memory (nothing special to do)
[TS_WM_NoRest_NoGSR,TS_WM_NoRest_GSR,Scrubbing_WM_NoRest,Paradigm1_WM,Paradigm2_WM] = ...
    AoT_Create_NoRest_Run(TS_WM_NoGSR,TS_WM_GSR,Scrubbing_WM,task_paradigm.WM,n_samples_WM,0,0,n_subjects,TR,is_plot);

% Motor task (we need to remove one sample from the second run to have 
% matched numbers across runs)
[TS_MOTOR_NoRest_NoGSR,TS_MOTOR_NoRest_GSR,Scrubbing_MOTOR_NoRest,Paradigm1_MOTOR,Paradigm2_MOTOR] = ...
    AoT_Create_NoRest_Run(TS_MOTOR_NoGSR,TS_MOTOR_GSR,Scrubbing_MOTOR,task_paradigm.Motor,n_samples_MOTOR,0,1,n_subjects,TR,is_plot);

% Social task (we remove no sample)
[TS_SOCIAL_NoRest_NoGSR,TS_SOCIAL_NoRest_GSR,Scrubbing_SOCIAL_NoRest,Paradigm1_SOCIAL,Paradigm2_SOCIAL] = ...
    AoT_Create_NoRest_Run(TS_SOCIAL_NoGSR,TS_SOCIAL_GSR,Scrubbing_SOCIAL,task_paradigm.Social,n_samples_SOCIAL,0,0,n_subjects,TR,is_plot);

% Language task (we remove no sample)
[TS_LANGUAGE_NoRest_NoGSR,TS_LANGUAGE_NoRest_GSR,Scrubbing_LANGUAGE_NoRest,Paradigm1_LANGUAGE,Paradigm2_LANGUAGE] = ...
    AoT_Create_NoRest_Run(TS_LANGUAGE_NoGSR,TS_LANGUAGE_GSR,Scrubbing_LANGUAGE,task_paradigm.Language,n_samples_LANGUAGE,0,0,n_subjects,TR,is_plot);

% Emotion task (no sample removed)
[TS_EMOTION_NoRest_NoGSR,TS_EMOTION_NoRest_GSR,Scrubbing_EMOTION_NoRest,Paradigm1_EMOTION,Paradigm2_EMOTION] = ...
    AoT_Create_NoRest_Run(TS_EMOTION_NoGSR,TS_EMOTION_GSR,Scrubbing_EMOTION,task_paradigm.Emotion,n_samples_EMOTION,0,0,n_subjects,TR,is_plot);

% Relational task (4 samples removed in first recording)
[TS_RELATIONAL_NoRest_NoGSR,TS_RELATIONAL_NoRest_GSR,Scrubbing_RELATIONAL_NoRest,Paradigm1_RELATIONAL,Paradigm2_RELATIONAL] = ...
    AoT_Create_NoRest_Run(TS_RELATIONAL_NoGSR,TS_RELATIONAL_GSR,Scrubbing_RELATIONAL,task_paradigm.Relational,n_samples_RELATIONAL,4,0,n_subjects,TR,is_plot);



%% 3. What is the optimal number of samples?
% In this section, we want to determine the optimal number of data samples
% to use for proper estimation, using regional delta kurtosis (Delta_K) as
% metric of interest
% We perform bootstrapping in order to assess how much the estimates
% fluctuate when different subsets of time points/subjects are used
% We perform this process for all available paradigms

% Numbers of samples that we wish to consider (in total) for estimation
n_samples_tot = [1000:500:10000];

% How many runs do we wish to concatenate in total (max amount
% available)
n_runs_toreach = length(run_ids)*n_subjects;

% This is our loop for the number of folds
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

    % Emotion task
    [Data_EMOTION,SL_EMOTION,TL_EMOTION,RL_EMOTION] = ...
    AoT_ConcatenateRuns(TS_EMOTION_GSR,run_ids,n_subjects,100,...
    is_random,ridx(f,:),'None',n_regions,Scrubbing_EMOTION,...
    'None',n_runs_toreach,[],[]);

    % Motor task
    [Data_MOTOR,SL_MOTOR,TL_MOTOR,RL_MOTOR] = ...
    AoT_ConcatenateRuns(TS_MOTOR_GSR,run_ids,n_subjects,100,...
    is_random,ridx(f,:),'None',n_regions,Scrubbing_MOTOR,...
    'None',n_runs_toreach,[],[]);

    % Social task
    [Data_SOCIAL,SL_SOCIAL,TL_SOCIAL,RL_SOCIAL] = ...
    AoT_ConcatenateRuns(TS_SOCIAL_GSR,run_ids,n_subjects,100,...
    is_random,ridx(f,:),'None',n_regions,Scrubbing_SOCIAL,...
    'None',n_runs_toreach,[],[]);

    % Working memory task
    [Data_WM,SL_WM,TL_WM,RL_WM] = ...
    AoT_ConcatenateRuns(TS_WM_GSR,run_ids,n_subjects,100,...
    is_random,ridx(f,:),'None',n_regions,Scrubbing_WM,...
    'None',n_runs_toreach,[],[]);

    % Language task
    [Data_LANGUAGE,SL_LANGUAGE,TL_LANGUAGE,RL_LANGUAGE] = ...
    AoT_ConcatenateRuns(TS_LANGUAGE_GSR,run_ids,n_subjects,100,...
    is_random,ridx(f,:),'None',n_regions,Scrubbing_LANGUAGE,...
    'None',n_runs_toreach,[],[]);

    % Gambling task
    [Data_GAMBLING,SL_GAMBLINGE,TL_GAMBLING,RL_GAMBLING] = ...
    AoT_ConcatenateRuns(TS_GAMBLING_GSR,run_ids,n_subjects,100,...
    is_random,ridx(f,:),'None',n_regions,Scrubbing_GAMBLING,...
    'None',n_runs_toreach,[],[]);

    % Relational task
    [Data_RELATIONAL,SL_RELATIONALE,TL_RELATIONAL,RL_RELATIONAL] = ...
    AoT_ConcatenateRuns(TS_RELATIONAL_GSR,run_ids,n_subjects,100,...
    is_random,ridx(f,:),'None',n_regions,Scrubbing_RELATIONAL,...
    'None',n_runs_toreach,[],[]);
 
    % Running a loop for computations across "number of samples" values
    idx_subj = 1;
    
    % We consider an increasing amount of samples per run; note that n_TP
    % is an equivalent number of samples per run (owing to the way the
    % functions are implemented, this is what we feed as input argument)
    for n_TP = n_samples_tot/n_subjects/length(run_ids)
    
        n_TP
    
        % The metrics of interest are
        % computed for each paradigm; in addition to Delta kurtosis (DK*),
        % we also compute an alternative measure (DKL*), functional
        % connectivity (FC*) and effective connectivity as inferred through
        % MAR-1 coefficients (EC*)
        [DK_RS(idx_subj,f,:),DKL_RS(idx_subj,f,:),FC_RS(idx_subj,f,:),...
            EC_RS(idx_subj,f,:),NATP_RS{idx_subj,f}] = ...
            AoT_ComputeKurtosis_Full_Shiney(Data_RS,Order,n_TP,...
                    SL_RS,TL_RS,RL_RS,'All');
                
        [DK_EMOTION(idx_subj,f,:),DKL_EMOTION(idx_subj,f,:),FC_EMOTION(idx_subj,f,:),...
            EC_EMOTION(idx_subj,f,:),NATP_EMOTION{idx_subj,f}] = ...
            AoT_ComputeKurtosis_Full_Shiney(Data_EMOTION,Order,n_TP,...
                    SL_EMOTION,TL_EMOTION,RL_EMOTION,'All');
                
        [DK_SOCIAL(idx_subj,f,:),DKL_SOCIAL(idx_subj,f,:),FC_SOCIAL(idx_subj,f,:),...
            EC_SOCIAL(idx_subj,f,:),NATP_SOCIAL{idx_subj,f}] = ...
            AoT_ComputeKurtosis_Full_Shiney(Data_SOCIAL,Order,n_TP,...
                    SL_SOCIAL,TL_SOCIAL,RL_SOCIAL,'All');
                
        [DK_MOTOR(idx_subj,f,:),DKL_MOTOR(idx_subj,f,:),FC_MOTOR(idx_subj,f,:),...
            EC_MOTOR(idx_subj,f,:),NATP_MOTOR{idx_subj,f}] = ...
            AoT_ComputeKurtosis_Full_Shiney(Data_MOTOR,Order,n_TP,...
                    SL_MOTOR,TL_MOTOR,RL_MOTOR,'All');
                
        [DK_WM(idx_subj,f,:),DKL_WM(idx_subj,f,:),FC_WM(idx_subj,f,:),...
            EC_WM(idx_subj,f,:),NATP_WM{idx_subj,f}] = ...
            AoT_ComputeKurtosis_Full_Shiney(Data_WM,Order,n_TP,...
                    SL_WM,TL_WM,RL_WM,'All');
                
        [DK_LANGUAGE(idx_subj,f,:),DKL_LANGUAGE(idx_subj,f,:),FC_LANGUAGE(idx_subj,f,:),...
            EC_LANGUAGE(idx_subj,f,:),NATP_LANGUAGE{idx_subj,f}] = ...
            AoT_ComputeKurtosis_Full_Shiney(Data_LANGUAGE,Order,n_TP,...
                    SL_LANGUAGE,TL_LANGUAGE,RL_LANGUAGE,'All');
                
        [DK_GAMBLING(idx_subj,f,:),DKL_GAMBLING(idx_subj,f,:),FC_GAMBLING(idx_subj,f,:),...
            EC_GAMBLING(idx_subj,f,:),NATP_GAMBLING{idx_subj,f}] = ...
            AoT_ComputeKurtosis_Full_Shiney(Data_GAMBLING,Order,n_TP,...
                    SL_GAMBLING,TL_GAMBLING,RL_GAMBLING,'All');
                
        [DK_RELATIONAL(idx_subj,f,:),DKL_RELATIONAL(idx_subj,f,:),FC_RELATIONAL(idx_subj,f,:),...
            EC_RELATIONAL(idx_subj,f,:),NATP_RELATIONAL{idx_subj,f}] = ...
            AoT_ComputeKurtosis_Full_Shiney(Data_RELATIONAL,Order,n_TP,...
                    SL_RELATIONAL,TL_RELATIONAL,RL_RELATIONAL,'All');
                
        idx_subj = idx_subj + 1;
    end
end

% The above is used to plot Figure 2
if is_plot
    
    
end

% We define the optimal number of samples from the above results, used in
% subsequent analyses
n_samples_opt = 8000;



%% 4. Generation of null data
% To be able to infer significance of the arrow of time, we also generate
% null data

% We generate n_null sets of surrogates
for n = 1:n_null

    n
    
    % Random parameters for surrogate data generation: this includes the ones
    % for phase randomization as well as for the normal realizations necessary
    % in the amplitude-adjusted phase randomization process. Owing to the
    % distinct sample numbers across paradigms, we need to generate separate
    % parameters
    PHI_RS = rand([(n_samples_RS-2)/2 1]);
    GAUSS_RS = sort(randn(n_samples_RS,n_regions));

    PHI_MOTOR = rand([(n_samples_MOTOR-2)/2 1]);
    GAUSS_MOTOR = sort(randn(n_samples_MOTOR,n_regions));
    
    n_samples_MOTOR2 = 170;
    PHI_MOTOR2 = rand([(n_samples_MOTOR2-2)/2 1]);
    GAUSS_MOTOR2 = sort(randn(n_samples_MOTOR2,n_regions));

    PHI_SOCIAL = rand([(n_samples_SOCIAL-2)/2 1]);
    GAUSS_SOCIAL = sort(randn(n_samples_SOCIAL,n_regions));

    PHI_LANGUAGE = rand([(n_samples_LANGUAGE-2)/2 1]);
    GAUSS_LANGUAGE = sort(randn(n_samples_LANGUAGE,n_regions));

    PHI_GAMBLING = rand([(n_samples_GAMBLING-1)/2 1]);
    GAUSS_GAMBLING = sort(randn(n_samples_GAMBLING,n_regions));

    PHI_RELATIONAL = rand([(n_samples_RELATIONAL-2)/2 1]);
    GAUSS_RELATIONAL = sort(randn(n_samples_RELATIONAL,n_regions));

    PHI_WM = rand([(n_samples_WM-1)/2 1]);
    GAUSS_WM = sort(randn(n_samples_WM,n_regions));

    PHI_EMOTION = rand([(n_samples_EMOTION-2)/2 1]);
    GAUSS_EMOTION = sort(randn(n_samples_EMOTION,n_regions));

    % We loop across folds...
    for f = 1:n_folds

        f

        % We generate the null data first: it will always be done with the
        % same surrogate parameters regardless of the folds: the data from
        % each subject is randomized, then concatenated with a random
        % subject ordering
        [DN_RS,subj_label,time_label,run_label] = ...
        AoT_ConcatenateRuns(TS_RS_GSR,run_ids,n_subjects,n_samples_opt/n_subjects/length(run_ids),...
        is_random,ridx(f,:),'AAFT',n_regions,Scrubbing_RS,...
        'None',n_runs_toreach,PHI_RS,GAUSS_RS);

        % In a subsequent step, the metrics of interest are then
        % computed
        [DKN_RS(f,n,:),DKLN_RS(f,n,:),FCN_RS(f,n,:),...
            ECN_RS(f,n,:),NATPN_RS{f,n}] = ...
            AoT_ComputeKurtosis_Full_Shiney(DN_RS,Order,n_samples_opt/n_subjects/length(run_ids),...
                    subj_label,time_label,run_label,'All');

        % Motor paradigm (full)
        [DN_MOTOR,subj_label,time_label,run_label] = ...
        AoT_ConcatenateRuns(TS_MOTOR_GSR,run_ids,n_subjects,n_samples_opt/n_subjects/length(run_ids),...
        is_random,ridx(f,:),'AAFT',n_regions,Scrubbing_MOTOR,...
        'None',n_runs_toreach,PHI_MOTOR,GAUSS_MOTOR);

        [DKN_MOTOR(f,n,:),DKLN_MOTOR(f,n,:),FCN_MOTOR(f,n,:),...
            ECN_MOTOR(f,n,:),NATPN_MOTOR{f,n}] = ...
            AoT_ComputeKurtosis_Full_Shiney(DN_MOTOR,Order,n_samples_opt/n_subjects/length(run_ids),...
                    subj_label,time_label,run_label,'All');
                
        % Motor paradigm (only no rest bits)
        [DN_MOTOR2,subj_label,time_label,run_label] = ...
        AoT_ConcatenateRuns(TS_MOTOR_NoRest_GSR,run_ids,n_subjects,n_samples_opt/n_subjects/length(run_ids),...
        is_random,ridx(f,:),'AAFT',n_regions,Scrubbing_MOTOR_NoRest,...
        'None',n_runs_toreach,PHI_MOTOR2,GAUSS_MOTOR2);

        [DKN_MOTOR2(f,n,:),DKLN_MOTOR2(f,n,:),FCN_MOTOR2(f,n,:),...
            ECN_MOTOR2(f,n,:),NATPN_MOTOR2{f,n}] = ...
            AoT_ComputeKurtosis_Full_Shiney(DN_MOTOR2,Order,n_samples_opt/n_subjects/length(run_ids),...
                    subj_label,time_label,run_label,'All');

        % Emotion paradigm
        [DN_EMOTION,subj_label,time_label,run_label] = ...
        AoT_ConcatenateRuns(TS_EMOTION_GSR,run_ids,n_subjects,n_samples_opt/n_subjects/length(run_ids),...
        is_random,ridx(f,:),'AAFT',n_regions,Scrubbing_EMOTION,...
        'None',n_runs_toreach,PHI_EMOTION,GAUSS_EMOTION);

        [DKN_EMOTION(f,n,:),DKLN_EMOTION(f,n,:),FCN_EMOTION(f,n,:),...
            ECN_EMOTION(f,n,:),NATPN_EMOTION{f,n}] = ...
            AoT_ComputeKurtosis_Full_Shiney(DN_EMOTION,Order,n_samples_opt/n_subjects/length(run_ids),...
                    subj_label,time_label,run_label,'All');
        
        % Language paradigm
        [DN_LANGUAGE,subj_label,time_label,run_label] = ...
        AoT_ConcatenateRuns(TS_LANGUAGE_GSR,run_ids,n_subjects,n_samples_opt/n_subjects/length(run_ids),...
        is_random,ridx(f,:),'AAFT',n_regions,Scrubbing_LANGUAGE,...
        'None',n_runs_toreach,PHI_LANGUAGE,GAUSS_LANGUAGE);

        [DKN_LANGUAGE(f,n,:),DKLN_LANGUAGE(f,n,:),FCN_LANGUAGE(f,n,:),...
            ECN_LANGUAGE(f,n,:),NATPN_LANGUAGE{f,n}] = ...
            AoT_ComputeKurtosis_Full_Shiney(DN_LANGUAGE,Order,n_samples_opt/n_subjects/length(run_ids),...
                    subj_label,time_label,run_label,'All');

        % Social paradigm
        [DN_SOCIAL,subj_label,time_label,run_label] = ...
        AoT_ConcatenateRuns(TS_SOCIAL_GSR,run_ids,n_subjects,n_samples_opt/n_subjects/length(run_ids),...
        is_random,ridx(f,:),'AAFT',n_regions,Scrubbing_SOCIAL,...
        'None',n_runs_toreach,PHI_SOCIAL,GAUSS_SOCIAL);

        [DKN_SOCIAL(f,n,:),DKLN_SOCIAL(f,n,:),FCN_SOCIAL(f,n,:),...
            ECN_SOCIAL(f,n,:),NATPN_SOCIAL{f,n}] = ...
            AoT_ComputeKurtosis_Full_Shiney(DN_SOCIAL,Order,n_samples_opt/n_subjects/length(run_ids),...
                    subj_label,time_label,run_label,'All');

        % Working memory paradigm
        [DN_WM,subj_label,time_label,run_label] = ...
        AoT_ConcatenateRuns(TS_WM_GSR,run_ids,n_subjects,n_samples_opt/n_subjects/length(run_ids),...
        is_random,ridx(f,:),'AAFT',n_regions,Scrubbing_WM,...
        'None',n_runs_toreach,PHI_WM,GAUSS_WM);

        [DKN_WM(f,n,:),DKLN_WM(f,n,:),FCN_WM(f,n,:),...
            ECN_WM(f,n,:),NATPN_WM{f,n}] = ...
            AoT_ComputeKurtosis_Full_Shiney(DN_WM,Order,n_samples_opt/n_subjects/length(run_ids),...
                    subj_label,time_label,run_label,'All');

        % Relational paradigm
        [DN_RELATIONAL,subj_label,time_label,run_label] = ...
        AoT_ConcatenateRuns(TS_RELATIONAL_GSR,run_ids,n_subjects,n_samples_opt/n_subjects/length(run_ids),...
        is_random,ridx(f,:),'AAFT',n_regions,Scrubbing_RELATIONAL,...
        'None',n_runs_toreach,PHI_RELATIONAL,GAUSS_RELATIONAL);

        [DKN_RELATIONAL(f,n,:),DKLN_RELATIONAL(f,n,:),FCN_RELATIONAL(f,n,:),...
            ECN_RELATIONAL(f,n,:),NATPN_RELATIONAL{f,n}] = ...
            AoT_ComputeKurtosis_Full_Shiney(DN_RELATIONAL,Order,n_samples_opt/n_subjects/length(run_ids),...
                    subj_label,time_label,run_label,'All');

        % Gambling paradigm
        [DN_GAMBLING,subj_label,time_label,run_label] = ...
        AoT_ConcatenateRuns(TS_GAMBLING_GSR,run_ids,n_subjects,n_samples_opt/n_subjects/length(run_ids),...
        is_random,ridx(f,:),'AAFT',n_regions,Scrubbing_GAMBLING,...
        'None',n_runs_toreach,PHI_GAMBLING,GAUSS_GAMBLING);

        [DKN_GAMBLING(f,n,:),DKLN_GAMBLING(f,n,:),FCN_GAMBLING(f,n,:),...
            ECN_GAMBLING(f,n,:),NATPN_GAMBLING{f,n}] = ...
            AoT_ComputeKurtosis_Full_Shiney(DN_GAMBLING,Order,n_samples_opt/n_subjects/length(run_ids),...
                    subj_label,time_label,run_label,'All');
    end
end

% From the regional null data, we can compute a whole-brain null as well.
% For Figure 2 (bottom panel), we aggregate the data across null folds
% together
NullDist_RS = squeeze(median(DKN_RS,1));
NullDist_RS = NullDist_RS(:);



%% 5. Assessment of regional AoT patterns
% Here, we consider regional patterns across tasks, using the null data for
% proper thresholding

% For the resting state case below:

% We need to determine the significance threshold for each region...
for r = 1:n_regions
    
    % Actual values (n_regions x 1)
    tmp_actual = median(squeeze(DK_RS(n_samples_tot == n_samples_opt,:,:)),1)';
    
    % Null values (n_regions x n_null)
    tmp_null = squeeze(median(DKN_RS,1))';
    
    % Estimation of a Gaussian fit and extraction of its parameters
    
    % Extraction of significance thresholds
end



%% 6. Dynamic tracking of the AoT
% We leverage a sliding window approach to track fluctuations of AoT over
% time, and compare them to the ones of activation time courses and of
% dynamic functional connectivity outcomes
 
% Motor paradigm
[DK_dyn_MOTOR,DKL_dyn_MOTOR,FC_dyn_MOTOR,EC_dyn_MOTOR,Act_dyn_MOTOR] = ...
    AoT_Compute_Dynamic_Evolution(TS_MOTOR_GSR,W,Delta,run_ids,Order);

% Emotion paradigm
[DK_dyn_EMOTION,DKL_dyn_EMOTION,FC_dyn_EMOTION,EC_dyn_EMOTION,Act_dyn_EMOTION] = ...
    AoT_Compute_Dynamic_Evolution(TS_EMOTION_GSR,W,Delta,run_ids,Order);

% Working memory paradigm
[DK_dyn_WM,DKL_dyn_WM,FC_dyn_WM,EC_dyn_WM,Act_dyn_WM] = ...
    AoT_Compute_Dynamic_Evolution(TS_WM_GSR,W,Delta,run_ids,Order);

% Social paradigm
[DK_dyn_SOCIAL,DKL_dyn_SOCIAL,FC_dyn_SOCIAL,EC_dyn_SOCIAL,Act_dyn_SOCIAL] = ...
    AoT_Compute_Dynamic_Evolution(TS_SOCIAL_GSR,W,Delta,run_ids,Order);

% Language paradigm
[DK_dyn_LANGUAGE,DKL_dyn_LANGUAGE,FC_dyn_LANGUAGE,EC_dyn_LANGUAGE,Act_dyn_LANGUAGE] = ...
    AoT_Compute_Dynamic_Evolution(TS_LANGUAGE_GSR,W,Delta,run_ids,Order);

% Time frame of interest to zoom on the motor paradigm
time_OI = 47:66;
time_OI2 = time_OI + 83;

% How many regions to examine (we select the top 4 with the largest 
% AoT "energy" over the time interval of interest)
Top_number = 20;

% The time courses of these top regions, as well as their indices and their
% hemisphere, are extracted. Note that we focus on the left hemisphere
% regions for now
[TS_most_causal,idx_most_causal,hemi] = AoT_ExamineTimeInterval(DK_dyn_MOTOR(1:200,:),time_OI,Top_number);
[TS_most_causal2,idx_most_causal2,hemi] = AoT_ExamineTimeInterval(DK_dyn_MOTOR(1:200,:),time_OI2,Top_number);

% Plots for the generation of Figure 4
if is_plot
    
    % Heatmaps for dynamic changes in activation (smoothed, dynamic FC and
    % AoT) for the motor task
    figure;
    imagesc(1:size(Act_dyn_MOTOR,2),1:n_regions,Act_dyn_MOTOR);
    colormap(CM_RB);
    caxis([-max(abs(Act_dyn_MOTOR(:))),max(abs(Act_dyn_MOTOR(:)))]);
    title('Smoothed activation');
    xlabel('Time [window index]');
    ylabel('Regions');
    set(gca,'Box','off');
    
    figure;
    imagesc(1:size(Act_dyn_MOTOR,2),1:n_regions,FC_dyn_MOTOR);
    colormap(CM_RB);
    caxis([-max(abs(FC_dyn_MOTOR(:))),max(abs(FC_dyn_MOTOR(:)))]);
    title('Dynamic regional degree');
    xlabel('Time [window index]');
    ylabel('Regions');
    set(gca,'Box','off');
    
    figure;
    imagesc(1:size(Act_dyn_MOTOR,2),1:n_regions,DK_dyn_MOTOR);
    colormap(CM_RB);
    caxis(0.5*[-max(abs(DK_dyn_MOTOR(:))),max(abs(DK_dyn_MOTOR(:)))]);
    title('Dynamic arrow-of-time');
    xlabel('Time [window index]');
    ylabel('Regions');
    set(gca,'Box','off');
    
    % Generation of the paradigm time course (convoluted with the HRF) for
    % the motor task
    figure;
    imagesc(Paradigm1_MOTOR');
    colormap(flipud(cbrewer('seq','Greys',1000)));
    caxis([-0.6,1.6]);
    set(gca,'Box','off');
    
    % Zoom on the four most influential areas
    figure;
    imagesc(TS_most_causal);
    colormap(CM_RB);
    caxis([-1,1]);
end

