%% This function runs supplementary analyses for the data of a given
% paradigm. We want to assess similarity in the obtained AoT patterns when
% varying run index, GSR inclusion/exclusion, scrubbing scheme, AoT measure
% used, way to sample the time points to incorporate, and type of recording
% (only task blocks vs full) to consider
function [DK_SM,DKL_SM,NATP_SM,Final_Data,SimMat] = ...
    AOT_ComputeSupplementaryAnalyses_RS(TS_MOTOR_GSR,...
    TS_MOTOR_NoGSR,Scrubbing_MOTOR,...
    n_samples_opt,n_subjects,n_folds,n_regions,...
    n_runs_toreach,is_random,Order)

    % We will start with the motor case. In total, we have 128 different
    % combinations to assess...
    Run_choices = {'LR','RL'};
    GSR_choices = {'With','Without'};
    Motion_choices = {'None','Mild','Moderate','Stringent'};
    Measure_choices = {'Kurtosis','KL divergence'};
    Sampling_choices = {'All','Random','First','Block'};

    % Number of samples we want per subject
    n_samples_needed = n_samples_opt/n_subjects;

    for run = 1:2
        for gsr = 1:2
            for motion = 1:4
                for sampling = 1:4

                    % We run these for n_folds folds...
                    for f = 1:n_folds

                        f

                        if gsr == 1

                            [Data_SM,SL_SM,TL_SM,RL_SM] = ...
                                AoT_ConcatenateRuns(TS_MOTOR_GSR,run,n_subjects,100,...
                                is_random,[],'None',n_regions,Scrubbing_MOTOR,...
                                Motion_choices{motion},n_runs_toreach,[],[]);

                            [DK_SM(f,run,gsr,motion,sampling,:),DKL_SM(f,run,gsr,motion,sampling,:),...
                                NATP_SM{f,run,gsr,motion,sampling}] = ...
                                AoT_ComputeKurtosis_Full_Shiney(Data_SM,Order,n_samples_needed,...
                                        SL_SM,TL_SM,RL_SM,Sampling_choices{sampling});

                        elseif gsr == 2

                            [Data_SM,SL_SM,TL_SM,RL_SM] = ...
                                AoT_ConcatenateRuns(TS_MOTOR_NoGSR,run,n_subjects,100,...
                                is_random,[],'None',n_regions,Scrubbing_MOTOR,...
                                Motion_choices{motion},n_runs_toreach,[],[]);

                            [DK_SM(f,run,gsr,motion,sampling,:),DKL_SM(f,run,gsr,motion,sampling,:),...
                                NATP_SM{f,run,gsr,motion,sampling}] = ...
                                AoT_ComputeKurtosis_Full_Shiney(Data_SM,Order,n_samples_needed,...
                                        SL_SM,TL_SM,RL_SM,Sampling_choices{sampling});
                        end
                    end
                end
            end
        end
    end

    [Final_Data,SimMat] = AOT_ComputeSimilarity_RS(squeeze(median(DK_SM,1)),squeeze(median(DKL_SM,1)));
end