%% This function computes the dynamic evolution of the arrow of time
%
% Inputs:
% - TS contains the time courses
% - W is the window length and Delta the shift
% - run_ids summarizes which runs to consider
% - Order is the order of the model
%
% Outputs:
% - deltak is the dynamic evolution of the AoT (kurtosis estimation)
% - deltakl is the same for the Kullback-Leibler divergence-based approach
% - DEG-FC is the dynamic evolution of functional connectivity, where we
% sum across all incoming edges to get a regional measure
% - EC is the effective connectivity information across connections
% - Act contains the smoothed activation time courses of the regions
function [deltak,deltakl,DEG_FC,EC,Act] = AoT_Compute_Dynamic_Evolution(TS,W,Delta,run_ids,Order)

    % Start and end of the window
    t_start = 1;
    t_end = W;

    % Sliding window indices
    idx_dyn = 1;
    
    % Max time sample at disposal
    t_max = size(TS{1},2);
    
    % Number of regions
    n_regions = size(TS{1},1);

    % Window computations while we can still go further in time...
    while t_start <= t_max-W+1

        t_start

        Data = [];
        subj_label = [];
        run_label = [];
        time_label = [];

        % Concatenating the data across subjects and runs
        for s = 1:length(TS)
            for run = 1:length(run_ids)
                tmp_data = squeeze(TS{s}(:,t_start:t_end,run));
                Data = [Data,tmp_data];
                subj_label = [subj_label,s*ones(1,size(tmp_data,2))];
                run_label = [run_label,run*ones(1,size(tmp_data,2))];
                time_label = [time_label,t_start:t_end];
            end
        end

        % Computation of the metrics of interest
        [deltak(:,idx_dyn),deltakl(:,idx_dyn),FC(:,idx_dyn),...
            EC(:,idx_dyn)] = ...
            AoT_ComputeKurtosis_Full_Shiney(Data',Order,W-1,...
                    subj_label,time_label,run_label,'All');
                
        tmp = jVecToUpperTriMat(FC(:,idx_dyn),n_regions)+jVecToUpperTriMat(FC(:,idx_dyn),n_regions)';
        DEG_FC(:,idx_dyn) = sum((tmp))/2;
                
        Act(:,idx_dyn) = mean(Data,2);

        idx_dyn = idx_dyn + 1;
        t_start = t_start + Delta;
        t_end = t_end + Delta;
    end
end