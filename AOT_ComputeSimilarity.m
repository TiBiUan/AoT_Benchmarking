%% This function computes the similarity in regional AoT patterns across
% several possible paradigm choices
%
% Inputs:
% - X and X2 contain the kurtosis and Kullback-Leibler-based data (all
% subcases as individual dimensions x n_regions)
%
% Outputs:
% - Final_Data contains the data after reordering (after concatenating, in
% meaningful fashion, the input dimensions into only one)
% - SimMat is the correlation matrix
function [Final_Data,SimMat] = AOT_ComputeSimilarity(X,X2)

    % Will contain the reordered data
    Data = zeros(128,419);
    Data2 = zeros(128,419);
    
    % Defines the conditions for each of the factors
    runs_vec = [ones(64,1);2*ones(64,1)];
    gsr_vec = [ones(32,1);2*ones(32,1)];
    gsr_vec = repmat(gsr_vec,128/64,1);
    mot_vec = [ones(8,1);2*ones(8,1);3*ones(8,1);4*ones(8,1)];
    mot_vec = repmat(mot_vec,128/32,1);
    sampling_vec = [ones(2,1);2*ones(2,1);3*ones(2,1);4*ones(2,1)];
    sampling_vec = repmat(sampling_vec,128/8,1);
    epochs_vec = [1;2];
    epochs_vec = repmat(epochs_vec,64,1);
    
    % Does the reordering
    for run = 1:2
        for gsr = 1:2
            for motion = 1:4
                for sampling = 1:4
                    for epochs = 1:2
                        idx_OI = find((runs_vec == run) & (gsr_vec == gsr) & (mot_vec == motion) & (sampling_vec == sampling) & (epochs_vec == epochs));
                        Data(idx_OI,:) = X(run,gsr,motion,sampling,epochs,:);
                        Data2(idx_OI,:) = X2(run,gsr,motion,sampling,epochs,:);
                    end
                end
            end
        end
    end
    
    % Merges kurtosis and Kullback-Leibler results
    Final_Data = [Data;Data2];
    
    % Computes correlation
    SimMat = corr(Final_Data');
end

