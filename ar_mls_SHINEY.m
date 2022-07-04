%% This function performs the estimation of a multivariate autoregressive
% model on the data at hand
%
% Inputs:
% - TS contains the time samples at hand
% - p is the model order
% - n_TP is the number of samples to retain per run
% - SL, TL and RL are the subject, time and run labels
% - seltype is the type of sample selection scheme to use: can be 'Random'
% (samples are randomly selected for each subject within available ones),
% 'Block' (a contiguous block of samples is selected for each subject,
% with random starting point), 'First' (the first available samples are
% selected for each subject), or 'All' (in this case, we instead use all
% the available samples for each subject, even if this may yield biases
% across paradigms, and a final retained subject with less included
% samples)
%
% Outputs:
% - Y, B, Z and E are the model matrices, such that Y = B*Z + E. Y and Z
% contain the regional time points (for order 1, at t and t+1), B contains
% the coefficients of the model, and E the errors
% - n_available_TPs summarizes how many samples per run are available in
% total
function [Y,B,Z,E,n_available_TPs] = ar_mls_SHINEY(TS,p,n_TP,SL,TL,RL,seltype)

    % Identified parameters
    T   = size(TS,1);
    k   = size(TS,2);

    % Rewriting the problem as Y=BZ+E
    Y   = zeros(k,T-p);
    Z   = zeros(k*p+1,T-p);

    % Filling up Y from the data TS
    Y = TS';
    Y(:, 1:p) = [];

    % Filling up Z

    % First row (intercept)
    Z(1,:)=ones(1,T-p);

    % Other rows
    for j = 1:p
        Z((j-1)*k+2 : j*k+1, :) = TS(p-j+1:T-j, :)';
    end

    % This will contain the labels for the subjects (NSL) and runs (RSL) across
    % all the data that must be kept
    NSL = [];
    RSL = [];

    % Number of distinct runs to consider
    n_runs = length(unique(RL));

    % I want to trim the temporal dimensions of the matrices Y and Z, so that I
    % remove the time points for which the labels are not consistent
    for i = 1:T-p

        % The trimming will be different as a function of the order of the
        % model, but in all cases, we want to remove the data points for which
        % we do not have a succession (t and t +1 in same subject)
        switch p

            % Order 1
            case 1

                if SL(i) ~= SL(i+1) || abs(TL(i)-TL(i+1)) ~= 1 || RL(i) ~= RL(i+1)
                    toremove(i) = 1;
                else
                    toremove(i) = 0;
                    NSL = [NSL,SL(i)];
                    RSL = [RSL,RL(i)];
                end

            % Order 2  
            case 2

                if SL(i) ~= SL(i+1) || SL(i) ~= SL(i+2) || SL(i+1) ~= SL(i+2) || abs(TL(i)-TL(i+1)) ~= 1 || abs(TL(i)-TL(i+2)) ~= 2 || abs(TL(i+2)-TL(i+1)) ~= 1 || RL(i) ~= RL(i+1) || RL(i) ~= RL(i+2) || RL(i+1) ~= RL(i+2) 
                    toremove(i) = 1;
                else
                    toremove(i) = 0;
                    NSL = [NSL,SL(i)];
                    RSL = [RSL,RL(i)];
                end

            otherwise

        end

    end

    % We remove the "wrong" pairs of time points
    Z(:,logical(toremove)) = [];
    Y(:,logical(toremove)) = [];

    % Will contain the final set of data points
    Y_new = [];
    Z_new = [];

    % The label values that are retained
    RetSubj = unique(NSL);
    RetRun = unique(RSL);

    % We remove the excessive time points to reach our desired size
    for idx_subj = 1:length(RetSubj)
        for r = 1:n_runs

            % How many data points are available?
            n_available_TPs(idx_subj,r) = sum(NSL==RetSubj(idx_subj) & RSL==RetRun(r));

            % Contain the data (all available)
            tmpZ = Z(:,NSL==RetSubj(idx_subj) & RSL==RetRun(r));
            tmpY = Y(:,NSL==RetSubj(idx_subj) & RSL==RetRun(r));

            % If we have less than the required amount of data points for the
            % estimates, we can only take as many as we have available...
            if n_available_TPs(idx_subj,r) < n_TP
                tmp_TP = n_available_TPs(idx_subj,r);
            else
                tmp_TP = n_TP;
            end

            % How do we pick the amount of time points that we need within the
            % final matrices?
            switch seltype

                % We may randomly pick a subset of them
                case 'Random'
                    idx = randperm(size(tmpZ,2));
                    Y_new = [Y_new,tmpY(:,idx(1:tmp_TP))];
                    Z_new = [Z_new,tmpZ(:,idx(1:tmp_TP))];

                % We may pick a block as contiguous as possible, starting
                % somewhere random
                % NOTE: for now, the blocks are chosen without considering
                % how "continuous" they truly are with regard to scrubbing
                case 'Block'
                    idx = randi([1,size(tmpZ,2)-tmp_TP+1]);
                    Y_new = [Y_new,tmpY(:,idx:idx+tmp_TP-1)];
                    Z_new = [Z_new,tmpZ(:,idx:idx+tmp_TP-1)];

                % We may just pick the first available ones
                case 'First'
                    Y_new = [Y_new,tmpY(:,1:tmp_TP)];
                    Z_new = [Z_new,tmpZ(:,1:tmp_TP)];
                    
                % We may also just keep all available data (in this case,
                % trimming is done below on the full data matrix)
                case 'All'
                    Y_new = [Y_new,tmpY];
                    Z_new = [Z_new,tmpZ];
                otherwise 
            end
        end
    end
    
    % Trimming to keep only the desired number in the 'All' case
    if strcmp(seltype,'All')        
        Y_new = Y_new(:,1:n_TP*length(RetSubj)*n_runs);
        Z_new = Z_new(:,1:n_TP*length(RetSubj)*n_runs);
    end

    % Identifying AR parameters
    B   = (Y_new*Z_new')/(Z_new*Z_new');

    % Computing residuals
    E   = Y_new-B*Z_new;
end