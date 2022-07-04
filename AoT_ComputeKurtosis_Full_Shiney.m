%% This function computes the metrics of interest on the data at hand
% 
% Inputs:
% - tmp_data contains the concaternated data from all the subjects (regions
% x [time x subjects]); only scrubbed samples are not included
% - Order is the order of the multivariate autoregressive model to use
% - n_TP is the number of time points that we wish to reach per run in our
% computations; it will be used to trim tmp_data properly
% - SL, TL and RL are label vectors for subject, time and run
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
% - deltak contains the estimated kurtosis-based statistic
% - deltaKL contains the estimated Kullback-Leibler divergence-based stat
% - FC and EC contain the estimated group-wise functional connectivity and
% effective connectivity (i.e., MAR parameters)
% - n_av_TPs is the number of retained time points for each subject and
% run; this is after having removed both scrubbed and non-contiguous
% samples (e.g., last time point from a subject paired with first from next
% subject)
function [deltak,deltaKL,FC,EC,n_av_TPs] = AoT_ComputeKurtosis_Full_Shiney(tmp_data,Order,n_TP,SL,TL,RL,seltype)

    % Number of regions at hand
    n_regions = size(tmp_data,2);

    try
        % Forward and backward MAR fitting
        [~,B_f,~,E_f,n_av_TPs] = ar_mls_SHINEY(tmp_data,Order,n_TP,SL,TL,RL,seltype);
        [~,~,~,E_b] = ar_mls_SHINEY(tmp_data(end:-1:1,:),Order,n_TP,SL(end:-1:1),TL(end:-1:1),RL(end:-1:1),seltype);

        % For each region, we can compute a univariate statistic
        for r = 1:n_regions

            % Kurtosis
            deltak(r) = (kurtosis(E_f(r,:))-3).^2 - (kurtosis(E_b(r,:))-3).^2;
        end

        % Creation of a standard normal distribution (compared below to the
        % actual error distributions for forward and backward fitting)
        M = 20;
        Nbin = 1000;
        Norm_distr = normpdf(linspace(-M,M,Nbin));
        
        % KB divergence
        for r = 1:n_regions
            
            % The error data is z-scored for each region, then converted 
            % into a distribution
            Ef_temp = zscore(E_f(r,:));
            [Ef_temp2,~] = histcounts(Ef_temp,Nbin,'BinLimits',[-M,M]);
            
            % Same process for the backward fitting case
            Eb_temp = zscore(E_b(r,:));
            [Eb_temp2,~] = histcounts(Eb_temp,Nbin,'BinLimits',[-M,M]);
            
            % KL divergence is computed
            deltaKL(r) = KLDiv(Ef_temp2,Norm_distr)-KLDiv(Eb_temp2,Norm_distr);

        end
        
        % Computation of FC and EC
        FC = jUpperTriMatToVec(corr(tmp_data));
        EC = B_f(:,2:end);
        EC = EC(:);
        
    catch
        
    end
end