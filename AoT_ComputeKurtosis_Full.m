function [deltak,deltam,deltaKL,E_b,Params_b] = AoT_ComputeKurtosis_Full(tmp_data,Order)

    deltam = NaN;

    n_regions = size(tmp_data,2);

    try
        % Forward and backward MAR fitting
        [~,Params_f,~,E_f] = ar_mls(tmp_data,Order);
        [~,Params_b,~,E_b] = ar_mls(tmp_data(end:-1:1,:),Order);

        % For each region, we can compute a univariate statistic
        for r = 1:n_regions

            % Kurtosis
            deltak(r) = (kurtosis(E_f(r,:))-3).^2 - (kurtosis(E_b(r,:))-3).^2;
        end

        % Multivariate kurtosis
        %deltam = (compute_mardia_coef(E_f')-(n_regions*(n_regions+2)))^2-(compute_mardia_coef(E_b')-(n_regions*(n_regions+2)))^2;
        
        M = 20;
        Nbin = 1000;
        Norm_distr = normpdf(linspace(-M,M,Nbin));
        
        % KB divergence
        for r = 1:n_regions
            Ef_temp = zscore(E_f(r,:)); %normalizing before binning
            [Ef_temp2,~] = histcounts(Ef_temp,Nbin,'BinLimits',[-M,M]);
            Eb_temp = zscore(E_b(r,:)); %normalizing before binning
            [Eb_temp2,~] = histcounts(Eb_temp,Nbin,'BinLimits',[-M,M]);
            deltaKL(r) = KLDiv(Ef_temp2,Norm_distr)-KLDiv(Eb_temp2,Norm_distr);

        end
        
    catch
        1
        deltak = NaN(n_regions,1);
        deltam = NaN;
        deltaKL = NaN(n_regions,1);
    end
end