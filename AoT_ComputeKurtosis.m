function [deltak,deltam] = AoT_ComputeKurtosis(tmp_data,Order)

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
        % deltam = (compute_mardia_coef(E_f')-(n_regions*(n_regions+2)))^2-(compute_mardia_coef(E_b')-(n_regions*(n_regions+2)))^2;
    catch
        1
        deltak = NaN(n_regions,1);
        deltam = NaN;
    end
end