function [deltak , deltam , deltaKL, Ef , Eb] = get_arrow(TS)

%getting delta kurtosis from multivariate time series TS of size T(nb. time
%points) x N(nb. variables)

T = size(TS,1);
N = size(TS,2);

[~,~,~,Ef] = ar_mls(TS,1);
[~,~,~,Eb] = ar_mls(TS(end:-1:1,:),1);

deltak = zeros(1,N);
deltaKL = zeros(1,N);
Nbin = 100;
M = 20; %Max and min values of the distrubutions to be compared to
Norm_distr = normpdf(linspace(-M,M,Nbin));
for i=1:N
    deltak(i) = (kurtosis(Ef(i,:))-3)^2-(kurtosis(Eb(i,:))-3)^2; %expected to be positive
    Ef_temp = zscore(Ef(i,:)); %normalizing before binning
    [Ef_temp2,~] = histcounts(Ef_temp,Nbin,'BinLimits',[-M,M]);
    Eb_temp = zscore(Eb(i,:)); %normalizing before binning
    [Eb_temp2,~] = histcounts(Eb_temp,Nbin,'BinLimits',[-M,M]);
    deltaKL(i) = KLDiv(Ef_temp2,Norm_distr)-KLDiv(Eb_temp2,Norm_distr);
end

deltam = (compute_mardia_coef(Ef')-(N*(N+2)))^2-(compute_mardia_coef(Eb')-(N*(N+2)))^2;
