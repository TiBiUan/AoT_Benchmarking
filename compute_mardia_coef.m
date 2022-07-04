function mardia_kurtosis = compute_mardia_coef(X)

% compute the mardia's kurtosis coefficient,
% the exepect mardia's kurtosis is n_var * (n_var + 2) for a normal distribution

[n_sample,n_var] = size(X);

difT = [];

for	j = 1:n_var
   eval(['difT=[difT,(X(:,j)-mean(X(:,j)))];']);
end

S = cov(X,1);

D = difT * inv(S) * difT';  %squared-Mahalanobis' distances matrix
mardia_kurtosis = trace(D.^2)/n_sample;  %multivariate kurtosis coefficient

% if one need the mardia's skewness coefficient:
% mardia_skewness = (sum(sum(D.^3)))/n_sample^2;  %multivariate skewness coefficient

end

