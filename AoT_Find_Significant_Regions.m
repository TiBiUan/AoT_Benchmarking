%% This function compares actual data to two types of null distributions
% and determines the regions that show significant directionality
% Since we know that we only expect positive-valued statistics, we will use
% a one-tailed statistical assessment
%
% Inputs:
% - Actual is the actual data for a given paradigm (size 1 x R)
% - Null_* are both null data sets (each n_nulls x R)
function [idx_shuffle,idx_AAFT,pval_shuffle,pval_AAFT] = AoT_Find_Significant_Regions(Actual,Null_Shuffle,Null_AAFT)

    % This is the number of regions at hand, i.e., how many tests we have
    % to correct for in order to carry out Bonferroni correction
    n_regions = length(Actual);
    n_nulls = size(Null_Shuffle,1);

    % We also enable a similar type of assessment for bootstrapped data if
    % the first argument is empty
    if isempty(Actual)
        MIN = prctile(Actual,5);
        
        % In this case, only output argument 1 will be meaningful
        idx_shuffle = find(MIN > 0);
        idx_AAFT = [];
        pval_shuffle = [];
        pval_AAFT = [];
        
    % Null-based assessments
    else
        for r = 1:n_regions
            pval_shuffle(r) = n_regions*sum(Null_Shuffle(:,r) > Actual(r))/n_nulls;
            pval_AAFT(r) = n_regions*sum(Null_AAFT(:,r) > Actual(r))/n_nulls;           
        end
        
        idx_shuffle = find(pval_shuffle < 0.05);
        idx_AAFT = find(pval_AAFT < 0.05);
    end
end