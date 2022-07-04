%% This function extracts the AoT time courses for the regions that show
% the largest energy (biggest values on the whole) within a time interval
%
% Inputs:
% - TS contains all the AoT time courses (regions x time)
% - time_OI is a vector with the time samples to zoom on
% - top_number is how many regions we want to extract
function [TS_most_causal,idx_most_causal,hemi] = AoT_ExamineTimeInterval(TS,time_OI,top_number)
    
    % Samples the subset of interest
    X_OI = TS(:,time_OI);
    
    % Computes the overall absolute intensity over the time frame
    Intensity = sum(abs(X_OI),2);
    
    % Extracts the indices of the most influential areas
    [~,idx_causal] = sort(Intensity,'descend');

    % Most causal regions and their associated time courses
    idx_most_causal = idx_causal(1:top_number);
    TS_most_causal = X_OI(idx_most_causal,:);

    % Normalizes so that the largest value is equal to +1 or -1
    for r = 1:size(TS_most_causal,1)
        TS_most_causal(r,:) = TS_most_causal(r,:)/max(abs(TS_most_causal(r,:)));
    end
    
    % Determination of the hemisphere
    hemi = idx_most_causal;
    hemi(hemi <= 200) = 1;
    hemi(hemi > 200) = 2;
end