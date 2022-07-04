%% This function creates "no rest" versions of the HCP task paradigms, which
% are also convoluted with the canonical HRF
% 
% Inputs:
% - TS contains the time courses
% - Scrubbing contains the scrubbing information
% - Paradigm contains the paradigm information
% - cut is the final sample to cut at in order to subdivide the
% "Paradigm" array into its two halves
% - rem1 and rem2 are the number of samples to remove at the end of each
% half to have a balanced number of data points (manually changes with each
% paradigm)
% - n_subjects is the total number of subjects
% - TR is the TR of the data (in s)
function [TS_NR,TS_NR2,Scrubbing_NR,P1_smooth,P2_smooth] = AoT_Create_NoRest_Run(TS,TS2,Scrubbing,Paradigm,cut,rem1,rem2,n_subjects,TR,is_plot)

    % Raw paradigm time courses for both halves ("no rest" is set to 1, we do
    % not consider the different sub-conditions separately owing to the limited
    % amount of samples)
    Paradigm1 = Paradigm(1:cut);
    Paradigm1_full = Paradigm1;
    Paradigm1(Paradigm1>0) = 1;
    Paradigm2 = Paradigm(cut+1:end);
    Paradigm2_full = Paradigm2;
    Paradigm2(Paradigm2>0) = 1;
    
    % At this stage, each of the paradigm vectors is filled with only 0
    % (rest moment) and 1 (task moment)
    
    % Creates the HRF to use
    A.dt = TR;
    A.name = 'hrf';
    HRF = spm_get_bf(A);

    % Creating paradigms by convolving with the HRF; we set the frames to
    % retain as the ones with signal > 0.5
    U.u = Paradigm1';
    U.name = {'Whatever'};
    P1 = spm_Volterra(U, HRF.bf);
    P1_smooth = P1;
    P1(P1>0.5) = 1;
    P1(P1<=0.5) = 0;
    P1 = logical(P1);
    
    U.u = Paradigm2';
    P2 = spm_Volterra(U, HRF.bf);
    P2_smooth = P2;
    P2(P2>0.5) = 1;
    P2(P2<=0.5) = 0;
    P2 = logical(P2);
    
    if is_plot
        figure;
        set(gca,'Box','off');
        plot(Paradigm1,'k');
        hold on;
        plot(P1_smooth,'b');
        plot(Paradigm1_full,'m');
        plot(P1,'r');

        figure;
        set(gca,'Box','off');
        plot(Paradigm2,'k');
        hold on;
        plot(P2_smooth,'b');
        plot(Paradigm2_full,'m');
        plot(P2,'r');
    end
    
    % Now goes to use the paradigm information for each subject...
    for s = 1:n_subjects
        
        % Data from first run
        tmp = squeeze(TS{s}(:,:,1));
        tmpp = squeeze(TS2{s}(:,:,1));
        
        % Samples only the frames to retain
        tmp2 = tmp(:,P1);
        tmp22 = tmpp(:,P1);
        
        % Creates the output; to have a matching length between both runs,
        % we remove rem1 frames at the end if needed
        TS_NR{s}(:,:,1) = tmp2(:,1:end-rem1);
        TS_NR2{s}(:,:,1) = tmp22(:,1:end-rem1);

        % Same for the second run
        tmp = squeeze(TS{s}(:,:,2));
        tmp2 = tmp(:,P2);
        TS_NR{s}(:,:,2) = tmp2(:,1:end-rem2);
        
        tmpp = squeeze(TS2{s}(:,:,2));
        tmp22 = tmpp(:,P2);
        TS_NR2{s}(:,:,2) = tmp22(:,1:end-rem2);
        
        % Same process with scrubbing
        tmp = Scrubbing{s}(P1,1);
        Scrubbing_NR{s}(:,1) = tmp(1:end-rem1);

        tmp = Scrubbing{s}(P2,2);
        Scrubbing_NR{s}(:,2) = tmp(1:end-rem2);
    end
end