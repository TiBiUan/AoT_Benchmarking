%% This function concatenates the data from different HCP subjects together
% prior to computations
%
% Inputs:
% - TC contains the time courses for all the subjects (each cell has size
% R x T x n_runs, with R the number of regions, T of time points)
% - runs is a vector containing the indices of the runs to concatenate
% - n_tot_subjects is the total number of subjects to our disposal
% - n_TP is the number of time points to sample per run per subject
% - is_random defines whether we want to select the subjects randomly or in
% order
% - rand_sequence is a provided order of subject indices that we want to
% sample (if we wish to use the same random order across several calls, we
% want to pregenerate it outside of the function); keep empty if there is
% no predefined order
% - null_type depicts whether we want to distort the data according to a
% given null scheme. Available options are 'None', 'Shuffle' (for shuffling
% the time points), 'PR' (phase randomisation) or 'AAFT' 
% (amplitude-adjusted phase randomisation)
% - region_end is the index of the last region to retain; this is because
% sometimes, we want to discard subcortical areas (the final indices). Use
% R if all regions should be kept
% - scrubbing contains the scrubbing information as a cell array (one cell
% per subject), with each cell of size T x n_runs
% - scrubbing_level defines the level of rigor of scrubbing; it can be
% 'none' (no scrubbing), 'mild' (only scrubbed time points are removed), 
% 'moderate' (1 time point before and 2 after) or 'rigorous' (3 time
% points before and 6 time points after)
% - n_runs_toreach is the number of runs that we want to extract
% information from
%
% Outputs:
% - Data contains the output concatenated data (n_subjects*max_run*n_TP x 
% n_regions)
% - subj_label contains the index of the subject that each retained sample
% belonged to. The time_label and run_label vectors contain the same
% information for time point and run index
% - n_removed_runs is the number of runs that we could not properly process
% across the analyzed ones (because not enough data points were available
% upon scrubbing)
% - removed_DPs summarises the removed data points across time, runs and 
% subjects
% - n_added_runs is the final number of runs that were concatenated with
% success
function [Data,subj_label,time_label,run_label,n_removed_runs,removed_DPs,n_added_runs] = ...
    AoT_ConcatenateRuns(TC,runs,n_tot_subjects,n_TP,...
    is_random,rand_sequence,null_type,region_end,scrubbing,...
    scrubbing_level,n_runs_toreach,PHI,GAUSS)
    
    % Subjects to concatenate, either in order or randomly picked, are in
    % subjects_OI
    if is_random
        if isempty(rand_sequence)
            subjects_OI = randperm(n_tot_subjects);
        else
            subjects_OI = rand_sequence;
        end
    else
        subjects_OI = 1:n_tot_subjects;
    end
    
    % Ensures that there is no wrong combination of input parameters
    if (strcmp(scrubbing_level,'Mild') || strcmp(scrubbing_level,'Medium') || strcmp(scrubbing_level,'Rigorous')) && ~strcmp(null_type,'None')
        errordlg('Invalid parameter combination: scrubbing cannot be done together with surrogate generation!');
    end
        
        
    % This will contain our outputs; the data will be added progressively
    Data = [];
    subj_label = [];
    time_label = [];
    run_label = [];
    
    % This indexes how many subjects have been considered so far
    idx_subject = 1;
    
    % This indexes the number of removed runs in total
    n_removed_runs = 0;
    n_added_runs = 0;
    
    % We perform the same process for each of the available runs
    while n_added_runs < n_runs_toreach
        for run = runs  
            
            % Reading the data for the run in question and the subject at
            % hand (both time courses and scrubbing)
            tmp_data = squeeze(TC{subjects_OI(idx_subject)}(1:region_end,:,run));
            tmp_scrub = logical(scrubbing{idx_subject}(:,run))';
            
            % This contains the current time labels before removing
            % anything
            tmp_time_label = 1:size(tmp_data,2);
            
            % Creates the scrubbing mask, at a predefined level of
            % stringency
            switch scrubbing_level
                
                % Mild: only removing scrubbed time points
                case 'Mild'
                    tmp_scrub_final = tmp_scrub;
                    
                % Medium: also removing two time points later and one
                % beforehand
                case 'Medium'
                    
                    tmp_scrub1 = circshift(tmp_scrub,1);
                    tmp_scrub1 = [0,tmp_scrub1(2:end)];
                    tmp_scrub2 = circshift(tmp_scrub,2);
                    tmp_scrub2 = [0,0,tmp_scrub2(3:end)];
                    
                    tmp_scrubm1 = circshift(tmp_scrub,-1);
                    tmp_scrubm1 = [tmp_scrubm1(1:end-1),0];
                    
                    % Tags all the time points to scrub out
                    tmp_scrub_final = tmp_scrub | tmp_scrub1 | ...
                        tmp_scrub2 | tmp_scrubm1;
                    
                    % This removes "single time points"
                    %tmp_scrub_final = (~(bwareaopen(~tmp_scrub_final,2)));
                    
                % Rigorous: removing 6 time points before and 3 after
                case 'Rigorous'
                    tmp_scrub1 = circshift(tmp_scrub,1);
                    tmp_scrub1 = [0,tmp_scrub1(2:end)];
                    tmp_scrub2 = circshift(tmp_scrub,2);
                    tmp_scrub2 = [0,0,tmp_scrub2(3:end)];
                    tmp_scrub3 = circshift(tmp_scrub,3);
                    tmp_scrub3 = [0,0,0,tmp_scrub3(4:end)];
                    tmp_scrub4 = circshift(tmp_scrub,4);
                    tmp_scrub4 = [0,0,0,0,tmp_scrub4(5:end)];
                    tmp_scrub5 = circshift(tmp_scrub,5);
                    tmp_scrub5 = [0,0,0,0,0,tmp_scrub5(6:end)];
                    tmp_scrub6 = circshift(tmp_scrub,6);
                    tmp_scrub6 = [0,0,0,0,0,0,tmp_scrub6(7:end)];
                    tmp_scrubm1 = circshift(tmp_scrub,-1);
                    tmp_scrubm1 = [tmp_scrubm1(1:end-1),0];
                    tmp_scrubm2 = circshift(tmp_scrub,-2);
                    tmp_scrubm2 = [tmp_scrubm2(1:end-2),0,0];
                    tmp_scrubm3 = circshift(tmp_scrub,-3);
                    tmp_scrubm3 = [tmp_scrubm3(1:end-3),0,0,0];
                    
                    tmp_scrub_final = tmp_scrub | tmp_scrub1 | ...
                        tmp_scrub2 | tmp_scrub3 | tmp_scrub4 | ...
                        tmp_scrub5 | tmp_scrub6 | tmp_scrubm1 | ...
                        tmp_scrubm2 | tmp_scrubm3;
                    
                    %tmp_scrub_final = (~(bwareaopen(~tmp_scrub_final,2)));
                    
                % Else, we just do not remove anything!
                otherwise
                    
                    tmp_scrub_final = logical(zeros(size(tmp_scrub)));
            end
                
            % Removes the excess data points
            tmp_data(:,tmp_scrub_final) = [];
            tmp_time_label(tmp_scrub_final) = [];
            
            % This matrix contains the removal information for all subjects
            removed_DPs(idx_subject,run,:) = tmp_scrub_final;
            
            % If we have enough data points, we add the subject data
            if size(tmp_data,2) >= n_TP
                
                % Which null type do we want?
                switch null_type

                    % This retains the data as they are
                    case 'None'
                        Data = [Data,tmp_data];
                        subj_label = [subj_label,idx_subject*ones(1,size(tmp_data,2))];
                        run_label = [run_label,run*ones(1,size(tmp_data,2))];
                        time_label = [time_label,tmp_time_label];

                    % This randomly shuffles the time points
                    case 'Shuffle'
                        Data = [Data,tmp_data(:,randperm(size(tmp_data,2)))];
                        subj_label = [subj_label,idx_subject*ones(1,size(tmp_data,2))];
                        run_label = [run_label,run*ones(1,size(tmp_data,2))];
                        time_label = [time_label,tmp_time_label];
                        
                    case 'PR'
                        tmp = squeeze(TC{subjects_OI(idx_subject)}(1:region_end,:,run))';
                        
                        if ~isempty(PHI)
                            surr = CBIG_RL2017_get_PR_surrogate_RL(tmp,1,PHI);
                        else
                            surr = CBIG_RL2017_get_PR_surrogate_RL(tmp,1,[]);
                        end
                        
                        Data = [Data,surr'];
                        subj_label = [subj_label,idx_subject*ones(1,size(surr,1))];
                        run_label = [run_label,run*ones(1,size(surr,1))];
                        time_label = [time_label,tmp_time_label];

                    % This computes amplitude-adjusted phase-randomised null
                    case 'AAFT'
                        tmp = squeeze(TC{subjects_OI(idx_subject)}(1:region_end,:,run))';
                        
                        if ~isempty(PHI)
                            surr = generate_AAFT_RL(tmp,1,PHI,GAUSS);
                        else
                            surr = generate_AAFT_RL(tmp,1,[],[]);
                        end
                        
                        Data = [Data,surr'];
                        subj_label = [subj_label,idx_subject*ones(1,size(surr,1))];
                        run_label = [run_label,run*ones(1,size(surr,1))];
                        time_label = [time_label,tmp_time_label];

                    otherwise
                        errordlg('Which null scheme is this? Not recognized...');
                end
                
                n_added_runs = n_added_runs + 1;
            else
                
                disp(['PROBLEM WITH SUBJECT ',num2str(idx_subject)]);
                n_removed_runs = n_removed_runs + 1;
                
            end
        end
        
        idx_subject = idx_subject + 1;
    end
    
    % We want time as the first dimension for the output
    Data = Data';
end