X_gsr = [];
X_measure = [];
X_runs = [];

X_mot = [];
X_sampling = [];



% All similar except for the run
for a = 1:2

        for c = 1:2
            for d = 1:4
                for e = 1:4
                    
                    X_runs = [X_runs,SM(gsr_vec==a & measure_vec==d & mot_vec == e & runs_vec == 1 & sampling_vec == c, gsr_vec==a & measure_vec==d & mot_vec == e & runs_vec == 2 & sampling_vec == c)];
                    
                    
                end
            end
        end
    
end


% All similar except for the GSR
for a = 1:2

        for c = 1:2
            for d = 1:4
                for e = 1:4
                    
                    X_gsr = [X_gsr,SM(gsr_vec==1 & measure_vec==d & mot_vec == e & runs_vec == a & sampling_vec == c, gsr_vec==2 & measure_vec==d & mot_vec == e & runs_vec == a & sampling_vec == c)];
                    
                    
                end
            end
        end
    
end

% All similar except for the measure
for a = 1:2

        for c = 1:2
            for d = 1:4
                for e = 1:4
                    
                    X_measure = [X_measure,SM(gsr_vec==d & measure_vec==1 & mot_vec == e & runs_vec == a & sampling_vec == c, gsr_vec==d & measure_vec==2 & mot_vec == e & runs_vec == a & sampling_vec == c)];
                    
                    
                end
            end
        end
    
end



% All similar except for the motion correction
for a = 1:2

        for c = 1:2
            for d = 1:2
                for e = 1:4
                    
                    X_mot = [X_mot,[SM(gsr_vec==d & measure_vec==c & mot_vec == 1 & runs_vec == a & sampling_vec == e, gsr_vec==d & measure_vec==c & mot_vec == 2 & runs_vec == a & sampling_vec == e);...
                        SM(gsr_vec==d & measure_vec==c & mot_vec == 1 & runs_vec == a & sampling_vec == e, gsr_vec==d & measure_vec==c & mot_vec == 3 & runs_vec == a & sampling_vec == e);...
                        SM(gsr_vec==d & measure_vec==c & mot_vec == 1 & runs_vec == a & sampling_vec == e, gsr_vec==d & measure_vec==c & mot_vec == 4 & runs_vec == a & sampling_vec == e);...
                        SM(gsr_vec==d & measure_vec==c & mot_vec == 2 & runs_vec == a & sampling_vec == e, gsr_vec==d & measure_vec==c & mot_vec == 3 & runs_vec == a & sampling_vec == e);...
                        SM(gsr_vec==d & measure_vec==c & mot_vec == 2 & runs_vec == a & sampling_vec == e, gsr_vec==d & measure_vec==c & mot_vec == 4 & runs_vec == a & sampling_vec == e);...
                        SM(gsr_vec==d & measure_vec==c & mot_vec == 3 & runs_vec == a & sampling_vec == e, gsr_vec==d & measure_vec==c & mot_vec == 4 & runs_vec == a & sampling_vec == e)]];
                end
            end
        end
    
end


% All similar except for the sampling scheme
for a = 1:2
        for c = 1:2
            for d = 1:2
                for e = 1:4
                    
                    X_sampling = [X_sampling,[SM(gsr_vec==d & measure_vec==c & mot_vec == e & runs_vec == a & sampling_vec == 1,gsr_vec==d & measure_vec==c & mot_vec == e & runs_vec == a & sampling_vec == 2);...
                        SM(gsr_vec==d & measure_vec==c & mot_vec == e & runs_vec == a & sampling_vec == 1,gsr_vec==d & measure_vec==c & mot_vec == e & runs_vec == a & sampling_vec == 3);...
                        SM(gsr_vec==d & measure_vec==c & mot_vec == e & runs_vec == a & sampling_vec == 1,gsr_vec==d & measure_vec==c & mot_vec == e & runs_vec == a & sampling_vec == 4);...
                        SM(gsr_vec==d & measure_vec==c & mot_vec == e & runs_vec == a & sampling_vec == 2,gsr_vec==d & measure_vec==c & mot_vec == e & runs_vec == a & sampling_vec == 3);...
                        SM(gsr_vec==d & measure_vec==c & mot_vec == e & runs_vec == a & sampling_vec == 2,gsr_vec==d & measure_vec==c & mot_vec == e & runs_vec == a & sampling_vec == 4);...
                        SM(gsr_vec==d & measure_vec==c & mot_vec == e & runs_vec == a & sampling_vec == 3,gsr_vec==d & measure_vec==c & mot_vec == e & runs_vec == a & sampling_vec == 4)]];
                end
            end
        end
end
          
                        
                        