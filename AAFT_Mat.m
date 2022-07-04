function [Data2] = AAFT_Mat(Data)

    for r = 1:size(Data,1)
        Data2(r,:) = generate_AAFT_RL(Data(r,:)',1);
    end

end