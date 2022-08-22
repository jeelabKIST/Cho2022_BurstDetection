function [matV_best, matC_best, matV_poor, matC_poor] = get_best_and_poor_heatmaps(listNoise, listCycle, dataset, convert_opt, round_opt)
    %% Function 'get_best_and_poor_heatmaps'
    % DESCRIPTION
    % Constructs heatmaps with the values of best and worst perfomring
    % algorithms
    
    % INPUT
    %    Variable      Data Type           Description
    % 1. listNoise     [1 x N vector]    : list of SNR levels (dB)
    % 2. listCycle     [1 x N vector]    : list of burst durations (cycles)
    % 3. dataset       [struct / cell]   : heatmaps of individual algorithms constructed with varying statistical measures
    % 4. convert_opt   [boolean]         : convert values in a percentage scale
    %                                      1) true
    %                                      2) false
    % 5. round_opt     [string]          : round values to integer or first digit decimal
    %                                      1) 'integer'
    %                                      2) 'decimal'
    
    % OUTPUT
    %    Variable      Data Type           Description
    % 1. matV_best     [N x N matrix]    : heatmaps with values of best algorithms
    % 2. matC_best     [N x N matrix]    : color matrix of best algorithms
    % 3. matV_poor     [N x N matrix]    : heatmaps with values of worst algorithms
    % 4. matC_poor     [N x N matrix]    : color matrix of worst algorithms
    
    % Written by SungJun Cho, October 11, 2021
    % Last Modified on October 29, 2021
    %% Get Data and Preallocate Outputs
    nNoise = length(listNoise);
    nCycle = length(listCycle);
    if isa(dataset,'struct')
        nMethod = length(fieldnames(dataset));
    else
        nMethod = length(dataset);
    end
    hp = helper;
    [V_bp,V_ev,V_stp,V_mtp,V_cwt] = hp.unpack_data(dataset);
    if convert_opt
        [V_bp,V_ev,V_stp,V_mtp,V_cwt] = hp.convert_to_percent(V_bp,V_ev,V_stp,V_mtp,V_cwt);
    end
    if strcmp(round_opt,'integer')
        [V_bp,V_ev,V_stp,V_mtp,V_cwt] = hp.dec2int(V_bp,V_ev,V_stp,V_mtp,V_cwt);
    elseif strcmp(round_opt,'decimal')
        [V_bp,V_ev,V_stp,V_mtp,V_cwt] = hp.round_dec(V_bp,V_ev,V_stp,V_mtp,V_cwt);
    end
    [matV_best, matC_best, matV_poor, matC_poor] = deal(zeros(nNoise,nCycle));
    %% Compute Values of Best and Worst Performing Algorithms
    for i = 1:nNoise
        for j = 1:nCycle
            val_set = [V_cwt(i,j),V_mtp(i,j),V_stp(i,j),V_ev(i,j),V_bp(i,j)];
            % Get Best Performances
            [max_val, max_idx] = max(val_set);
            matV_best(i,j) = max_val;
            if length(unique(val_set)) < 5 && sum(val_set == max_val) > 1
                if length(unique(val_set)) == 1
                    matC_best(i,j) = nMethod + 2; % same values for all matrices
                else
                    matC_best(i,j) = nMethod + 1; % same values for multiple matrices
                end
            else
                matC_best(i,j) = max_idx;
            end
            % Get Poor Performances
            [min_val, min_idx] = min(val_set);
            matV_poor(i,j) = min_val;
            if length(unique(val_set)) < 5 && sum(val_set == min_val) >1
                if length(unique(val_set)) == 1
                    matC_poor(i,j) = nMethod + 2; % same values for all matrices
                else
                    matC_poor(i,j) = nMethod + 1; % same values for multiple matrices
                end
            else
                matC_poor(i,j) = min_idx;
            end
        end
    end
end