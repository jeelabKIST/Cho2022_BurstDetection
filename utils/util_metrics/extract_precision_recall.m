function [prf_per_snr] = extract_precision_recall(statistics, listCycle)
    %% Function 'extract_precision_recall'
    % DESCRIPTION
    % Computes precision, recall, and F1-scores for bursts detected from a
    % set of synthetic signals with a particular SNR level.
    
    % INPUT
    %    Variable         Data Type           Description
    % 1. statistics       [1 x N cell]      : burst statistics
    % 2. listCycle        [1 x N double]    : list of burst durations (cycles)
    
    % OUTPUT
    %    Variable         Data Type           Description
    % 1. prf_per_snr      [1 x N double]    : precision, recall, and F1-score values 
    %                                         for each burst duration and SNR noise level
    
    % Written by SungJun Cho, October 11, 2021
    % Last Modified on February 11, 2023
    %% Compute Precision, Recall, and F1-Score
    statistics = [statistics{:}];
    statistics = sortrows(statistics(3:4,:)',2);
    prf_per_snr = zeros(3,length(listCycle));
    uniqueCycle = unique(cell2mat(statistics(:,2)));
    if sum(uniqueCycle < 3) ~= 0 || sum(uniqueCycle > 12) ~= 0
        out_of_range = uniqueCycle(~ismember(uniqueCycle,listCycle));
        for i = 1:length(out_of_range)
            c = out_of_range(i);
            idx = cellfun(@(x) x==c, statistics(:,2));
            if c < listCycle(1)
                error('ComputationError: Some detected bursts did not satisfy the duration threshold.');
            end
            if c > listCycle(end)
                statistics(idx,2) = {listCycle(end)};
            end
        end
    end
    for i = 1:length(listCycle)
        c = listCycle(i);
        idx = cellfun(@(x) x==c, statistics(:,2));
        stat_per_cyc = statistics(idx,1);
        TP = sum(ismember(stat_per_cyc,'tp'));
        FP = sum(ismember(stat_per_cyc,'fp'));
        FN = sum(ismember(stat_per_cyc,'fn'));
        P = TP/(TP+FP);    % precision
        R = TP/(TP+FN);    % recall
        F = (2*P*R)/(P+R); % F1-score
        prf_per_snr(1,i) = P;
        prf_per_snr(2,i) = R;
        prf_per_snr(3,i) = F;
    end
end