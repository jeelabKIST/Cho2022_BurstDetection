function [tc_per_snr] = extract_temporal_concurrence(statistics, listCycle)
    %% Function 'extract_temporal_concurrence'
    % DESCRIPTION
    % Computes temporal concurrence between the simulated and detected bursts 
    % from a set of synthetic signals with a particular SNR level. Detected
    % bursts should be true positives.
    
    % INPUT
    %    Variable         Data Type           Description
    % 1. statistics       [1 x N cell]      : burst statistics
    % 2. listCycle        [1 x N double]    : list of burst durations (cycles)
    
    % OUTPUT
    %    Variable         Data Type           Description
    % 1. tc_per_snr       [1 x N double]    : temporal concurrence values for each
    %                                         burst duration and SNR noise level
    
    % Written by SungJun Cho, October 11, 2021
    % Last Modified on February 11, 2023
    %% Compute Temporal Concurrence
    statistics = [statistics{:}];
    tpIdx = strcmp(statistics(3,:),'tp');
    statistics = statistics(:,tpIdx);
    statistics = sortrows(statistics(4:end,:)',1);
    tc_per_snr = zeros(1,length(listCycle));
    for i = 1:length(listCycle)
        c = listCycle(i);
        idx = cellfun(@(x) x==c, statistics(:,1));
        if ~any(idx)
            tc = NaN;
        else
            tc = mean(cell2mat(statistics(idx,end)));
        end
        tc_per_snr(1,i) = tc;
    end
end