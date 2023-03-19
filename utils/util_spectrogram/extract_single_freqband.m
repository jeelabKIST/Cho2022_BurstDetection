function [btime, burst] = extract_single_freqband(btime, burst, idx)
    %% Function: 'extract_single_freqband'
    % DESCRIPTION
    % This function extracts the time periods and spectral powers of the bursts
    % occurred in a single-frequency band.
    
    % INPUT
    %    Variable    Data Type            Description
    % 1. btime       [N x N cell ]      : time periods of the bursts detected from a spectrogram
    % 2. burst       [N x N cell ]      : spectral power of the bursts detected from a spectrogram
    % 3. idx         [integer Z >= 0]   : index of a specific frequency band
    
    % OUTPUT
    %    Variable    Data Type            Description
    % 1. btime       [1 x N cell]       : time periods of the bursts at a single-frequency band
    % 2. burst       [1 x N cell]       : spectral power of the bursts at a specific frequency band
    
    % Written by SungJun Cho, August 25, 2021
    % Last modified on October 29, 2021
    %% Extract burst information from single frequency-band
    btime = btime(idx,:);
    burst = burst(idx,:);
    btime = btime(~cellfun('isempty',btime));
    burst = burst(~cellfun('isempty',burst));
end