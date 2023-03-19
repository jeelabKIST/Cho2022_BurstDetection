function [freqs] = get_multifreq(f, Fs)
    %% Function: 'get_multifreq'
    % DESCRIPTION
    % Computes all the frequencies used in the burst simluation

    % USAGE
    % Full Input : get_multifreq(f, Fs)
    % Example    : get_multifreq(9, 512)

    % INPUT
    %    Variable       Data Type            Description
    % 1. f              [number N]         : lowest frequency for building the burst components
    % 2. Fs             [number N]         : sampling frequency

    % OUTPUT
    % INPUT
    %    Variable       Data Type            Description
    % 1. freqs          [N x 1 vector]     : list of frequencies used to simluate the bursts

    % Written by SungJun Cho, January 27, 2023
    % Last modified on March 14, 2023
    %% Compute Frequencies
    steps = 4:5; % constants used for defining higher frequencies
    nOsc = length(steps)+1; % number of burst types
    freqs = zeros(1,nOsc);
    for i = 1:length(steps)
        n = steps(i);
        f_prime = f+(Fs/(2^n));
        freqs(i) = f_prime;
    end
    freqs(end) = f;
    freqs = sort(freqs);
end