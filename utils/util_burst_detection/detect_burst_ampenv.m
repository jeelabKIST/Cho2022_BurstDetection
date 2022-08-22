function [btime, burst, benvelope] = detect_burst_ampenv(Fs, time, signal, low, high, prcthrA)
    %% Function: 'detect_burst_ampenv'
    % DESCRIPTION
    % Detects bursts in the time domain by imposing an amplitude and a duration
    % threshold on an amplitude envelope of the input bandpass-filtered signal.

    % USAGE
    % Full Input : detect_burst_ampenv(Fs, time, filtered_signal, f_low, f_high, prcthrA)
    % Example    : detect_burst_ampenv(512, tvec, filt_sig, 20, 30, 95)

    % INPUT
    %    Variable       Data Type                 Description
    % 1. Fs             [number N]              : sampling rate
    % 2. time           [1 x N vector]          : time array for the input signal
    % 3. signal         [1 x N vector]          : input signal
    %                                             Note) Must be bandpass-filtered before using this function
    % 4. low            [number N]              : low cutoff frequency of the input signal
    % 5. high           [number N]              : high cutoff frequency of the input signal
    % 6. prcthrA        [number N]              : amplitude threshold in percentile

    % OUTPUT
    % INPUT
    %    Variable       Data Type                 Description
    % 1. btime          [1 x N cell vector]     : set of time arrays of the detected bursts
    % 2. burst          [1 x N cell vector]     : set of amplitude arrays of the detected bursts
    % 3. benvelope      [1 x N cell vector]     : set of amplitude envelope arrays of the detected bursts

    % Written by SungJun Cho, April 4, 2021
    % Modified on October 29, 2021
    %% Get Amplitude Envelope Using Hilbert Transform
    signal_env = abs(hilbert(signal));
    if length(signal) ~= length(signal_env)
        error('ERROR: The lengths of the signal and its envelope do not match.');
    end
    %% Set Parameters
    % [1] Amplitude Thresholds
    % Define amplitude threshold and subthreshold by percentile
    thr = prctile(signal_env,prcthrA);    % set amplitude threshold
    abvthr = find(signal_env > thr);      % find indices above the set threshold

    % [2] Duration Thresholds
    distinction_thr = Fs/((low+high)/2);  % Threshold for separation of distinct bursts
                                          % Or alternatively, use: 2 cycle
                                          % difference with (Fs/low)*2 OR Fs/((low+high)/2)*2
    duration_thr = 3*(Fs/((low+high)/2)); % Threhsold for minimum duration of each burst
    %% Identify Separate Bursts
    burstend = zeros(1,length(signal_env));
    for i = 1:length(abvthr)
        if i == length(abvthr)
            break
        end
        % If indices above the threshold exceeds at least 1 cycle difference,
        % then regard them as separate bursts.
        if abvthr(i+1) - abvthr(i) > distinction_thr
            burstend(1,i) = abvthr(i);
        end
        % Consider the cases where there is only one burst
        if sum(burstend) == 0 && ~isempty(abvthr)
            burstend(1,length(abvthr)+1) = abvthr(end);
        end
        % Include the last detected burst if not already included
        if ~ismember(burstend,abvthr(end))
            burstend(1,length(abvthr)+1) = abvthr(end);
        end
    end
    burstend = nonzeros(burstend)';
    %% Screen for Preliminary Bursts
    burst = cell(1,length(burstend));
    btime = cell(1,length(burstend));
    benvelope = cell(1,length(burstend));
    for i = 1:length(burstend)
        if i == 1
            burst{1,i} = signal(abvthr(i):burstend(i));
            btime{1,i} = time(abvthr(i):burstend(i));
            benvelope{1,i} = signal_env(abvthr(i):burstend(i));
        elseif i > 1
            start = find(abvthr == burstend(i-1),1,'last') + 1;
            burst{1,i} = signal(abvthr(start):burstend(i));
            btime{1,i} = time(abvthr(start):burstend(i));
            benvelope{1,i} = signal_env(abvthr(start):burstend(i));
        end
    end
    %% Eliminate Incorrect Detections and Bursts Shorter Than 3 Cycles
    elimIdx = zeros(1,length(benvelope));
    for i = 1:length(benvelope)
        if length(benvelope{1,i}) <= duration_thr
            elimIdx(1,i) = i;
        elseif length(benvelope{1,i}) > duration_thr
            continue;
        end
    end
    elimIdx = nonzeros(elimIdx)';
    burst(elimIdx) = [];
    btime(elimIdx) = [];
    benvelope(elimIdx) = [];

    if length(burst) ~= length(btime) && length(burst) ~= length(benvelope)
        error('ERROR: Sizes of time and amplitude arrays of detected bursts do not match.');
    end
end