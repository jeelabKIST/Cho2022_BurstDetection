function [btime, burst] = detect_burst_timeseries(Fs, times, signal, low, high, thrA_factor)
    %% Function: 'detect_burst_timeseries'
    % DESCRIPTION
    % Detects bursts in the time domain by imposing an amplitude and a duration
    % threshold on a bandpass-filtered signal

    % USAGE
    % Full Input : detect_burst_timeseries(Fs, times, filtered_signal, f_low, f_high, thrA_factor)
    % Example    : detect_burst_timeseries(512, tvec, filt_sig, 20, 30, 2.0)

    % INPUT
    %    Variable       Data Type                 Description
    % 1. Fs             [number N]              : sampling rate
    % 2. times          [1 x N vector]          : time array for the input signal
    % 3. signal         [1 x N vector]          : input signal
    %                                             Note) Must be bandpass-filtered before using this function
    % 4. low            [number N]              : low cutoff frequency of the input signal
    % 5. high           [number N]              : high cutoff frequency of the input signal
    % 6. thrA_factor    [number N]              : amplitude threshold factor (i.e., multiplier of standard deviation)

    % OUTPUT
    % INPUT
    %    Variable       Data Type                 Description
    % 1. btime          [1 x N cell vector]     : set of time arrays of the detected bursts
    % 2. burst          [1 x N cell vector]     : set of amplitude arrays of the detected bursts

    % Written by SungJun Cho, April 3, 2021
    % Modified on October 29, 2021
    %% Set Parameters
    % [1] Amplitude Thresholds
    % Define amplitude threshold and subthreshold by first and second moments
    thr = mean(signal) + thrA_factor*std(signal);
    subthr = mean(signal) + (thrA_factor-1)*std(signal);
    abvthr = find(signal > thr);          % find indices above the set threshold
    belsubthr = find(signal <= subthr);   % find indices below or equal to the set subthreshold

    % [2] Duration Thresholds
    distinction_thr = Fs/((low+high)/2);  % Threshold for separation of distinct bursts
                                          % Or alternatively, use: 2 cycle difference with (Fs/low)*2 OR Fs/((low+high)/2)*2
    duration_thr = 3*(Fs/((low+high)/2)); % Threhsold for minimum duration of each burst
    error_thr = Fs/high;                  % Threshold for multiple detections of the same peak
    %% Identify Separate Bursts
    burstend = zeros(1,length(signal));
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
    for i = 1:length(burstend)
        if i == 1
            burst{1,i} = signal(abvthr(i):burstend(i));
            btime{1,i} = times(abvthr(i):burstend(i));
        elseif i > 1
            start = find(abvthr == burstend(i-1),1,'last') + 1;
            burst{1,i} = signal(abvthr(start):burstend(i));
            btime{1,i} = times(abvthr(start):burstend(i));
        end
    end
    %% Eliminate Incorrect Detections and Bursts Shorter Than 3 Cycles
    elimIdx = zeros(1,length(burst));
    for i = 1:length(burst)
        if length(burst{1,i}) <= duration_thr
            elimIdx(1,i) = i;
        elseif length(burst{1,i}) > duration_thr
            % Minimum requirements the burst waveforms must suffice
            [pks,locs] = findpeaks(burst{1,i}); % function 'findpeaks' needs at least 2 samples
            pks = pks(pks > thr);  % for cases of incorrect detection
            for j = 1:length(locs) % for cases of multiple detection of same peak
                if j == length(locs)
                    break
                end
                if locs(j+1) - locs(j) < error_thr
                    locs(j) = 0;
                end
            end
            locs = nonzeros(locs)';
            if length(pks) < 3 || length(locs) < 3 % If less than 3 peaks are detected
                elimIdx(1,i) = i;                  % then eliminate the detected burst
            end
        end
    end
    elimIdx = nonzeros(elimIdx)';
    burst(elimIdx) = [];
    btime(elimIdx) = [];

    if length(burst) ~= length(btime)
        error('ERROR: Sizes of time and amplitude arrays of detected bursts do not match.');
    end
    %% Extract Bursts with Subthreshold-Defined Starting and End Points
    for i = 1:length(btime)
        blogic_temp = ismember(times,btime{i});
        blogic_temp = find(blogic_temp ~= 0);
        blogic_start = find(belsubthr < blogic_temp(1), 1, 'last');
        blogic_end = find(belsubthr > blogic_temp(end), 1, 'first');
        blogic_start = belsubthr(blogic_start);
        blogic_end = belsubthr(blogic_end);
        % Consider the cases in which bursts start or end at the beginning and end of a time vector
        if btime{i}(1) == 0
            blogic_start = 1;
        end
        if btime{i}(end) == times(end)
            blogic_end = length(times);
        end
        % Consider the cases in which there are no values below subthreshold at the beginning and end of a time vector
        if isempty(blogic_start)
            blogic_start = find(times == btime{i}(1));
        end
        if isempty(blogic_end)
            blogic_end = find(times == btime{i}(end));
        end
        btime{i} = times(blogic_start:blogic_end);
        burst{i} = signal(blogic_start:blogic_end); 
    end
end