function [bscoreMat] = map_burst_metric(scoreMat, btimeMat, dtMat, tvec, signal, Fs, lo_f, hi_f)
    %% Function 'map_burst_metric'
    % DESCRIPTION
    % Maps burst detections to their corresponding performance metric
    % values for each algorithm
    
    % INPUT
    %    Variable       Data Type           Description
    % 1. scoreMat       [1 x N cell]      : set of heatmaps of a chosen performance metric for each algorithm
    % 2. btimeMat       [1 x N cell]      : set of time arrays of detected bursts for each algorithm
    % 3. dtMat          [1 x N cell]      : set of a sampling rate or sliding window time step for each algorithm
    % 4. tvec           [1 x N double]    : time array for the input signal
    % 5. signal         [1 x N double]    : input signal
    % 6. Fs             [number N]        : sampling frequency
    % 7. lo_f           [number N]        : low cutoff frequency of the input signal
    % 8. hi_f           [number N]        : high cutoff frequency of the input signal
    
    % OUTPUT
    %    Variable       Data Type           Description
    % 1. bscoreMat     [1 x N cell]       : performance metric values mapped from burst detections of each algorithm
    
    % Written by SungJun Cho, April 30, 2023
    % Last Modified on April 30, 2023
    %% Set Parameters
    nMethod = length(btimeMat); % number of algorithms
    %% Compute Burst Cycle Lengths
    foi = (lo_f + hi_f) / 2;
    wrapper = @(x,y) sec2cyc(x,foi,y);
    blenMat = cellfun(wrapper, btimeMat, dtMat, 'UniformOutput', false); % length (in cycles) of all bursts
    % Since bursts longer than 12 cycles or shorter than 3 cycles can be 
    % detected, we bound maximum length of bursts to 12 cycles and its 
    % minimum to 3 cycles in order to compare them with decision matrices.
    for n = 1:nMethod
        for i = 1:length(blenMat{n})
            if blenMat{n}(i) > 12
                blenMat{n}(i) = 12;
            end
            if blenMat{n}(i) < 3
                blenMat{n}(i) = 3;
            end
        end
    end
    %% Compute Burst Signal-to-Noise Ratio
    snrMat = get_snr(btimeMat,tvec,signal,Fs,lo_f,hi_f);
    % Since bursts with SNR higher than 10 dB or lower than -10 dB can be 
    % detected, we bound maximum SNR of bursts to 10 dB and minimum to -10 
    % dB in order to compare them with decision matrices.
    for n = 1:nMethod
        for i = 1:length(snrMat{n})
            if snrMat{n}(i) > 10
                snrMat{n}(i) = 10;
            end
            if snrMat{n}(i) < -10
                snrMat{n}(i) = -10;
            end
        end
    end
    %% Get Metric Scores for Each Burst
    cycList = 3:12;
    snrList = -10:2:10;
    bscoreMat = cell(size(btimeMat)); % metric scores for each burst
    for n = 1:nMethod
        score = scoreMat{n};
        blen = blenMat{n};
        bsnr = snrMat{n};
        if length(blen) < 2 % if no bursts or only one burst is detected
            bscoreMat{n} = NaN;
            continue;
        end
        if length(blen) == length(bsnr)
            bscore = zeros(size(blen));
            for i = 1:length(bscore)
                cycIdx = find(blen(i) == cycList);
                if mod(bsnr(i),2) == 0
                    bscore(i) = score(bsnr(i) == snrList,cycIdx);
                else
                    bscore(i) = mean([score(bsnr(i)-1 == snrList,cycIdx), ... 
                        score(bsnr(i)+1 == snrList,cycIdx)]);
                end
            end
        end
        bscoreMat{n} = bscore;
    end
end
%% Appendix: In-Script Functions
% Function #1: Converts burst lengths from samples to cycles
function [num_cycles] = sec2cyc(btime,foi,dt)
    num_cycles = cellfun(@(b) round(length(b)*(foi*dt)), btime);
end

% Function #2: Compute SNR (dB) of each burst
function [snrMat] = get_snr(btimeMat,tvec,signal,Fs,lo_f,hi_f)
    % [1] Get bandpass-filtered signals with and without bursts
    filt_sig = eegfilt(signal,Fs,lo_f,hi_f);
    % [2] Measure SNR for each burst
    snrMat = cell(size(btimeMat));
    for i = 1:length(btimeMat)
        btime = btimeMat{i};
        bsig = cell(size(btime));
        qnoise = ones(1,length(tvec));
        for n = 1:length(btime)
            % (A) For time domain methods
            if i < 3
                bIdx = ismember(tvec,btime{n});
                bsig{n} = filt_sig(bIdx);
                qnoise(bIdx) = 0;
            % (B) For time-frequency domain methods
            else
                binterval = [btime{n}(1) btime{n}(end)];
                bIdx = find(tvec >= binterval(1) & tvec <= binterval(end));
                bsig{n} = filt_sig(bIdx);
                qnoise(bIdx) = 0;
            end
        end
        amp_sig = cell2mat(cellfun(@(x) rms(x), bsig, 'UniformOutput', false));
        amp_noise = rms(filt_sig(logical(qnoise)));
        snrMat{i} = round(10*log10((amp_sig/amp_noise).^2));
    end
end