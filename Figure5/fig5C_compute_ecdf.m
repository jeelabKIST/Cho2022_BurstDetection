%% Configure Library Paths
util_path = genpath('/Users/jeelab/Desktop/Cho2022_BurstDetection/utils');
data_path = genpath('/Users/jeelab/Desktop/Cho2022_BurstDetection/data/experimental_data/ESCAPE');
addpath(util_path);
addpath(data_path);
%% Load Data
% [1] Get Theoretical Results
HEATMAP = load('/Users/jeelab/Desktop/Cho2022_BurstDetection/data/simulation_data/HM_gamma.mat').HEATMAP;
DET_CONF = load('/Users/jeelab/Desktop/Cho2022_BurstDetection/data/simulation_data/DC_gamma.mat');
scoreMat = {
    struct2cell(HEATMAP.precision)';
    struct2cell(HEATMAP.recall)';
    struct2cell(HEATMAP.f1_score)';
    struct2cell(HEATMAP.concurrence)';
    DET_CONF.detection_confidence
};
metric_names = {'Precision','Recall','F1-Score','Temporal Concurrence','Detection Confidence'};
method_names = {'BP','ENV','S-STFT','MTP','CWT'};
nMethod = length(method_names);
% [2] Get LFP Signal Data
load('lfps/escape_eeg_all.mat');
anatomy_idx = 1; % 1-BLA; 2-PFC; 3-HPC
Fs = srate;      % sampling rate
Nyq = Fs/2;      % Nyquist frequency
tvec = times;    % time vector
nTrial = size(data,3);
% (***) NOTE: Set nTrial=1 if you want to visualize eCDFs and their AUC values
% for Mouse 1 Day 1 Session 1 (Figure 5C).
signals = squeeze(data(anatomy_idx,:,:));
% [3] Grab Helper Functions
hp = helper;
%% Set Hyperparameters For Each Algorithm
% [1] Time Domain
f = 40;     % frequency of interest
lo_f = 35;  % lower frequency bound
hi_f = 45;  % upper frequency bound
mid_f =  (lo_f + hi_f) / 2;
% [2] Time-Frequency Domain
freq_range = unique([15:5:30 30:10:70 70:15:100]); % frequencies of interest
dt = 0.01;                                         % time step for sliding windows
nOscCycle_stp = 6; nOscCycle_mtp = 6;              % frequency-dependent window size
nw = 2; ntp = 3;                                   % DPSS slepian sequences
wname = 'amor'; scale_opt = 'default';             % type of and scaling option for the wavelet
interp_opt = false;
check_opt = true;
plot_psd = false;
plot_tp = false;
verbose = false;
% For `freq_range`, you can alternatively use:
%   freq_range = 2.^((1/5).*(1:50))*0.1; % f_min = 0.1; num_octave = 5;
%   freq_range = freq_range(freq_range > 5);
% [3] Amplitude Thresholds
thrA_factor = 2.5;       % for bandpass filtering method
prcthrA = 95;            % for envelope-based method
specthrA_factor = 1.8;   % for spectral analysis methods
% [4] Window Time Steps
dt_stp = dt;
dt_mtp = dt;
dt_cwt = 1/Fs;
%% Empirical Cumulative Distribution Functions (eCDF) and Area Under Curve (AUC)
tic; fprintf('Computing eCDF and AUC ... \n');
AUC = struct();
[AUC.precision, AUC.recall,AUC.f1_score,AUC.concurrence,AUC.detection_confidence] = deal(zeros(nTrial,nMethod));
for n = 1:nTrial
    % [1] Get Burst Time Intervals Using Each Algorithm (Burst Detections)
    signal = double(signals(:,n)');
    % (A) Apply IIR Notch Filter
    w0 = 60/Nyq;
    q_factor = 35;
    bw = w0/q_factor;
    [b,a] = iirnotch(w0,bw);
    signal = filtfilt(b,a,signal);
    % (B) Bandpass Filtering
    filt_sig = eegfilt(signal,Fs,lo_f,hi_f);
    [btime_bp,burst_bp] = detect_burst_timeseries(Fs,tvec,filt_sig,lo_f,hi_f,thrA_factor);
    % (C) Amplitude Envelope
    [btime_ev,burst_ev,~] = detect_burst_ampenv(Fs,tvec,filt_sig,lo_f,hi_f,prcthrA);
    % (D) Single-Tapered Short Time Fourier Transform
    [Spec_f_stp,Spec_t_stp,Spec_stp,~] = spectrogram_stp(tvec,signal,Fs,freq_range,dt_stp,nOscCycle_stp,verbose,interp_opt);
    [btime_stp,burst_stp,~] = detect_burst_spectrogram(Spec_f_stp,Spec_t_stp,Spec_stp,dt_stp,specthrA_factor,verbose);
    minIdx_stp = find(abs(Spec_f_stp-mid_f) == min(abs(Spec_f_stp-mid_f)));
    [btime_stp,burst_stp] = extract_single_freqband(btime_stp,burst_stp,minIdx_stp);
    % (E) Multitaper Spectrogram
    [Spec_f_mtp,Spec_t_mtp,Spec_mtp,~] = spectrogram_mtp(tvec,signal,Fs,nw,ntp,freq_range,dt_mtp,nOscCycle_mtp,interp_opt,plot_psd,check_opt,plot_tp,verbose);
    [btime_mtp,burst_mtp,~] = detect_burst_spectrogram(Spec_f_mtp,Spec_t_mtp,Spec_mtp,dt_mtp,specthrA_factor,verbose);
    minIdx_mtp = find(abs(Spec_f_mtp-mid_f) == min(abs(Spec_f_mtp-mid_f)));
    [btime_mtp,burst_mtp] = extract_single_freqband(btime_mtp,burst_mtp,minIdx_mtp);
    % (F) Continuous Wavelet Spectrogram
    [wav_freq,wav_time,Spec_cwt,~] = spectrogram_cwt(tvec,signal,Fs,wname,scale_opt,verbose);
    [btime_cwt,burst_cwt,~] = detect_burst_spectrogram(wav_freq,wav_time,Spec_cwt',dt_cwt,specthrA_factor,verbose);
    minIdx_cwt = find(abs(wav_freq-mid_f) == min(abs(wav_freq-mid_f)));
    [btime_cwt,burst_cwt] = extract_single_freqband(btime_cwt,burst_cwt,minIdx_cwt);
    % [2] Parse Bursts (Experiment Stage 2)
    btime_bp = parse_data(btime_bp);
    btime_ev = parse_data(btime_ev);
    btime_stp = parse_data(btime_stp);
    btime_mtp = parse_data(btime_mtp);
    btime_cwt = parse_data(btime_cwt);
    % [3] Compute eCDF and AUC
    btimeMat = hp.pack_data(btime_bp, btime_ev, btime_stp, btime_mtp, btime_cwt);
    dtMat = hp.pack_data(1/Fs, 1/Fs, dt_stp, dt_mtp, dt_cwt);
    % [4] Plot Intersecting Bursts and eCDFs
    verbose_cdf = true;
    indicator = {'precision','recall','f1_score','concurrence','detection_confidence'};
    for m = 1:length(scoreMat)
        if n > 1
            verbose_cdf = false;
        end
        [bscoreMat, auc] = select_best_algorithm(n,metric_names{m},scoreMat{m},btimeMat,dtMat, ... 
                           tvec,signal,Fs,lo_f,hi_f,verbose_cdf);
        AUC.(indicator{m})(n,:) = auc;
    end
    % [5] Report Progress
    if mod(n,10) == 0
        fprintf(['Processed: ' num2str(n) '/' num2str(nTrial) ' trials \n']);
    end
end
elapsed_time = toc;
fprintf(['Computation complete. (computational time: ' num2str(elapsed_time) 's) \n']);
%% Save AUC Results
file_name = 'AUC_stage2.mat';
if ~isfile(file_name)
    save(file_name,'AUC');
else
    error('Save Error: The specified file already exists.');
end
%% Appendix: Accessory Functions
% Function #1: Plots time intervals of the detected bursts
function plot_burst_raster_plot(btimeMat)
    nMethod = length(btimeMat);
    figure(); hold on;
    for n = 1:nMethod
        btime = btimeMat{n};
        for b = 1:length(btime)
            scatter(btime{b}, ones(size(btime{b}))*n, 'k', 'filled');
        end
    end
    ylim([0 nMethod+1]);
    yticks(1:nMethod); yticklabels({'BP','ENV','STP','MTP','CWT'});
    xlabel('Time (s)');
    ylabel('Algorithm');
    title('Rater Plots');
    set(gca,'FontSize',15,'TickDir','out','LineWidth',2);
    set(gcf,'Color','w');
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

% Function #3: Parse Bursts by Time
function [btime] = parse_data(btime)
    idx = cellfun(@(b) b(1) >= 60 && b(end) <=120, btime);
    btime = btime(idx);
end

% Function #4: Converts burst lengths from samples to cycles
function [num_cycles] = sec2cyc(btime,foi,dt)
    num_cycles = cellfun(@(b) round(length(b)*(foi*dt)), btime);
end

% Function #5: Plots empirical cumulative distribution functions (eCDF) and calculates their AUC
function [bscoreMat, auc] = select_best_algorithm(trialID,metric_name,scoreMat,btimeMat,dtMat,tvec,signal,Fs,lo_f,hi_f,verbose)
    if nargin < 11
        verbose = true;
    end
    nMethod = length(btimeMat); % number of algorithms
    okabe_ito = [[0,114,178];
                 [213,94,0];
                 [0,158,115];
                 [204,121,167];
                 [240,228,66]]./255; % color map
    % [1] Visualize Burst Time Points
    if verbose
        plot_burst_raster_plot(btimeMat);
    end
    % [2] Compute Burst Cycle Lengths
    foi = (lo_f + hi_f) / 2;
    wrapper = @(x,y) sec2cyc(x,foi,y);
    blenMat = cellfun(wrapper, btimeMat, dtMat, 'UniformOutput', false); % length (in cycles) of all bursts
    % Since bursts longer than 12 cycles can be detected, we bound maximum
    % length of bursts to 12 cycles in order to compare them with decision
    % matrices.
    for n = 1:nMethod
        for i = 1:length(blenMat{n})
            if blenMat{n}(i) > 12
                blenMat{n}(i) = 12;
            end
        end
    end
    % [3] Compute Burst Signal-to-Noise Ratio
    snrMat = get_snr(btimeMat,tvec,signal,Fs,lo_f,hi_f);
    % Since bursts with SNR higher than 10 dB can be detected, we bound maximum
    % SNR of bursts to 10 dB in order to compare them with decision
    % matrices.
    for n = 1:nMethod
        for i = 1:length(snrMat{n})
            if snrMat{n}(i) > 10
                snrMat{n}(i) = 10;
            end
        end
    end
    % [4] Get Detection Confidence for Each Burst
    cycList = 3:12;
    snrList = -10:2:10;
    bscoreMat = cell(size(btimeMat)); % detection confidence scores for each burst
    for n = 1:nMethod
        score = scoreMat{n};
        blen = blenMat{n};
        bsnr = snrMat{n};
        if length(blen) == length(bsnr)
            bscore = zeros(size(blen));
            for i = 1:length(bscore)
                cycIdx = find(blen(i) == cycList);
                if mod(bsnr(i),2) == 0
                    snrIdx = find(bsnr(i) == snrList);
                    bscore(i) = score(snrIdx,cycIdx);
                else
                    snrIdx_lo = find(bsnr(i)-1 == snrList);
                    snrIdx_hi = find(bsnr(i)+1 == snrList);
                    bscore(i) = mean([score(snrIdx_lo,cycIdx),score(snrIdx_hi,cycIdx)]);
                end
            end
        end
        bscoreMat{n} = bscore;
    end
    % [5] Compute and Plot Empirical Cumulative Distribution Function
    [cdf_x, cdf_y] = deal(cell(1,nMethod));
    auc = zeros(1,nMethod);
    figure(); hold on;
    xmin = 0;
    xmax = 0;
    h_set = zeros(1,nMethod);
    for n = 1:nMethod
        [h, ~] = cdfplot(sort(bscoreMat{n}));
        h.Color = okabe_ito(n,:);
        h.Color(4) = 0.9;
        h.LineWidth = 4.5;
        h.XData(1) = 0;
        h.XData(end) = max(cat(2,bscoreMat{:}),[],'all');
        if xmax < h.XData(end)
            xmax = h.XData(end);
        end
        cdf_x{n} = h.XData;
        cdf_y{n} = h.YData;
        auc(n) = trapz(cdf_x{n},cdf_y{n});
        h_set(n) = h;
    end
    baseline = plot(linspace(xmin, xmax, 100),linspace(0, 1, 100),'Color','k','LineWidth',3.5,'LineStyle','--');
    baseline.Color(4) = 0.5;
    uistack(baseline,'bottom');
    h.Parent.XLabel.String = metric_name;
    h.Parent.YLabel.String = 'Cumulative Probability';
    title({['Example: Trial #' num2str(trialID)], [metric_name ' eCDF']});
    yticks(linspace(0,1,6));
    leg_labels = {'BP','ENV','S-STFT','MTP','CWT'};
    for n = 1:nMethod
        leg_labels{n} = [leg_labels{n} ' (AUC = ' num2str(auc(n),'%.3f') ')'];
    end
    leg = legend(h_set,leg_labels,'Location','best');
    leg.ItemTokenSize = [15,18];
    grid off;
    set(leg,'FontSize',24,'FontWeight','bold');
    set(gca,'FontSize',35,'FontWeight','bold','LineWidth',5,'YLim',[-0.01,1.01],'TickDir','Out','TickLength',[0.015, 0.025], ...
        'XMinorTick','on','YMinorTick','on','Position',[0.191107069608701,0.151302987340527,0.718640222810071,0.712761615620442]);
    set(gcf,'Color','w','Position',[703,199,831,743]);
    if ~verbose
        close all;
    end
end