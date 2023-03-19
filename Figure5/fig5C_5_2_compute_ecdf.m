%% Configure Library Paths
util_path = genpath('/Users/scho/Neuroscience_KIST/Cho2022_BurstDetection/utils');
data_path = genpath('/Users/scho/Neuroscience_KIST/Cho2022_BurstDetection/data/experimental_data/ESCAPE');
addpath(util_path);
addpath(data_path);
%% Load Data
% [1] Get Theoretical Results
band_name = 'gamma';
HEATMAP = load(['/Users/scho/Neuroscience_KIST/Cho2022_BurstDetection/data/simulation_data/HM_' band_name '.mat']).HEATMAP;
DET_CONF = load(['/Users/scho/Neuroscience_KIST/Cho2022_BurstDetection/data/simulation_data/DC_' band_name '.mat']);
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
if strcmp(band_name, 'gamma')
    lo_f = 35; % lower frequency bound
    hi_f = 45; % upper frequency bound
elseif strcmp(band_name, 'beta')
    lo_f = 20;
    hi_f = 30;
end
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
[AUC.precision,AUC.recall,AUC.f1_score,AUC.concurrence,AUC.detection_confidence] = deal(zeros(nTrial,nMethod));
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
    verbose_cdf = true;
    indicator = {'precision','recall','f1_score','concurrence','detection_confidence'};
    for m = 1:length(scoreMat)
        if n > 1
            verbose_cdf = false;
        end
        [bscoreMat,auc] = compute_ecdf_auc(n,metric_names{m},scoreMat{m},btimeMat,dtMat, ... 
                           tvec,signal,Fs,lo_f,hi_f,verbose_cdf);
        AUC.(indicator{m})(n,:) = auc;
    end
    % [4] Report Progress
    if mod(n,10) == 0
        fprintf(['Processed: ' num2str(n) '/' num2str(nTrial) ' trials \n']);
    end
end
elapsed_time = toc;
fprintf(['Computation complete. (computational time: ' num2str(elapsed_time) 's) \n']);
%% Save AUC Results
file_name = ['AUC_stage2_' band_name '.mat'];
if ~isfile(file_name)
    save(file_name,'AUC');
else
    error('Save Error: The specified file already exists.');
end
%% Appendix: Accessory Functions
% Function #1: Parse Bursts by Time
function [btime] = parse_data(btime)
    idx = cellfun(@(b) b(1) >= 60 && b(end) <=120, btime);
    btime = btime(idx);
end
