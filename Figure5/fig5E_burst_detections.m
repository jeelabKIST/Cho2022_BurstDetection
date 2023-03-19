%% Configure Library Paths
util_path = genpath('/Users/scho/Neuroscience_KIST/Cho2022_BurstDetection/utils');
data_path = genpath('/Users/scho/Neuroscience_KIST/Cho2022_BurstDetection/data/experimental_data/ESCAPE');
addpath(util_path);
addpath(data_path);
%% Load Data
method_names = {'BP','ENV','S-STFT','MTP','CWT'};
nMethod = length(method_names);
% [1] Get LFP Signal Data
load('lfps/escape_eeg_all.mat');
anatomy_idx = 1; % 1-BLA; 2-PFC; 3-HPC
Fs = srate;      % sampling rate
Nyq = Fs/2;      % Nyquist frequency
tvec = times;    % time vector
nTrial = size(data,3);
signals = squeeze(data(anatomy_idx,:,:));
%% Read Annotations
sheet_name = 'Gamma';
HUMAN_ANNOT = read_annotation_file(sheet_name,nTrial,tvec,util_path,data_path);
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
%% Preprocess Signal
trial_idx = 1;
signal = double(signals(:,trial_idx)');
% [1] Apply IIR Notch Filter
w0 = 60/Nyq;
q_factor = 35;
bw = w0/q_factor;
[b,a] = iirnotch(w0,bw);
signal = filtfilt(b,a,signal);
% [2] Filter Signal
filt_sig = eegfilt(signal,Fs,lo_f,hi_f);
%% Burst Detections
% [1] Bandpass Filtering
[btime_bp,burst_bp] = detect_burst_timeseries(Fs,tvec,filt_sig,lo_f,hi_f,thrA_factor);
% [2] Amplitude Envelope
[btime_ev,burst_ev,~] = detect_burst_ampenv(Fs,tvec,filt_sig,lo_f,hi_f,prcthrA);
% [3] Single-Tapered Short Time Fourier Transform
[Spec_f_stp,Spec_t_stp,Spec_stp,~] = spectrogram_stp(tvec,signal,Fs,freq_range,dt_stp,nOscCycle_stp,verbose,interp_opt);
[btime_stp,burst_stp,~] = detect_burst_spectrogram(Spec_f_stp,Spec_t_stp,Spec_stp,dt_stp,specthrA_factor,verbose);
minIdx_stp = find(abs(Spec_f_stp-mid_f) == min(abs(Spec_f_stp-mid_f)));
[btime_stp,burst_stp] = extract_single_freqband(btime_stp,burst_stp,minIdx_stp);
% [4] Multitaper Spectrogram
[Spec_f_mtp,Spec_t_mtp,Spec_mtp,~] = spectrogram_mtp(tvec,signal,Fs,nw,ntp,freq_range,dt_mtp,nOscCycle_mtp,interp_opt,plot_psd,check_opt,plot_tp,verbose);
[btime_mtp,burst_mtp,~] = detect_burst_spectrogram(Spec_f_mtp,Spec_t_mtp,Spec_mtp,dt_mtp,specthrA_factor,verbose);
minIdx_mtp = find(abs(Spec_f_mtp-mid_f) == min(abs(Spec_f_mtp-mid_f)));
[btime_mtp,burst_mtp] = extract_single_freqband(btime_mtp,burst_mtp,minIdx_mtp);
% [5] Continuous Wavelet Spectrogram
[wav_freq,wav_time,Spec_cwt,~] = spectrogram_cwt(tvec,signal,Fs,wname,scale_opt,verbose);
[btime_cwt,burst_cwt,~] = detect_burst_spectrogram(wav_freq,wav_time,Spec_cwt',dt_cwt,specthrA_factor,verbose);
minIdx_cwt = find(abs(wav_freq-mid_f) == min(abs(wav_freq-mid_f)));
[btime_cwt,burst_cwt] = extract_single_freqband(btime_cwt,burst_cwt,minIdx_cwt);
% [6] Parse Bursts (Experiment Stage 2)
btime_bp = parse_data(btime_bp);
btime_ev = parse_data(btime_ev);
btime_stp = parse_data(btime_stp);
btime_mtp = parse_data(btime_mtp);
btime_cwt = parse_data(btime_cwt);
btimeMat = {btime_bp,btime_ev,btime_stp,btime_mtp,btime_cwt};
% [7] Get Human Annotations
time_intervals = HUMAN_ANNOT.(['trial' num2str(trial_idx)]);
binary_trace = zeros(1,length(tvec));
for i = 1:size(time_intervals,1)
    [~,idx_start] = min(abs(tvec - time_intervals(i,1)));
    [~,idx_end] = min(abs(tvec - time_intervals(i,2)));
    binary_trace(idx_start:idx_end) = 1;
end
%% Plot Signals and Burst Masks
% [1] Parse Signals
start = find(tvec == 60);
stop = find(tvec == 120);
tvec_s2 = tvec(start:stop);
signal_s2 = signal(start:stop);
filt_sig_s2 = filt_sig(start:stop);
ground_truth = binary_trace(start:stop);
truth_idx = logical(ground_truth);
% [2] Set Visualization Parameters
okabe_ito = [[0,114,178];
             [213,94,0];
             [0,158,115];
             [204,121,167];
             [240,228,66]]./255;
grey = [97,97,97]./255;
red = [214,39,40]./255;
light_green = [157, 204, 191]./255;
scale_factor = 0.4;
y_factor = 0.23;
hmin = 78.62;
hmax = 80.9;
% [3] Visualize Plot
figure(); hold on;
plot(tvec_s2,signal_s2*scale_factor,'Color',grey,'LineWidth',2);
scatter(tvec_s2(truth_idx),zeros(1,length(tvec_s2(truth_idx))) + 0.15,100,red,'filled','Marker','s');
ylbl = cell(1,nMethod+1);
ylbl{1} = text(hmin-0.02,0,'Raw');
for n = 1:nMethod
    plot(tvec_s2, filt_sig_s2 - (n*y_factor),'Color','k','LineWidth',3);
    btime = btimeMat{n};
    for i = 1:length(btime)
        bstart_idx = find(tvec_s2 >= btime{i}(1),1,'first');
        bstop_idx = find(tvec_s2 <= btime{i}(end),1,'last');
        pb = plot(tvec_s2(bstart_idx:bstop_idx),filt_sig_s2(bstart_idx:bstop_idx)-(n*y_factor),'Color',okabe_ito(n,:));
        pb.LineWidth = 3;
    end
    ylbl{n+1} = text(hmin-0.03,-(n*y_factor),method_names{n});
end
line([80.48,80.68],[-1.3,-1.3],'LineWidth',6,'Color','k');
txt = text(80.7,-1.3,'200 ms');
xlim([hmin,hmax]);
xlabel('Time (s)');
% (A) Figure Settings
set([ylbl{:}],'FontSize',30,'FontName','Helvetica','FontWeight','bold','HorizontalAlignment','right');
set(txt,'FontSize',30,'FontName','Helvetica','FontWeight','bold','VerticalAlignment','middle');
set(gca,'Color','none','Box','off','TickDir','out','XColor','none','YColor','none');
set(gcf,'Color','w','Position',[145,150,1469,794]);
%% Zoom-In Visualization
method_idx = 3;
lw = 5;
figure(); hold on;
plot(tvec_s2,filt_sig_s2,'Color','k','LineWidth',lw);
for i = 1:length(btimeMat{method_idx})
    bstart_idx = find(tvec_s2 >= btimeMat{method_idx}{i}(1),1,'first');
    bstop_idx = find(tvec_s2 <= btimeMat{method_idx}{i}(end),1,'last');
    bdur = tvec_s2(bstop_idx)-tvec_s2(bstart_idx);
    rectangle('Position',[tvec_s2(bstart_idx),-0.05,bdur,0.1],'EdgeColor',light_green,'FaceColor',light_green);
    pb = plot(tvec_s2(bstart_idx:bstop_idx),filt_sig_s2(bstart_idx:bstop_idx),'Color',okabe_ito(3,:));
    pb.LineWidth = lw;
end
xlim([80.2,80.68]);
ylim([-0.07,0.07]);
set(gca,'Color','none','Box','off','TickDir','out','XColor','none','YColor','none');
set(gcf,'Color','w','Position',[430,555,947,356]);
%% Appendix: Accessory Functions
% Function #1: Parse Bursts by Time
function [btime] = parse_data(btime)
    idx = cellfun(@(b) b(1) >= 60 && b(end) <=120, btime);
    btime = btime(idx);
end