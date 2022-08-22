%% Configure Library Paths
util_path = genpath('/Users/jeelab/Desktop/Cho2022_BurstDetection/utils');
data_path = genpath('/Users/jeelab/Desktop/Cho2022_BurstDetection/data/experimental_data');
addpath(util_path);
addpath(data_path);
%% Load LFP Signal Data
% [1] Preprocess Data
load('lfps/escape_eeg_all.mat');
anatomy_idx = 1; % 1-BLA; 2-PFC; 3-HPC
Fs = srate;      % sampling rate
Nyq = Fs/2;      % Nyquist frequency
tvec = times;    % time vector
signals = double(squeeze(data(anatomy_idx,:,:)));
% [2] Select Specific Trial
trialID = 57; % selected example (Mouse 8 x Day 1 x Session 1)
signal = signals(:,trialID)';
%% Set Hyperparameters For Each Algorithm
nMethod = 5;
% [1] Time Domain
f = 40;        % frequency of interest
lo_f = 35;  % lower frequency bound
hi_f = 45;  % upper frequency bound
mid_f =  (lo_f + hi_f) / 2;
% [2] Time-Frequency Domain
freq_range = unique([15:5:30 30:10:60]); % frequencies of interest
dt = 0.01;                               % time step for sliding windows
nOscCycle_stp = 6; nOscCycle_mtp = 6;    % frequency-dependent window size
nw = 2; ntp = 3;                         % DPSS slepian sequences
wname = 'amor'; scale_opt = 'default';   % type of and scaling option for the wavelet
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
%% Get Burst Time Intervals Using Each Algorithm (Burst Detections)
% [1] Bandpass Filtering
filt_sig = eegfilt(signal,Fs,lo_f,hi_f);
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
%% Get Binary Burst Time Traces
TRACE = struct();
algorithms = {'BP','ENV','STP','MTP','CWT'};
btime = {btime_bp,btime_ev,btime_stp,btime_mtp,btime_cwt};
for i = [1,2,5]
    data = btime{i};
    trace = zeros(1,length(tvec));
    for n = 1:length(data)
        trace = trace + ismember(tvec, data{n});
    end
    TRACE.(algorithms{i}) = trace;
end
for i = [3,4]
    data = btime{i};
    trace = zeros(1,length(tvec));
    for n = 1:length(data)
        bstart = data{n}(1);
        bend = data{n}(end);
        [~,idx_start] = min(abs(tvec - bstart));
        [~,idx_end] = min(abs(tvec - bend));
        trace(idx_start:idx_end) = 1;
    end
    TRACE.(algorithms{i}) = trace;
end
%% Get Ground Truth Bursts (Human Observations)
bstart_human = [85.8779, 86.3662];
bend_human = [85.9893, 86.4785];
nBurst = length(bstart_human);
trace = zeros(1,length(tvec));
for i = 1:nBurst
    [~,idx_start] = min(abs(tvec - bstart_human(i)));
    [~,idx_end] = min(abs(tvec - bend_human(i)));
    trace(idx_start:idx_end) = 1;
end
TRACE.HUMAN = trace;
%% Visualize Time Traces
% [1] Set Visualization Parameters
okabe_ito = [[0,114,178];
             [213,94,0];
             [0,158,115];
             [204,121,167];
             [240,228,66]]./255;
red = [219,77,77]./255;
grey = [81,81,81]./255;
hmin = 85.65; hmax = 86.75;
scale = 15;
% [2] Plot Burst Time Traces
figure(); hold on;
plot(tvec,filt_sig*scale+8.6,'Color',[0,0,0,0.8],'LineWidth',3.6);
plot(tvec,TRACE.HUMAN+8.6,'Color',red,'LineWidth',4);
k = nMethod+1;
for n = 1:nMethod
    plot(tvec,TRACE.(algorithms{n}) + k,'Color',okabe_ito(n,:),'LineWidth',4);
    k = k - 1.5;
end
xline(85.9326,'--','LineWidth',3.4);
line([hmax,hmax],[0.00*scale+7.2, 0.05*scale+7.2],'LineWidth',3.0,'Color',grey);
line([hmax-0.15,hmax],[7.2,7.2],'LineWidth',3.0,'Color',grey);
txt1 = text(hmax+0.01,7.58,'50µV');
txt2 = text(hmax-0.15,6.95,'150ms');
xlim([hmin, hmax]);
ylim([-0.5, 10]);
yticks([0:1.5:6, 8.5]);
yticklabels([fliplr(algorithms), {'Human'}]);
ytickangle(90);
xlabel('Time (s)');
title({'Example Time Traces of', 'Detected Bursts'});
set([txt1, txt2],'FontSize',25,'FontName','Helvetica','FontWeight','bold');
set(gca,'TickDir','out','TickLength',[0.015, 0.025],'LineWidth',4,'FontWeight','bold','FontSize',32,'Position',[0.163005780346821,0.137434554973822,0.683236994219653,0.738219895287958]);
set(gcf,'Color','w','Position',[481,112,865,764]);
%% Visualize Temporal Concurrences
% [1] Compute Temporal Concurrences
tstart = zeros(nMethod,nBurst);
tend = zeros(nMethod,nBurst);
for n = 1:nMethod
    data = btime{n};
    k = 1;
    for i = 1:length(data)
        bstart = data{i}(1);
        bend = data{i}(end);
        if bstart > hmin && bend < hmax
            if bstart - mean(tstart(1:n-1,k)) > 0.3
                k = k + 1;
            end
            tstart(n,k) = bstart;
            tend(n,k) = bend;
            k = k + 1;
        end
    end
end
temporal_concurrence = zeros(nMethod,nBurst);
for n = 1:nMethod
    for i = 1:nBurst
        nmr = min([bend_human(i), tend(n,i)]) - max([bstart_human(i), tstart(n,i)]);
        den = max([bend_human(i), tend(n,i)]) - min([bstart_human(i), tstart(n,i)]);
        temporal_concurrence(n,i) = nmr/den;
        if tstart(n,i) == 0 && tend(n,i) == 0
            temporal_concurrence(n,i) = 0;
        end
    end
end
temporal_concurrence = temporal_concurrence .* 100; % convert to percentage scale
temporal_concurrence = rot90(temporal_concurrence,2); % rearrange values
% [2] Set Visualization Parameters
alpha_list = [0.6,0.4];
% [3] Plot Temporal Concurrences
figure();
bh = barh(temporal_concurrence);
for i = 1:length(bh)
    bh(i).FaceColor = 'k';
    bh(i).FaceAlpha = alpha_list(i);
end
xticks(0:20:100);
yticks(1:nMethod);
yticklabels(fliplr(algorithms));
ytickangle(90);
xlabel('Temporal Concurrence');
ylabel('Algorithms');
title({'Single-Burst', 'Temporal Concurrence'});
% (A) Figure Settings
leg = legend([bh(2) bh(1)], {'Burst #1','Burst #2'});
set(leg,'Position',[0.7555,0.7559,0.1803,0.1001]);
set(bh,'EdgeColor','none','BaseValue',0.004,'ShowBaseLine','off');
set(gca,'Box','off','TickDir','out','TickLength',[0.015, 0.025],'LineWidth',4,'FontWeight','bold','FontSize',32, ... 
    'Position',[0.163005780346821,0.137434554973822,0.683236994219653,0.738219895287958]);
set(gcf,'Color','w','Position',[481,112,865,764]);