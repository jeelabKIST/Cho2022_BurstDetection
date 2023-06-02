%% Configure Library Paths
util_path = genpath('/Users/scho/Neuroscience_KIST/Cho2022_BurstDetection/utils');
data_path = genpath('/Users/scho/Neuroscience_KIST/Cho2022_BurstDetection/data/experimental_data/ESCAPE');
addpath(util_path);
addpath(data_path);
%% Load Data
% [1] Get Theoretical Results
band_name = 'theta';
HEATMAP = load(['/Users/scho/Neuroscience_KIST/Cho2022_BurstDetection/data/simulation_data/HM_' band_name '.mat']).HEATMAP;
DET_CONF = load(['/Users/scho/Neuroscience_KIST/Cho2022_BurstDetection/data/simulation_data/DC_' band_name '.mat']);
scoreMat = {
    struct2cell(HEATMAP.precision)';
    struct2cell(HEATMAP.recall)';
    struct2cell(HEATMAP.f1_score)';
    struct2cell(HEATMAP.concurrence)';
    DET_CONF.detection_confidence
};
metric_names = {'precision','recall','f1_score','concurrence','detection_confidence'};
method_names = {'BP','ENV','S-STFT','MTP','CWT'};
nMetric = length(metric_names);
nMethod = length(method_names);
% [2] Get LFP Signal Data
load('lfps/escape_eeg_all.mat');
anatomy_idx = 1; % 1-BLA; 2-PFC; 3-HPC
Fs = srate;      % sampling rate
Nyq = Fs/2;      % Nyquist frequency
tvec = times;    % time vector
nTotTrial = size(data,3);
sel_trials = sort(cat(2,1:8:nTotTrial,5:8:nTotTrial)); % 1st trials of Day 1 & Day 2 for 8 mice
signals = squeeze(data(anatomy_idx,:,:));
% [3] Get AUC Values
auc = load(['../Figure5/AUC_stage2_' band_name '.mat']).AUC;
inverse_auc = auc;
for n = 1:nMetric
    metric = metric_names{n};
    inverse_auc.(metric) = inverse_auc.(metric)(sel_trials,:);
    inverse_auc.(metric) = 1./inverse_auc.(metric);
end
% [4] Grab Helper Functions
hp = helper;
%% Set Hyperparameters For Each Algorithm
% [1] Time Domain
if strcmp(band_name, 'theta')
    lo_f = 8;  % lower frequency bound
    hi_f = 10; % upper frequency bound
elseif strcmp(band_name, 'beta')
    lo_f = 20;
    hi_f = 30;
elseif strcmp(band_name, 'gamma')
    lo_f = 35; 
    hi_f = 45; 
end
mid_f =  (lo_f + hi_f) / 2;
% [2] Time-Frequency Domain
freq_range = unique([5:2:15 15:5:30 30:10:70 70:15:100]); % frequencies of interest
dt = 0.01;                                                % time step for sliding windows
nOscCycle_stp = 6; nOscCycle_mtp = 6;                     % frequency-dependent window size
nw = 2; ntp = 3;                                          % DPSS slepian sequences
wname = 'amor'; scale_opt = 'default';                    % type of and scaling option for the wavelet
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
specthrA_stp = 1.8;      % for S-STFT method
specthrA_mtp = 1.8;      % for MTP method
specthrA_cwt = 1.8;      % for CWT methods
% [4] Window Time Steps
dt_stp = dt;
dt_mtp = dt;
dt_cwt = 1/Fs;
% [5] Alter Hyperparameters for Theta Band % EDIT
if strcmp(band_name, 'theta')
    ntp = 2;                 % number of DPSS tapers
    thrA_factor = 1.8;       % for bandpass filtering method
    prcthrA = 90;            % for envelope-based method
    specthrA_stp = 1.2;      % for S-STFT method
    specthrA_mtp = 1.3;      % for MTP method % EDIT
    specthrA_cwt = 1.0;      % for CWT method
    dt_mtp = dt_mtp / 2;     % time step for sliding windows
end
%% Extract Burst Metrics For Each Algorithm Per Trial
tic; fprintf('Mapping burst metrics ... \n');
MAPPING = struct();
[MAPPING.precision,MAPPING.recall,MAPPING.f1_score,MAPPING.concurrence,MAPPING.detection_confidence] = deal(cell(length(sel_trials),nMethod));
for t = 1:length(sel_trials)
    % [1] Get Single-Trial Signal
    signal = double(signals(:,sel_trials(t))');
    w0 = 60/Nyq;
    q_factor = 35;
    bw = w0/q_factor;
    [b,a] = iirnotch(w0,bw);
    signal = filtfilt(b,a,signal); % apply IIR notch filter
    % [2] Detect Bursts
    % (A) Bandpass Filtering
    filt_sig = eegfilt(signal,Fs,lo_f,hi_f);
    [btime_bp,burst_bp] = detect_burst_timeseries(Fs,tvec,filt_sig,lo_f,hi_f,thrA_factor);
    % (B) Amplitude Envelope
    [btime_ev,burst_ev,~] = detect_burst_ampenv(Fs,tvec,filt_sig,lo_f,hi_f,prcthrA);
    % (C) Single-Tapered Short Time Fourier Transform
    [Spec_f_stp,Spec_t_stp,Spec_stp,~] = spectrogram_stp(tvec,signal,Fs,freq_range,dt_stp,nOscCycle_stp,verbose,interp_opt);
    [btime_stp,burst_stp,~] = detect_burst_spectrogram(Spec_f_stp,Spec_t_stp,Spec_stp,dt_stp,specthrA_stp,verbose);
    minIdx_stp = find(abs(Spec_f_stp-mid_f) == min(abs(Spec_f_stp-mid_f)));
    [btime_stp,burst_stp] = extract_single_freqband(btime_stp,burst_stp,minIdx_stp);
    % (D) Multitaper Spectrogram
    [Spec_f_mtp,Spec_t_mtp,Spec_mtp,~] = spectrogram_mtp(tvec,signal,Fs,nw,ntp,freq_range,dt_mtp,nOscCycle_mtp,interp_opt,plot_psd,check_opt,plot_tp,verbose);
    [btime_mtp,burst_mtp,~] = detect_burst_spectrogram(Spec_f_mtp,Spec_t_mtp,Spec_mtp,dt_mtp,specthrA_mtp,verbose);
    minIdx_mtp = find(abs(Spec_f_mtp-mid_f) == min(abs(Spec_f_mtp-mid_f)));
    [btime_mtp,burst_mtp] = extract_single_freqband(btime_mtp,burst_mtp,minIdx_mtp);
    % (E) Continuous Wavelet Spectrogram
    [wav_freq,wav_time,Spec_cwt,~] = spectrogram_cwt(tvec,signal,Fs,wname,scale_opt,verbose);
    [btime_cwt,burst_cwt,~] = detect_burst_spectrogram(wav_freq,wav_time,Spec_cwt',dt_cwt,specthrA_cwt,verbose);
    minIdx_cwt = find(abs(wav_freq-mid_f) == min(abs(wav_freq-mid_f)));
    [btime_cwt,burst_cwt] = extract_single_freqband(btime_cwt,burst_cwt,minIdx_cwt);
    % [3] Parse Bursts (Experiment Stage 2)
    btime_bp = parse_data(btime_bp);
    btime_ev = parse_data(btime_ev);
    btime_stp = parse_data(btime_stp);
    btime_mtp = parse_data(btime_mtp);
    btime_cwt = parse_data(btime_cwt);
    % [4] Map Metric Values to Burst Detections
    btimeMat = hp.pack_data(btime_bp, btime_ev, btime_stp, btime_mtp, btime_cwt);
    dtMat = hp.pack_data(1/Fs, 1/Fs, dt_stp, dt_mtp, dt_cwt);
    for m = 1:length(scoreMat) % for each performance metric
        bscoreMat = map_burst_metric(scoreMat{m},btimeMat,dtMat,tvec,signal,Fs,lo_f,hi_f); % dimension: (nMethod)
        MAPPING.(metric_names{m})(t,:) = bscoreMat; % dimension: (nTrial, nMethod)
    end
    % [5] Report Progress
    if mod(t,4) == 0
        fprintf(['Processed: ' num2str(t) '/' num2str(length(sel_trials)) ' trials \n'])
    end
end
elapsed_time = toc;
fprintf(['Computation complete. (computational time: ' num2str(elapsed_time) 's) \n']);
%% Rank Algorithms By 1/AUC For Each Trial
RANK = struct();
for n = 1:nMetric
    metric = metric_names{n};
    inv_auc = inverse_auc.(metric);
    [~,sort_idx] = sort(inv_auc,2,'descend'); % sort by column (algorithms)
    RANK.(metric) = sort_idx;
end
%% Visualize Rank vs. Metric Plot
% [1] Rank Metric Mappings by 1/AUC
RANKED_MAPPING = MAPPING;
for m = 1:nMetric
    metric = metric_names{m};
    mapping = MAPPING.(metric);
    rank = RANK.(metric);
    ranked_mapping = cell(size(mapping));
    for r = 1:size(rank,1)
        ranked_mapping(r,:) = mapping(r,rank(r,:));
    end
    ranked_mapping = arrayfun(@(i) horzcat(ranked_mapping{:,i}),1:nMethod,'UniformOutput',false);
    RANKED_MAPPING.(metric) = ranked_mapping;
end
% [2] Save Ranked Mappings for Visualization
file_name = ['ranked_mapping_' band_name '.mat'];
if ~isfile(file_name)
    save(file_name,'RANKED_MAPPING');
else
    error('Save Error: The specified file already exists.')
end
% [3] Extract Ranked Detection Confidence Scores
ranked_metric = RANKED_MAPPING.detection_confidence;
% NOTE: Metric type can be changed for the specific use case.
mean_rm = cellfun(@(x) mean(x,'omitnan'),ranked_metric);
std_rm = cellfun(@(x) std(x,'omitnan'),ranked_metric);
% [4] Visualize Metric Values by Ranking
figure(); hold on;
x_axis = 1:nMethod;
b = bar(x_axis,mean_rm);
b.FaceColor = 'flat';
er = cell(size(x_axis)); % preallocate error bars
for n = 1:nMethod
    b.CData(n,:) = 'w';
    er{n} = errorbar(x_axis(n),mean_rm(n),std_rm(n),std_rm(n));
    er{n}.Color = 'k';
end
plot_scatter = false;
if plot_scatter == true
    red = [184,15,10]./255;
    s = cell(size(x_axis));
    for n = 1:nMethod
        rnd = unifrnd(-0.2,0.2,1,length(ranked_metric{n}));
        s{n} = scatter(ones(1,length(ranked_metric{n}))*n+rnd,ranked_metric{n},90,red,'filled');
        s{n}.MarkerFaceAlpha = 0.45;
        s{n}.LineWidth = 1.5;
    end
end
% (A) Axes Settings
xticks(x_axis);
yticks(0:0.2:1);
xlim([0.2,nMethod+0.8]);
ylim([0, 1.2]);
xlabel('Ranks');
ylabel('Detection Confidence');
% (B) Figure Settings
set(b,'EdgeColor','k','LineWidth',4,'BaseValue',-1e-4,'ShowBaseLine','off','FaceAlpha',0.3);
set([er{:}],'CapSize',20,'LineStyle','none','LineWidth',5);
set(gca,'FontSize',35,'LineWidth',5,'TickDir','out','TickLength',[0.015, 0.025],'XTickLabelRotation',0,'Box','off','FontWeight','bold',...
    'Position',[0.221648616102467,0.146554506483428,0.73623345369296,0.778445493516572]);
set(gcf,'Color','w','Position',[703,199,831,743]);
%% Plot Rank Distribution
% [1] Set Visualization Parameters
okabe_ito = [[0,114,178];
             [213,94,0];
             [0,158,115];
             [204,121,167];
             [240,228,66]]./255;
% [2] Visualize Distribution of Algorithms in Each Rank
titles = {'Precision','Recall','F1 Score','Temporal Concurrence','Deteciton Confidence'};
for m = 1:nMetric
    metric = metric_names{m};
    % (A) Count Each Algorithm Per Rank
    rank_counts = zeros(nMethod,nMethod);
    for n = 1:nMethod
        rank_counts(:,n) = histcounts(RANK.(metric)(:,n),1:nMethod+1);
    end
    % (B) Stacked Bar Plot
    figure();
    b = bar(rank_counts','stacked','FaceColor','flat');
    for c = 1:length(b)
        b(c).CData = okabe_ito(c,:);
    end
    ylim([0,length(sel_trials)]);
    yticks(0:8:length(sel_trials));
    xlabel('Rank');
    ylabel('Count');
    % (C) Figure Settings
    set(b,'EdgeColor','none','FaceAlpha',0.85);
    set(gca,'TickDir','out','FontSize',40,'FontWeight','bold','LineWidth',5,'Box','off');
    set(gcf,'Color','w','Position',[476,446,484,420]);
%     title(titles{m},'FontSize',14);
end
%% Appendix: Accessory Functions
% Function #1: Parse Bursts by Time
function [btime] = parse_data(btime)
    idx = cellfun(@(b) b(1) >= 60 && b(end) <=120, btime);
    btime = btime(idx);
end