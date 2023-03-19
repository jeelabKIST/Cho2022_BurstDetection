%% Configure Library Paths
util_path = genpath('/Users/scho/Neuroscience_KIST/Cho2022_BurstDetection/utils');
data_path = genpath('/Users/scho/Neuroscience_KIST/Cho2022_BurstDetection/data');
addpath(util_path);
addpath(data_path);
%% Load Data
% [1] Get Theoretical Results
band_name = 'gamma';
% NOTE: For Table 2, set `band_name` as 'gamma'. For Table 2-1, use 'beta'.
HEATMAP = load(['HM_' band_name '.mat']).HEATMAP;
DC = load(['DC_' band_name '.mat']);
scoreMat = {
    struct2cell(HEATMAP.precision)';
    struct2cell(HEATMAP.recall)';
    struct2cell(HEATMAP.f1_score)';
    struct2cell(HEATMAP.concurrence)';
    DC.detection_confidence
};
metric_names = {'Precision','Recall','F1-Score','Temporal Concurrence','Detection Confidence'};
method_names = {'BP','ENV','S-STFT','MTP','CWT'};
nMetric = length(metric_names); % number of metrics
nMethod = length(method_names); % number of algorithms
% [2] Get LFP Signal Data
load('lfps/escape_eeg_all.mat');
anatomy_idx = 1;       % 1-BLA; 2-PFC; 3-HPC
Fs = srate;            % sampling rate
Nyq = Fs/2;            % Nyquist frequency
tvec = times;          % time vector
nTrial = size(data,3); % number of trials
signals = squeeze(data(anatomy_idx,:,:));
% [3] Get Human Annotations
sheet_name = [upper(band_name(1)) band_name(2:end)];
HUMAN_ANNOT = read_annotation_file(sheet_name,nTrial,tvec,util_path,data_path);
completed_trials = cellfun(@(name) str2double(name), regexp(fieldnames(HUMAN_ANNOT),'\d*','Match')'); % trials with human annotations
nComplete = length(completed_trials);
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
%% Compute Experiment Performance By Each Mouse
% [1] Create Empty Tables To Store Results
[statistics_bp,statistics_ev,statistics_stp,statistics_mtp,statistics_cwt] = deal(cell(1,nTrial));
% [2] Measure Performance Metrics
for n = completed_trials
    signal = double(signals(:,n)');
    % (A) Apply IIR Notch Filter
    w0 = 60/Nyq;
    q_factor = 35;
    bw = w0/q_factor;
    [b,a] = iirnotch(w0,bw);
    signal = filtfilt(b,a,signal);
    %% Get Time Interval & Lengths of Human Observed Bursts
    human_obv = HUMAN_ANNOT.(['trial' num2str(n)]);
    btime_true = cell(1,length(human_obv));
    btrue_cycl = zeros(1,length(human_obv));
    for i = 1:length(human_obv)
        btime_true{i} = human_obv(i,:);
        btrue_cycl(i) = round(diff(btime_true{i})/(1/mid_f));
        % note: estimation of cycles can also be done by finding local 
        % maxima in given time intervals
    end
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
    %% Parse Bursts (Stage 2)
    btime_bp = parse_data(btime_bp);
    btime_ev = parse_data(btime_ev);
    btime_stp = parse_data(btime_stp);
    btime_mtp = parse_data(btime_mtp);
    btime_cwt = parse_data(btime_cwt);
    %% Measure Burst Statistics
    [bstat_bp] = get_burst_statistics(btime_true,btrue_cycl,btime_bp,1/Fs,mid_f); 
    [bstat_ev] = get_burst_statistics(btime_true,btrue_cycl,btime_ev,1/Fs,mid_f);
    [bstat_stp] = get_burst_statistics(btime_true,btrue_cycl,btime_stp,dt_stp,mid_f);
    [bstat_mtp] = get_burst_statistics(btime_true,btrue_cycl,btime_mtp,dt_mtp,mid_f);
    [bstat_cwt] = get_burst_statistics(btime_true,btrue_cycl,btime_cwt,dt_cwt,mid_f);
    %% Store Results: Statistics
    statistics_bp{n} = bstat_bp;
    statistics_ev{n} = bstat_ev;
    statistics_stp{n} = bstat_stp;
    statistics_mtp{n} = bstat_mtp;
    statistics_cwt{n} = bstat_cwt;
end
%% Compute Performance Metrics
[PRECISION,RECALL,F1_SCORE,TC,DC] = deal(zeros(nComplete,nMethod));
statistics = {statistics_bp,statistics_ev,statistics_stp,statistics_mtp,statistics_cwt};
for m = 1:nMethod
    for c = 1:nComplete
        statistic = statistics{m};
        stat = statistic{completed_trials(c)};
        [p,r,f1] = extract_prf(stat);
        tc = extract_temporal_concurrence(stat);
        dc = (f1 * exp(tc)) / exp(1);
        PRECISION(c,m) = p;
        RECALL(c,m) = r;
        F1_SCORE(c,m) = f1;
        TC(c,m) = tc;
        DC(c,m) = dc;
    end
end
%% Performance Statistics: Burst Detections against Human Annotations
% [1] Define Relevant Lambda Functions
lmean = @(x) mean(x,'omitnan');
lstd = @(x) std(x,'omitnan');
lsum = @(x) sum(x,1,'omitnan');
lvar = @(x) var(x,'omitnan');
% Note: NaN values can occur when no bursts are detected by an algorithm.
% [2] Algorithmic Moments across Metrics
avgs = string([lmean(PRECISION); lmean(RECALL); lmean(F1_SCORE); lmean(TC); lmean(DC)]);
stds = string([lstd(PRECISION); lstd(RECALL); lstd(F1_SCORE); lstd(TC); lstd(DC)]);
symbols = repmat("+",[nMetric, nMethod]);
moments = strcat(avgs,symbols,stds);
mnt_tbl = create_table(moments,metric_names,method_names);
% [3] Manhattan Distances across Metrics
distances = [lsum(abs(1-PRECISION)); lsum(abs(1-RECALL)); lsum(abs(1-F1_SCORE)); lsum(abs(1-TC)); lsum(abs(1-DC))];
dist_tbl = create_table(distances,metric_names,method_names);
sod = sum(distances,1); % sum of distances
% [4] Algorithmic Variances across Metrics
variances = [lvar(PRECISION); lvar(RECALL); lvar(F1_SCORE); lvar(TC); lvar(DC)];
var_tbl = create_table(variances,metric_names,method_names);
sov = sum(variances,1); % sum of variances
%% (Optional) Additional Statistical Tests
% In this section, you can perform One-Way ANOVA with Tukey's HSD on
% each type of metrics (e.g., precision). You first test the null
% hypothesis (i.e., "all population means are equal") and then address 
% which specific pair of group means shows a statistical difference.
stat_metric = struct('precision', PRECISION); % set the name and values of the metric
show_result = true;
% [1] One-Way ANOVA
group = method_names;
[p,tbl,stats] = anova1(stat_metric.precision,group,'off');
% [2] Multiple Comparisons Test: Tukey's HSD
[c,m_est,~,gnames] = multcompare(stats,'CType','tukey-kramer','Display','off');
mc_tbl = array2table(c,'VariableNames',{'Group','Control Group','Lower Limit','Difference','Upper Limit','p-value'});
mc_tbl.('Group') = gnames(mc_tbl.('Group'));
mc_tbl.('Control Group') = gnames(mc_tbl.('Control Group'));
if show_result
    fprintf('*** One-Way ANOVA & Tukey HSD: %s \n\n', upper(string(fieldnames(stat_metric))));
    disp(mc_tbl);
end
%% (Optional) Plot Performance Bar Graphs
% [1] Set Visualization Parameters
okabe_ito = [[0,114,178];
             [213,94,0];
             [0,158,115];
             [204,121,167];
             [240,228,66]]./255;
% [2] Visualize Performance Graphs
plot_bar_graph(PRECISION,okabe_ito,'Precision');
plot_bar_graph(RECALL,okabe_ito,'Recall');
plot_bar_graph(F1_SCORE,okabe_ito,'F1-Score');
plot_bar_graph(TC,okabe_ito,'Temporal Concurrence');
plot_bar_graph(DC,okabe_ito,'Detection Confidence');
%% Appendix: In-Script Functions
% Function #1: Compute Precision, Recall, F1-Score
function [P,R,F1] = extract_prf(statistics)
    classifications = statistics(3,:);
    TP = sum(ismember(classifications,'tp'));
    FN = sum(ismember(classifications,'fn'));
    FP = sum(ismember(classifications,'fp'));
    P = TP/(TP+FP); % precision
    R = TP/(TP+FN); % recall
    F1 = (2*P*R)/(P+R); % F1-score
end

% Function #2: Compute Temporal Concurrence
function [tc] = extract_temporal_concurrence(statistics)
    tpIdx = strcmp(statistics(3,:),'tp');
    tc_list = cell2mat(statistics(5,tpIdx));
    tc = mean(tc_list);
end

% Function #3: Parse Bursts by Time
function [btime] = parse_data(btime)
    idx = cellfun(@(b) b(1) >= 60 && b(end) <=120, btime);
    btime = btime(idx);
end

% Function #4: Format Data into Table
function [tbl] = create_table(data, row_name, col_name)
    tbl = array2table(data);
    tbl.Properties.RowNames = row_name;
    tbl.Properties.VariableNames = col_name;
end

% Function #5: Plot Performance Bar Graphs
function plot_bar_graph(data_table,color_palette,metric_name,size_opt)
    if nargin < 4
        size_opt = 'large';
    end
    % Organize Data
    nData = size(data_table,1);
    nMethod = size(data_table,2);
    data_avg = mean(data_table,1,'omitnan');
    data_std = std(data_table,1,'omitnan');
    grey = [97,97,97]./255;
    vmax = round(max(data_table, [], 'all'),1)+0.1;
    if mod(vmax,0.2) ~= 0
        vmax = vmax + 0.1;
    end
    % Visualize Bar Graphs
    figure(); hold on;
    x_axis = 1:nMethod;
    b = bar(x_axis, data_avg);
    b.FaceColor = 'flat';
    for n = 1:nMethod
        b.CData(n,:) = color_palette(n,:);
    end
    s = cell(size(x_axis));
    for n = 1:nMethod
        rnd = unifrnd(-0.2,0.2,1,nData);
        s{n} = scatter(ones(1,nData)*n+rnd,data_table(:,n),120,grey,'filled');
        s{n}.MarkerFaceAlpha = 0.65;
        s{n}.LineWidth = 1.5;
    end
    er = cell(size(x_axis));
    for n = 1:nMethod
        er{n} = errorbar(x_axis(n),data_avg(n),data_std(n),data_std(n));
        er{n}.Color = color_palette(n,:);
        if n == nMethod
            er{n}.Color = [209,206,23]./255; % darker yellow for visibility
        end
    end
    xticks(x_axis);
    xlim([0.2,nMethod+0.8]);
    ylim([0,vmax]);
    if strcmp(size_opt,'small')
        yticks(linspace(0,vmax,3));
        xticklabels(cell(1,nMethod));
    else
        xticklabels({'BP','ENV','S-STFT','MTP','CWT'});
    end
    xlabel('Algorithms');
    ylabel(metric_name);
    % Figure Settings
    set(b,'EdgeColor','none','BaseValue',0.004,'ShowBaseLine','off','FaceAlpha',0.4);
    set([er{:}],'CapSize',25,'LineStyle','none','LineWidth',10);
    if strcmp(size_opt,'small')
        axis square;
        set(gca,'Box','off','LineWidth',6,'TickDir','out','TickLength',[0.015, 0.025],'XTickLabelRotation',0, ... 
            'FontSize',48,'FontWeight','bold','Position',[0.1344,0.1455,0.8391,0.8115]);
        set(gcf,'Color','w','Position',[344,397,707,580]);
    else
        set(gca,'Box','off','LineWidth',6,'TickDir','out','TickLength',[0.015, 0.025],'XTickLabelRotation',0, ... 
            'FontSize',45,'FontWeight','bold','Position',[0.1755,0.1568,0.7206,0.7746]);
        set(gcf,'Color','w','Position',[344,1,1097,976]);
    end
end
