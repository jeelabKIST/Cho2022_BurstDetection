%% Configure Library Paths
util_path = genpath('/Users/scho/Neuroscience_KIST/Cho2022_BurstDetection/utils');
data_path = genpath('/Users/scho/Neuroscience_KIST/Cho2022_BurstDetection/data/experimental_data/ESCAPE');
addpath(util_path);
addpath(data_path);
%% Load Data
metric_names = {'Precision','Recall','F1-Score','Temporal Concurrence','Detection Confidence'};
method_names = {'BP','ENV','S-STFT','MTP','CWT'};
nMetric = length(metric_names); % number of metrics
nMethod = length(method_names); % number of algorithms
% [1] Get LFP Signal Data
load('lfps/escape_eeg_all.mat');
anatomy_idx = 1; % 1-BLA; 2-PFC; 3-HPC
Fs = srate;      % sampling rate
Nyq = Fs/2;      % Nyquist frequency
tvec = times;    % time vector
signals = squeeze(data(anatomy_idx,:,:));
nTotTrial = size(data,3);
% [2] Get Human Annotations
[HUMAN_ANNOT_THETA, ct_theta] = load_annotations('theta',nTotTrial,tvec,util_path,data_path);
[HUMAN_ANNOT_BETA, ct_beta] = load_annotations('beta',nTotTrial,tvec,util_path,data_path);
[HUMAN_ANNOT_GAMMA, ct_gamma] = load_annotations('gamma',nTotTrial,tvec,util_path,data_path);
if all(ct_beta == ct_gamma) && all(ct_theta == ct_gamma)
    completed_trials = ct_gamma;
else
    error('Datasets should have a same number of trials.');
end
% NOTE: Selected trials should be 1st trials of Day 1 & Day 2 for 8 mice.
nComplete = length(completed_trials);
% [3] Grab Helper Functions
hp = helper;
%% Compute Burst Durations of Entire Dataset
burst_duration_theta = get_burst_durations(HUMAN_ANNOT_THETA);
burst_duration_beta = get_burst_durations(HUMAN_ANNOT_BETA);
burst_duration_gamma = get_burst_durations(HUMAN_ANNOT_GAMMA);
% NOTE: We do not convert the burst durations (ms) to cycles, as the
% conversion of such type will make the variable discrete.
burst_durations = hp.pack_data(burst_duration_theta,burst_duration_beta,burst_duration_gamma);
%% Compute Burst SNRs of Entire Dataset
burst_snr_theta = get_burst_snr(HUMAN_ANNOT_THETA,tvec,signals,Fs,completed_trials,8,10);
burst_snr_beta = get_burst_snr(HUMAN_ANNOT_BETA,tvec,signals,Fs,completed_trials,20,30);
burst_snr_gamma = get_burst_snr(HUMAN_ANNOT_GAMMA,tvec,signals,Fs,completed_trials,35,45);
% NOTE: We do not round up the burst SNRs (dB) to integers, as the
% conversion of such type will make the variable discrete.
burst_snrs = hp.pack_data(burst_snr_theta,burst_snr_beta,burst_snr_gamma);
%% Save Data
file_names = {'burst_durations.mat', 'burst_snrs.mat'};
var_names = {'burst_durations', 'burst_snrs'};
for i = 1:length(file_names)
    if ~isfile(file_names{i})
        save(file_names{i},var_names{i});
    else
        warning(['"' file_names{i} '" already exists. Saving will be skipped.']);
        continue;
    end
end
%% Plot Joint Distributions of Burst Properties
plot_burst_properties(HUMAN_ANNOT_THETA,burst_duration_theta,burst_snr_theta);
plot_burst_properties(HUMAN_ANNOT_BETA,burst_duration_beta,burst_snr_beta);
plot_burst_properties(HUMAN_ANNOT_GAMMA,burst_duration_gamma,burst_snr_gamma);
%% Plot Probability Density Estimates of Burst Properties
% [1] Set Visualization Parameters
n_burst_theta = length(burst_duration_theta);
n_burst_beta = length(burst_duration_beta);
n_burst_gamma = length(burst_duration_gamma);
band_names = {
    ['Theta (n=' num2str(n_burst_theta) ')'], ...
    ['Beta (n=' num2str(n_burst_beta) ')'], ...
    ['Gamma (n=' num2str(n_burst_gamma) ')']
};
cmap = {'#459395','#EB7C69','#FDA638'};
% [2] Plot Distributions
plot_distributions(burst_durations,"Burst Durations (ms)",band_names,cmap,'small');
plot_distributions(burst_snrs,"Burst SNRs (dB)",band_names,cmap,'small');
%% Check Assumptions Prior To Statistical Tests
% [1] Normality
fprintf("\n* Checking normality for burst durations ...\n");
check_normality(burst_duration_theta,burst_duration_beta,burst_duration_gamma);
fprintf("\n* Chekcing normality for burst SNRs ...\n");
check_normality(burst_snr_theta,burst_snr_beta,burst_snr_gamma);
% [2] Equal Variances
% NOTE: If the distributions are normal, you can use the Bartlett's test.
%       The Levene's test is less susceptible to the violations of normality.
fprintf("\n* Checking homoscedasticity for burst durations ...\n");
check_equal_variance({burst_duration_theta,burst_duration_beta,burst_duration_gamma},'LeveneQuadratic');
fprintf("\n* Checking homoscedasticity for burst SNRs ...\n");
check_equal_variance({burst_snr_theta,burst_snr_beta,burst_snr_gamma},'LeveneQuadratic');
% [3] Distribution Shapes
fprintf("\n* Checking distribution shapes for burst durations ...\n");
check_distribution_shape({burst_duration_theta,burst_duration_beta,burst_duration_gamma},true,false);
fprintf("\n* Checking distribution shapes for burst SNRs ...\n");
check_distribution_shape({burst_snr_theta,burst_snr_beta,burst_snr_gamma},true,false);
%% Statistical Tests

% As the data distributions are mostly non-normal and have unequal
% variances, the Kruskal-Wallis (KW) with Dunn's test was used to compare
% their statistical differences. A possible alternative would be the
% pair-wise Mann-Witney tests with Bonferroni correction.

% Although the data distributions have unequal variances, they have a 
% similar shape of distribution according to the Kolmogorov-Smirnov test.
% Based on the observation of the raw data, the outliers were observed to 
% be coming naturally from the bursts.

% MATLAB does not provide the Dunn's test. For a reliable statistical test,
% R was implemented. Please see `fig5D_burst_stat_test.r`.
%% Appendix: Accessory Functions
% Function #1: Load Annotations from a Designated Excel Sheet
function [DATA, completed_trials] = load_annotations(band_name,nTrial,tvec,util_path,data_path)
    sheet_name = [upper(band_name(1)) band_name(2:end)];
    DATA = read_annotation_file(sheet_name,nTrial,tvec,util_path,data_path);
    completed_trials = cellfun(@(name) str2double(name), regexp(fieldnames(DATA),'\d*','Match')'); % trials with human annotations
end

% Function #2: Extract Burst Durations from the Entire Dataset
function [burst_duration] = get_burst_durations(DATA)
    DATA = struct2cell(DATA);
    burst_times = cat(1,DATA{:});
    burst_duration = diff(burst_times,1,2) * 1000; % unit: ms
end

% Function #3: Extract Burst SNRs from the Entire Dataset
function [burst_snr] = get_burst_snr(DATA, tvec, signals, Fs, trial_idx, lo_f, hi_f)
    nData = length(fieldnames(DATA)); % number of trials available
    snrMat = cell(nData,1);
    k = 1;
    for t = trial_idx % for each trial
        signal = double(signals(:,t)');
        filt_sig = eegfilt(signal,Fs,lo_f,hi_f);
        human_obv = DATA.(['trial' num2str(t)]);
        amp_sig = zeros(size(human_obv,1),1);
        qnoise = ones(length(tvec),1);
        for i = 1:length(human_obv) % for each burst
            binterval = human_obv(i,:);
            bIdx = find(tvec >= binterval(1) & tvec <= binterval(end));
            amp_sig(i) = rms(filt_sig(bIdx));
            qnoise(bIdx) = 0;
        end
        qnoise(tvec < 60 | tvec > 120) = 0;
        amp_noise = rms(filt_sig(logical(qnoise)));
        snrMat{k} = 10*log10((amp_sig/amp_noise).^2);
        k = k + 1;
    end
    burst_snr = cat(1,snrMat{:});
end

% Function #4: Visualize Scatter Plot of Burst Property Distributions
function plot_burst_properties(DATA, burst_duration, burst_snr)
    figure();
    % [1] Set Parameters
    nData = length(fieldnames(DATA));
    DATA = struct2cell(DATA);
    burst_counts = cellfun(@(x) size(x,1), DATA); % number of bursts per each trial
    % [2] Build Colormap that Assigns Same Color to Bursts from Same Trial
    colormap(parula(nData)); % for discrete colors
    cmap = linspace(1,nData,nData);
    colors = zeros(sum(burst_counts),1);
    cidx = horzcat([1;cumsum(cellfun(@(x) size(x,1), DATA(1:end-1))) + 1], [cumsum(cellfun(@(x) size(x,1), DATA))]);
    for i = 1:length(cidx)
        colors(cidx(i,1):cidx(i,2)) = cmap(i);
    end
    % [2] Visualize Burst Properties in a Scatter Plot
    scatter(burst_duration,burst_snr,80,colors,'filled');
    xlabel('Duration (ms)');
    ylabel('SNR (dB)');
    % (A) Colorbar
    cb = colorbar();
    cb.Label.String = 'Trials';
    % (B) Figure Settings
    axis square;
    set(cb,'YDir','reverse','TickDirection','out','LineWidth',1,'FontSize',14,'FontWeight','bold');
    set(gca,'TickDir','out','LineWidth',2,'FontSize',14,'FontWeight','bold');
    set(gcf,'Color','w');
end

% Function #5: Visualize Burst Property Distribution
function plot_distributions(dataset, x_lbl, lg_lbl, cmap, fig_size)
    % [1] Set Parameters
    if nargin < 5
        fig_size = 'large';
    end
    nData = length(dataset);
    hmin = min(vertcat(dataset{:}));
    hmax = max(vertcat(dataset{:}));
    range = hmax - hmin;
    % [2] Visualize Distributions
    figure(); hold on;
    lp = cell(1,nData);
    for n = 1:nData
        [pde,x] = ksdensity(dataset{n});
        lp{n} = plot(x,pde,'Color',cmap{n});
        area(x,pde,'FaceColor',cmap{n},'FaceAlpha',0.3);
    end
    xlim([hmin-range*0.1, hmax+range*0.1]);
    leg = legend([lp{:}],lg_lbl,'FontSize',14,'Location','best');
    % (A) Figure Settings
    set([lp{:}],'LineWidth',2);
    if strcmp(fig_size, 'large')
        set(leg,'FontSize',24,'FontWeight','bold');
        set(gca,'FontSize',35,'FontWeight','bold','TickDir','out','TickLength',[0.015, 0.025],'LineWidth',5, ...
            'Position',[0.1911,0.2150,0.7186,0.7128]);
        set(gcf,'Color','white','Position',[682,552,831,314]);
    else
        xlabel(x_lbl);
        ylabel('Probability Density');
        set(gca,'FontSize',14,'FontWeight','bold','TickDir','out','LineWidth',2);
        set(gcf,'Color','white');
    end
end