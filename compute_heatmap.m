%% Configure Library Path
util_path = genpath('/Users/jeelab/Desktop/Cho2022_BurstDetection/utils');
data_path = genpath('/Users/jeelab/Desktop/Cho2022_BurstDetection/data/simulation_data');
addpath(util_path);
addpath(data_path);
%% Set System Parameters
filepath_seed = fullfile(pwd,'data','simulation_data','randseed_beta.mat');
% To reproduce a previous result, save and call a set of random seeds that 
% was used to simulate the signals. The .mat file should be located in the
% directory corresponding to `data_path` above.
save_result = true;
filename_result = 'HM_beta.mat'; % file name for the results to be saved
filepath_result = fullfile(pwd,'data','simulation_data',filename_result); % path where results will be saved
%% Theoretical Simulation of the Signals
fprintf('Starting Signal Simulation ... \n');
% [1] Set Fixed Parameters
method = 'tukey';
T = 300;             % unit: s
f1 = 25; f2 = 40;    % unit: Hz
nBurst1 = 30;
nBurst2 = 30;
rng_pks = [3 12];
sim_verbose = false;
% [2] Set Varying Parameters
nSimulationPerNoise = 1000; % number of simulations per each SNR
listNoise = -10:2:10;
if isfile(filepath_seed)
    setSeed = load(filepath_seed).randseed;
    fprintf('*** Note: Random seeds imported for reproducibility. \n');
else
    setSeed = zeros(length(listNoise),nSimulationPerNoise);
    for n = 1:length(listNoise)
        setSeed(n,:) = randsample(1:2^15,nSimulationPerNoise);
    end
end
% [3] Preallocate Dataset
sz_dataset = cell(length(listNoise),nSimulationPerNoise);
[ds_signal,ds_burst,ds_length] = deal(sz_dataset);
% [4] Get Simulation Data
tic;
for n = 1:length(listNoise)
    % Set SNR (dB) Level
    SNR_wn = listNoise(n);
    for s = 1:size(setSeed,2)
        % Set Random Seed
        rnd_seed = setSeed(n,s);
        [~,simulated_signal,Fs,~,btime_true,burst_true,cyc_true] = generate_simulation_rnd(method,T,f1,f2,nBurst1,nBurst2,SNR_wn,rng_pks,rnd_seed,sim_verbose);
        ds_signal{n,s} = simulated_signal;
        ds_burst{n,s} = cat(2,btime_true,burst_true);
        ds_length{n,s} = cyc_true;
    end
    fprintf(['Progress: running ... [' num2str(n*s) '/' num2str(length(listNoise)*nSimulationPerNoise) '] \n']);
end
tvec = linspace(0,T,T*Fs)'; % time vector
Ps = 1/Fs; % sampling period
elapsed_time = toc;
clc; fprintf(['STEP1: Simulation complete! (computational time: ' num2str(elapsed_time) 's) \n']);
%% Set Parameters
fprintf('Starting parameter setting ... \n'); tic;
% [1] Time Domain: Bandpass Filtering
lf = f1-2; hf = f1+2;                                    % boundary of frequencies for filtering
% [2] Time-Frequency Domain: Power Spectrograms
freq_range = unique([15:5:30 30:10:60]);                 % frequencies of interest
dt = 0.01;                                               % time step for sliding windows
nOscCycle_stp = 6; nOscCycle_mtp = 6;                    % frequency-dependent window size
nw = 2; ntp = 3;                                         % DPSS slepian sequences
wname = 'amor'; scale_opt = 'default';                   % type of and scaling option for the wavelet
verbose = false; interp_opt = false;                     % visualization options
plot_psd = false; check_opt = true; plot_tp = false;     % visualization and optimization options
% [3] Heatmaps: Precision, Recall, F1-Score, and Temporal Accuracy
% (A) Amplitude Thresholds
thrA_factor = 2.5;                % for bandpass filtering method
prcthrA = 95;                     % for envelope-based method
specthrA_factor = 1.8;            % for spectral analysis methods
% (B) Window Time Step Used
dt_stp = dt;
dt_mtp = dt;
dt_cwt = Ps;
% For cwt.m, translation is fixed to 1.
% References: https://www.mathworks.com/matlabcentral/answers/283547-continuous-wavelet-transform-implementation-using-morlet?s_tid=srchtitle
%             https://www.mathworks.com/help/wavelet/gs/interpreting-continuous-wavelet-coefficients.html
%             https://www.mathworks.com/help/wavelet/gs/continuous-and-discrete-wavelet-transforms.html
% (C) Frequency and Time Vectors for Each Method
% Note: These vectors would be equivalent for every spectrogram, so it is
% unnecessary to output them every time.
input_signal = ds_signal{1,1};
[Spec_f_stp,Spec_t_stp,~,fft_interval_stp] = spectrogram_stp(tvec,input_signal,Fs,freq_range,dt,nOscCycle_stp,verbose,interp_opt);
[Spec_f_mtp,Spec_t_mtp,~,fft_interval_mtp] = spectrogram_mtp(tvec,input_signal,Fs,nw,ntp,freq_range,dt,nOscCycle_mtp,interp_opt,plot_psd,check_opt,plot_tp,verbose);
[wav_freq,wav_time,~,coi] = spectrogram_cwt(tvec,input_signal,Fs,wname,scale_opt,verbose);
% (D) Frequency of Interest
fft_dist_stp = abs(fft_interval_stp(:,2)-f1) + abs(fft_interval_stp(:,1)-f1);
fft_dist_mtp = abs(fft_interval_mtp(:,2)-f1) + abs(fft_interval_mtp(:,1)-f1);
[~,minIdx_stp] = min(fft_dist_stp);
[~,minIdx_mtp] = min(fft_dist_mtp);
parse_start = find(wav_freq<freq_range(end),1,'first');
parse_end  = find(wav_freq>freq_range(1),1,'last');
wav_freq = wav_freq(parse_start:parse_end);
[~,minIdx_cwt] = min(abs(wav_freq-f1));
% [4] Preallocate Heatmaps
listCycle = rng_pks(1):rng_pks(end);
[P_bp,P_ev,P_stp,P_mtp,P_cwt] = deal(zeros(length(listNoise),length(listCycle))); % precision
[R_bp,R_ev,R_stp,R_mtp,R_cwt] = deal(zeros(length(listNoise),length(listCycle))); % recall (sensitivity)
[F_bp,F_ev,F_stp,F_mtp,F_cwt] = deal(zeros(length(listNoise),length(listCycle))); % F1-score
[T_bp,T_ev,T_stp,T_mtp,T_cwt] = deal(zeros(length(listNoise),length(listCycle))); % temporal concurrence
% [5] Report Computational Time
elapsed_time = toc;
clc; fprintf(['STEP2: Parameters set! (computational time: ' num2str(elapsed_time) 's) \n']);
%% Construct Heatmaps
fprintf('Starting heatmap construction ... \n'); tic;
for rowIdx = 1:size(ds_signal,1)
    % [1] Preallocate Burst Statistics Outputs
    [statistics_bp, statistics_ev, statistics_stp, statistics_mtp, statistics_cwt] = deal(cell(1,size(ds_signal,2)));
    for colIdx = 1:size(ds_signal,2)
        %% Preprocessing of the Simulated Signals
        % [2] Time Domain: Bandpass Filtering
        filtered_signal = eegfilt(ds_signal{rowIdx,colIdx}',Fs,lf,hf);
        % [3] Time-Frequency Domain: Power Spectrograms
        input_signal = ds_signal{rowIdx,colIdx};
        % (A) Spectrogram: Single-Tapered Short Time Fourier Transform
        [~,~,Spec_stp,~] = spectrogram_stp(tvec,input_signal,Fs,freq_range,dt,nOscCycle_stp,verbose,interp_opt);
        % (B) Spectrogram: Multitaper
        [~,~,Spec_mtp,~] = spectrogram_mtp(tvec,input_signal,Fs,nw,ntp,freq_range,dt,nOscCycle_mtp,interp_opt,plot_psd,check_opt,plot_tp,verbose);
        % (C) Spectrogram: Continuous Wavelet Transform
        [~,~,Spec_cwt,~] = spectrogram_cwt(tvec,input_signal,Fs,wname,scale_opt,verbose);
        % (D) Parse CWT Spectrogram for Memory
        Spec_cwt = Spec_cwt(parse_start:parse_end,:); % note that dimension is transposed for cwt
        %% Burst Detection & Statistics
        % [4] Detect Bursts
        % (A) Get True Simulated Bursts
        btime_true = ds_burst{rowIdx,colIdx}(:,1);
        burst_true = ds_burst{rowIdx,colIdx}(:,2);
        btrue_cycl = ds_length{rowIdx,colIdx}';
        % (B) Detect Bursts: Time Domain
        [btime_bp,burst_bp] = detect_burst_timeseries(Fs,tvec,filtered_signal,lf,hf,thrA_factor);
        [btime_ev,burst_ev,benvelope] = detect_burst_ampenv(Fs,tvec,filtered_signal,lf,hf,prcthrA);
        % (C) Detect Bursts: Time-Frequency Domain
        [btime_stp,burst_stp] = detect_burst_spectrogram(Spec_f_stp,Spec_t_stp,Spec_stp,dt_stp,specthrA_factor,verbose);
        [btime_mtp,burst_mtp] = detect_burst_spectrogram(Spec_f_mtp,Spec_t_mtp,Spec_mtp,dt_mtp,specthrA_factor,verbose);
        [btime_cwt,burst_cwt] = detect_burst_spectrogram(wav_freq,wav_time,Spec_cwt',dt_cwt,specthrA_factor,verbose);
        [btime_stp,burst_stp] = extract_single_freqband(btime_stp,burst_stp,minIdx_stp);
        [btime_mtp,burst_mtp] = extract_single_freqband(btime_mtp,burst_mtp,minIdx_mtp);
        [btime_cwt,burst_cwt] = extract_single_freqband(btime_cwt,burst_cwt,minIdx_cwt);
        % (D) Compute Burst Statistics
        [bstat_bp] = get_burst_statistics(btime_true,btrue_cycl,btime_bp,Ps,f1);
        [bstat_ev] = get_burst_statistics(btime_true,btrue_cycl,btime_ev,Ps,f1);
        [bstat_stp] = get_burst_statistics(btime_true,btrue_cycl,btime_stp,dt_stp,f1);
        [bstat_mtp] = get_burst_statistics(btime_true,btrue_cycl,btime_mtp,dt_mtp,f1);
        [bstat_cwt] = get_burst_statistics(btime_true,btrue_cycl,btime_cwt,dt_cwt,f1);
        % (E) Store Burst Statistics
        statistics_bp{colIdx} = bstat_bp;
        statistics_ev{colIdx} = bstat_ev;
        statistics_stp{colIdx} = bstat_stp;
        statistics_mtp{colIdx} = bstat_mtp;
        statistics_cwt{colIdx} = bstat_cwt;
        % [5] Report Progress
        if mod(colIdx,round(size(ds_signal,2)/10)) == 0
            fprintf(['Progress: running ... [SNR: ' num2str(listNoise(rowIdx)) ' dB and trial #: ' num2str(colIdx) ' trials] \n']);
            disp(P_bp);
            disp(P_ev);
        end
    end
    %% Precision, Recall, F1-Score & Temporal Concurrence Heatmap
    % [6] Compute Precision, Recall, and F1-Score for Each SNR Value
    %     (1000 trials per one SNR value)
    prf_per_snr_bp = extract_precision_recall(statistics_bp,listCycle);
    prf_per_snr_ev = extract_precision_recall(statistics_ev,listCycle);
    prf_per_snr_stp = extract_precision_recall(statistics_stp,listCycle);
    prf_per_snr_mtp = extract_precision_recall(statistics_mtp,listCycle);
    prf_per_snr_cwt = extract_precision_recall(statistics_cwt,listCycle);
    % [7] Construct Heatmap
    P_bp(rowIdx,:) = prf_per_snr_bp(1,:);   P_ev(rowIdx,:) = prf_per_snr_ev(1,:);   P_stp(rowIdx,:) = prf_per_snr_stp(1,:);   P_mtp(rowIdx,:) = prf_per_snr_mtp(1,:);   P_cwt(rowIdx,:) = prf_per_snr_cwt(1,:);
    R_bp(rowIdx,:) = prf_per_snr_bp(2,:);   R_ev(rowIdx,:) = prf_per_snr_ev(2,:);   R_stp(rowIdx,:) = prf_per_snr_stp(2,:);   R_mtp(rowIdx,:) = prf_per_snr_mtp(2,:);   R_cwt(rowIdx,:) = prf_per_snr_cwt(2,:);
    F_bp(rowIdx,:) = prf_per_snr_bp(3,:);   F_ev(rowIdx,:) = prf_per_snr_ev(3,:);   F_stp(rowIdx,:) = prf_per_snr_stp(3,:);   F_mtp(rowIdx,:) = prf_per_snr_mtp(3,:);   F_cwt(rowIdx,:) = prf_per_snr_cwt(3,:);
    % [8] Compute Temporal Concurrence for Each SNR Value
    tc_per_snr_bp = extract_temporal_concurrence(statistics_bp,listCycle);
    tc_per_snr_ev = extract_temporal_concurrence(statistics_ev,listCycle);
    tc_per_snr_stp = extract_temporal_concurrence(statistics_stp,listCycle);
    tc_per_snr_mtp = extract_temporal_concurrence(statistics_mtp,listCycle);
    tc_per_snr_cwt = extract_temporal_concurrence(statistics_cwt,listCycle);
    % [9] Construct Heatmap
    T_bp(rowIdx,:) = tc_per_snr_bp;   T_ev(rowIdx,:) = tc_per_snr_ev;   T_stp(rowIdx,:) = tc_per_snr_stp;   T_mtp(rowIdx,:) = tc_per_snr_mtp;   T_cwt(rowIdx,:) = tc_per_snr_cwt;
end
% [10] Report Computing Time
elapsed_hmp = toc;
clc; fprintf(['STEP3: Heatmaps complete! (computational time: ' num2str(elapsed_hmp) 's) \n']);
%% Store Results in Structure
HEATMAP = struct();
HEATMAP.randseed = setSeed;
HEATMAP.precision.bp = P_bp;    HEATMAP.recall.bp = R_bp;    HEATMAP.f1_score.bp = F_bp;    HEATMAP.concurrence.bp = T_bp;
HEATMAP.precision.ev = P_ev;    HEATMAP.recall.ev = R_ev;    HEATMAP.f1_score.ev = F_ev;    HEATMAP.concurrence.ev = T_ev;
HEATMAP.precision.stp = P_stp;  HEATMAP.recall.stp = R_stp;  HEATMAP.f1_score.stp = F_stp;  HEATMAP.concurrence.stp = T_stp;
HEATMAP.precision.mtp = P_mtp;  HEATMAP.recall.mtp = R_mtp;  HEATMAP.f1_score.mtp = F_mtp;  HEATMAP.concurrence.mtp = T_mtp;
HEATMAP.precision.cwt = P_cwt;  HEATMAP.recall.cwt = R_cwt;  HEATMAP.f1_score.cwt = F_cwt;  HEATMAP.concurrence.cwt = T_cwt;
%% Save Results
if save_result
    if isfile(filepath_result)
        warning('Save attempt failed. A file with the specified name already exists.');
    else
        save(filepath_result,'HEATMAP');
    end
end
%% Appendix: In-Script Funcitons
% Function #1: Compute Precision, Recall, and F1-Score
function [prf_per_snr] = extract_precision_recall(statistics,listCycle)
    statistics = [statistics{:}];
    statistics = sortrows(statistics(3:4,:)',2);
    prf_per_snr = zeros(3,length(listCycle));
    uniqueCycle = unique(cell2mat(statistics(:,2)));
    if sum(uniqueCycle < 3) ~= 0 || sum(uniqueCycle > 12) ~= 0
        out_of_range = uniqueCycle(~ismember(uniqueCycle,listCycle));
        for i = 1:length(out_of_range)
            c = out_of_range(i);
            idx = cellfun(@(x) x==c, statistics(:,2));
            if c < listCycle(1)
                error('ComputationError: Some detected bursts did not satisfy the duration threshold.');
            end
            if c > listCycle(end)
                statistics(idx,2) = {listCycle(end)};
            end
        end
    end
    for i = 1:length(listCycle)
        c = listCycle(i);
        idx = cellfun(@(x) x==c, statistics(:,2));
        stat_per_cyc = statistics(idx,1);
        TP = sum(ismember(stat_per_cyc,'tp'));
        FP = sum(ismember(stat_per_cyc,'fp'));
        FN = sum(ismember(stat_per_cyc,'fn'));
        P = TP/(TP+FP);    % precision
        R = TP/(TP+FN);    % recall
        F = (2*P*R)/(P+R); % F1-score
        prf_per_snr(1,i) = P;
        prf_per_snr(2,i) = R;
        prf_per_snr(3,i) = F;
    end
end

% Function #2: Compute Temporal Concurrence between the Simulated and Detected Bursts
function [tc_per_snr] = extract_temporal_concurrence(statistics,listCycle)
    statistics = [statistics{:}];
    tpIdx = strcmp(statistics(3,:),'tp');
    statistics = statistics(:,tpIdx);
    statistics = sortrows(statistics(4:end,:)',1);
    tc_per_snr = zeros(1,length(listCycle));
    for i = 1:length(listCycle)
        c = listCycle(i);
        idx = cellfun(@(x) x==c, statistics(:,1));
        if ~any(idx)
            tc = NaN;
        else
            tc = mean(cell2mat(statistics(idx,end)));
        end
        tc_per_snr(1,i) = tc;
    end
end