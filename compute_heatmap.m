%% Configure Library Path
util_path = genpath('/Users/scho/Neuroscience_KIST/Cho2022_BurstDetection/utils');
data_path = genpath('/Users/scho/Neuroscience_KIST/Cho2022_BurstDetection/data/simulation_data');
addpath(util_path);
addpath(data_path);
%% Set System Parameters
filepath_seed = fullfile(pwd,'data','simulation_data','randseed.mat');
% To reproduce a previous result, save and call a set of random seeds that 
% was used to simulate the signals. The .mat file should be located in the
% directory corresponding to `data_path` above. 
filename_result = 'HM.mat'; % file name for the results to be saved
%% Theoretical Simulation of the Signals
fprintf('Starting Signal Simulation ... \n');
% [1] Set Fixed Parameters
method = 'tukey';
T = 300;               % unit: s
f = 9;                 % unit: Hz
f_spans = [1,2,4];
nBurst = 25;
rng_pks = [3 12];
sim_verbose = false;
listNoise = -10:2:10; % list of noise levels (dB)
listCycle = rng_pks(1):rng_pks(end); % list of burst cycle lengths
% [2] Set Varying Parameters
nSimulationPerNoise = 1000; % number of simulations per each SNR
if isfile(filepath_seed)
    setSeed = load(filepath_seed).setSeed;
    fprintf('*** Note: Random seeds imported for reproducibility. \n');
else
    setSeed = zeros(length(listNoise),nSimulationPerNoise);
    for n = 1:length(listNoise)
        setSeed(n,:) = randsample(1:2^15,nSimulationPerNoise);
    end
    save(filepath_seed,'setSeed');
end
% [3] Preallocate Dataset
sz_dataset = cell(length(listNoise),nSimulationPerNoise);
[ds_signal,ds_burst,ds_btime,ds_bfreq,ds_length] = deal(sz_dataset);
% [4] Get Simulation Data
tic;
for n = 1:length(listNoise)
    % Set SNR (dB) Level
    SNR_wn = listNoise(n);
    for s = 1:size(setSeed,2)
        % Set Random Seed
        rnd_seed = setSeed(n,s);
        [~,simulated_signal,Fs,~,btime_true,burst_true,bfreq_true,bcycl_true] ...
            = generate_simulation(method,T,f,f_spans,nBurst,SNR_wn,rng_pks,rnd_seed,sim_verbose);
        ds_signal{n,s} = simulated_signal;
        ds_burst{n,s} = burst_true;
        ds_btime{n,s} = btime_true;
        ds_bfreq{n,s} = bfreq_true;
        ds_length{n,s} = bcycl_true;
    end
    fprintf(['Progress: running ... [' num2str(n*s) '/' num2str(length(listNoise)*nSimulationPerNoise) '] \n']);
end
tvec = linspace(0,T,T*Fs)'; % time vector
Ps = 1/Fs; % sampling period
cnt_f = get_multifreq(f,Fs); % center frequencies of the frequency bands used in the simulation
nBand = length(cnt_f);
elapsed_time = toc;
clc; fprintf(['STEP1: Simulation complete! (computational time: ' num2str(elapsed_time) 's) \n']);
%% Set Parameters
fprintf('Starting parameter setting ... \n'); tic;
% [1] Time Domain: Bandpass Filtering
lfs = cnt_f - f_spans; hfs = cnt_f + f_spans;
% [2] Time-Frequency Domain: Power Spectrograms
% (A) Spectrogram Parameters
freq_range = unique([5:2:15 15:5:30 30:10:60]);          % frequencies of interest
dt = 0.01;                                               % time step for sliding windows
nOscCycle_stp = 6; nOscCycle_mtp = 6;                    % frequency-dependent window size
nw = 2; ntp = 3;                                         % DPSS slepian sequences
wname = 'amor'; scale_opt = 'default';                   % type of and scaling option for the wavelet
verbose = false; interp_opt = false;                     % visualization options
plot_psd = false; check_opt = true; plot_tp = false;     % visualization and optimization options
% (B) Amplitude Thresholds
thrA_factor = 2.5;                % for bandpass filtering method
prcthrA = 95;                     % for envelope-based method
specthrA_factor = 1.8;            % for spectral analysis methods
% (C) Window Time Step Used
dt_stp = dt;
dt_mtp = dt;
dt_cwt = Ps;
% For cwt.m, translation is fixed to 1.
% References: https://www.mathworks.com/matlabcentral/answers/283547-continuous-wavelet-transform-implementation-using-morlet?s_tid=srchtitle
%             https://www.mathworks.com/help/wavelet/gs/interpreting-continuous-wavelet-coefficients.html
%             https://www.mathworks.com/help/wavelet/gs/continuous-and-discrete-wavelet-transforms.html
% (D) Frequency and Time Vectors for Each Method
% Note: These vectors would be equivalent for every spectrogram, so it is
% unnecessary to output them every time.
input_signal = ds_signal{1,1};
[Spec_f_stp,Spec_t_stp,~,fft_interval_stp] = spectrogram_stp(tvec,input_signal,Fs,freq_range,dt,nOscCycle_stp,verbose,interp_opt);
[Spec_f_mtp,Spec_t_mtp,~,fft_interval_mtp] = spectrogram_mtp(tvec,input_signal,Fs,nw,ntp,freq_range,dt,nOscCycle_mtp,interp_opt,plot_psd,check_opt,plot_tp,verbose);
[wav_freq,wav_time,~,coi] = spectrogram_cwt(tvec,input_signal,Fs,wname,scale_opt,verbose);
% (E) Frequency of Interest
parse_start = find(wav_freq<freq_range(end),1,'first');
parse_end  = find(wav_freq>freq_range(1),1,'last');
wav_freq = wav_freq(parse_start:parse_end);
[minIdx_stp, minIdx_mtp, minIdx_cwt] = deal(zeros(1,nBand));
for i = 1:nBand
    fft_dist_stp = abs(fft_interval_stp(:,2)-cnt_f(i)) + abs(fft_interval_stp(:,1)-cnt_f(i));
    fft_dist_mtp = abs(fft_interval_mtp(:,2)-cnt_f(i)) + abs(fft_interval_mtp(:,1)-cnt_f(i));
    [~,minIdx_stp(i)] = min(fft_dist_stp);
    [~,minIdx_mtp(i)] = min(fft_dist_mtp);
    [~,minIdx_cwt(i)] = min(abs(wav_freq-cnt_f(i)));
end
fprintf(['Frequencies of interest selected from STP are ' num2str(Spec_f_stp(minIdx_stp)) '. \n']);
fprintf(['Frequencies of interest selected from MTP are ' num2str(Spec_f_mtp(minIdx_mtp)) '. \n']);
fprintf(['Frequencies of interest selected from CWT are ' num2str(wav_freq(minIdx_cwt)') '. \n']);
% [3] Report Computational Time
elapsed_time = toc;
clc; 
fprintf(['STEP2: Parameters set! (computational time: ' num2str(elapsed_time) 's) \n']);
%% Construct Heatmaps
fprintf('Starting heatmap construction ... \n'); tic;
% [1] Preallocate Heatmaps
[P_bp,P_ev,P_stp,P_mtp,P_cwt] = deal(zeros(nBand,length(listNoise),length(listCycle))); % precision
[R_bp,R_ev,R_stp,R_mtp,R_cwt] = deal(zeros(nBand,length(listNoise),length(listCycle))); % recall (sensitivity)
[F_bp,F_ev,F_stp,F_mtp,F_cwt] = deal(zeros(nBand,length(listNoise),length(listCycle))); % F1-score
[T_bp,T_ev,T_stp,T_mtp,T_cwt] = deal(zeros(nBand,length(listNoise),length(listCycle))); % temporal concurrence
for rowIdx = 1:size(ds_signal,1)
    % [2] Preallocate Burst Statistics Outputs
    [statistics_bp, statistics_ev, statistics_stp, statistics_mtp, statistics_cwt] = deal(cell(nBand,size(ds_signal,2)));
    for colIdx = 1:size(ds_signal,2)
        %% Preprocessing of the Simulated Signals
        input_signal = ds_signal{rowIdx,colIdx};
        % [3] Time-Frequency Domain: Power Spectrograms
        % (A) Spectrogram: Single-Tapered Short Time Fourier Transform
        [~,~,Spec_stp,~] = spectrogram_stp(tvec,input_signal,Fs,freq_range,dt_stp,nOscCycle_stp,verbose,interp_opt);
        % (B) Spectrogram: Multitaper
        [~,~,Spec_mtp,~] = spectrogram_mtp(tvec,input_signal,Fs,nw,ntp,freq_range,dt_mtp,nOscCycle_mtp,interp_opt,plot_psd,check_opt,plot_tp,verbose);
        % (C) Spectrogram: Continuous Wavelet Transform
        [~,~,Spec_cwt,~] = spectrogram_cwt(tvec,input_signal,Fs,wname,scale_opt,verbose);
        % (D) Parse CWT Spectrogram for Memory
        Spec_cwt = Spec_cwt(parse_start:parse_end,:); % note that dimension is transposed for cwt
        for bandIdx = 1:nBand
            % [4] Time Domain: Bandpass Filtering
            lf = lfs(bandIdx); hf = hfs(bandIdx);
            filtered_signal = eegfilt(input_signal',Fs,lf,hf);
            %% Burst Detection & Statistics
            % [5] Detect Bursts
            % (A) Get True Simulated Bursts
            btime_true = ds_btime{rowIdx,colIdx}(:,bandIdx)';
            burst_true = ds_burst{rowIdx,colIdx}(:,bandIdx)';
            bfreq_true = ds_bfreq{rowIdx,colIdx}(:,bandIdx)';
            bcycl_true = ds_length{rowIdx,colIdx}(:,bandIdx)';
            % (B) Detect Bursts: Time Domain
            [btime_bp,burst_bp] = detect_burst_timeseries(Fs,tvec,filtered_signal,lf,hf,thrA_factor);
            [btime_ev,burst_ev,benvelope] = detect_burst_ampenv(Fs,tvec,filtered_signal,lf,hf,prcthrA);
            % (C) Detect Bursts: Time-Frequency Domain
            [btime_stp,burst_stp] = detect_burst_spectrogram(Spec_f_stp,Spec_t_stp,Spec_stp,dt_stp,specthrA_factor,verbose);
            [btime_mtp,burst_mtp] = detect_burst_spectrogram(Spec_f_mtp,Spec_t_mtp,Spec_mtp,dt_mtp,specthrA_factor,verbose);
            [btime_cwt,burst_cwt] = detect_burst_spectrogram(wav_freq,wav_time,Spec_cwt',dt_cwt,specthrA_factor,verbose);
            [btime_stp,burst_stp] = extract_single_freqband(btime_stp,burst_stp,minIdx_stp(bandIdx));
            [btime_mtp,burst_mtp] = extract_single_freqband(btime_mtp,burst_mtp,minIdx_mtp(bandIdx));
            [btime_cwt,burst_cwt] = extract_single_freqband(btime_cwt,burst_cwt,minIdx_cwt(bandIdx));
            % (D) Compute Burst Statistics
            [bstat_bp] = get_burst_statistics(btime_true,bcycl_true,btime_bp,Ps,cnt_f(bandIdx));
            [bstat_ev] = get_burst_statistics(btime_true,bcycl_true,btime_ev,Ps,cnt_f(bandIdx));
            [bstat_stp] = get_burst_statistics(btime_true,bcycl_true,btime_stp,dt_stp,cnt_f(bandIdx));
            [bstat_mtp] = get_burst_statistics(btime_true,bcycl_true,btime_mtp,dt_mtp,cnt_f(bandIdx));
            [bstat_cwt] = get_burst_statistics(btime_true,bcycl_true,btime_cwt,dt_cwt,cnt_f(bandIdx));
            % (E) Store Burst Statistics
            statistics_bp{bandIdx,colIdx} = bstat_bp;
            statistics_ev{bandIdx,colIdx} = bstat_ev;
            statistics_stp{bandIdx,colIdx} = bstat_stp;
            statistics_mtp{bandIdx,colIdx} = bstat_mtp;
            statistics_cwt{bandIdx,colIdx} = bstat_cwt;
            % [6] Report Progress
            if mod(colIdx,round(size(ds_signal,2)/10)) == 0
                fprintf(['Progress: running ... [SNR: ' num2str(listNoise(rowIdx)) ' dB and trial #: ' num2str(colIdx) ' trials] \n']);
            end
        end
    end
    %% Precision, Recall, F1-Score & Temporal Concurrence Heatmap
    for bandIdx = 1:nBand
        % [7] Compute Precision, Recall, and F1-Score for Each SNR Value
        %     (1000 trials per one SNR value)
        prf_per_snr_bp = extract_precision_recall(statistics_bp(bandIdx,:),listCycle);
        prf_per_snr_ev = extract_precision_recall(statistics_ev(bandIdx,:),listCycle);
        prf_per_snr_stp = extract_precision_recall(statistics_stp(bandIdx,:),listCycle);
        prf_per_snr_mtp = extract_precision_recall(statistics_mtp(bandIdx,:),listCycle);
        prf_per_snr_cwt = extract_precision_recall(statistics_cwt(bandIdx,:),listCycle);
        % [8] Construct Heatmap
        P_bp(bandIdx,rowIdx,:) = prf_per_snr_bp(1,:);   P_ev(bandIdx,rowIdx,:) = prf_per_snr_ev(1,:);   P_stp(bandIdx,rowIdx,:) = prf_per_snr_stp(1,:);
        R_bp(bandIdx,rowIdx,:) = prf_per_snr_bp(2,:);   R_ev(bandIdx,rowIdx,:) = prf_per_snr_ev(2,:);   R_stp(bandIdx,rowIdx,:) = prf_per_snr_stp(2,:);
        F_bp(bandIdx,rowIdx,:) = prf_per_snr_bp(3,:);   F_ev(bandIdx,rowIdx,:) = prf_per_snr_ev(3,:);   F_stp(bandIdx,rowIdx,:) = prf_per_snr_stp(3,:);
        
        P_mtp(bandIdx,rowIdx,:) = prf_per_snr_mtp(1,:);   P_cwt(bandIdx,rowIdx,:) = prf_per_snr_cwt(1,:);
        R_mtp(bandIdx,rowIdx,:) = prf_per_snr_mtp(2,:);   R_cwt(bandIdx,rowIdx,:) = prf_per_snr_cwt(2,:);
        F_mtp(bandIdx,rowIdx,:) = prf_per_snr_mtp(3,:);   F_cwt(bandIdx,rowIdx,:) = prf_per_snr_cwt(3,:);
        % [9] Compute Temporal Concurrence for Each SNR Value
        tc_per_snr_bp = extract_temporal_concurrence(statistics_bp(bandIdx,:),listCycle);
        tc_per_snr_ev = extract_temporal_concurrence(statistics_ev(bandIdx,:),listCycle);
        tc_per_snr_stp = extract_temporal_concurrence(statistics_stp(bandIdx,:),listCycle);
        tc_per_snr_mtp = extract_temporal_concurrence(statistics_mtp(bandIdx,:),listCycle);
        tc_per_snr_cwt = extract_temporal_concurrence(statistics_cwt(bandIdx,:),listCycle);
        % [10] Construct Heatmap
        T_bp(bandIdx,rowIdx,:) = tc_per_snr_bp;     T_ev(bandIdx,rowIdx,:) = tc_per_snr_ev;     T_stp(bandIdx,rowIdx,:) = tc_per_snr_stp;
        T_mtp(bandIdx,rowIdx,:) = tc_per_snr_mtp;   T_cwt(bandIdx,rowIdx,:) = tc_per_snr_cwt;
    end
end
% [11] Report Computing Time
elapsed_hmp = toc;
clc; fprintf(['STEP3: Heatmaps complete! (computational time: ' num2str(elapsed_hmp) 's) \n']);
%% Save Results
P = {P_bp,P_ev,P_stp,P_mtp,P_cwt};
R = {R_bp,R_ev,R_stp,R_mtp,R_cwt};
F = {F_bp,F_ev,F_stp,F_mtp,F_cwt};
TC = {T_bp,T_ev,T_stp,T_mtp,T_cwt};
band_names = {'alpha','beta','gamma'};
for bandIdx = 1:nBand
    save_heatmap(P,R,F,TC,setSeed,band_names{bandIdx},bandIdx,filename_result);
end
%% Appendix: In-Script Funcitons
% Function #1: Save Band-Specific Heatmaps
function save_heatmap(P,R,F,TC,rndseed,name,bandIdx,filename_result)
    % Set Parameters
    filename_result = insertAfter(filename_result,'HM',['_' name]);
    filepath_result = fullfile(pwd,'data','simulation_data',filename_result); % path where results will be saved
    % Get Results
    h = helper;
    [P_bp,P_ev,P_stp,P_mtp,P_cwt] = h.unpack_data(P);
    [R_bp,R_ev,R_stp,R_mtp,R_cwt] = h.unpack_data(R);
    [F_bp,F_ev,F_stp,F_mtp,F_cwt] = h.unpack_data(F);
    [T_bp,T_ev,T_stp,T_mtp,T_cwt] = h.unpack_data(TC);
    % Select Band-Specific Heatmaps
    [P_bp,P_ev,P_stp,P_mtp,P_cwt] = extract_band_specific_heatmap(bandIdx,P_bp,P_ev,P_stp,P_mtp,P_cwt);
    [R_bp,R_ev,R_stp,R_mtp,R_cwt] = extract_band_specific_heatmap(bandIdx,R_bp,R_ev,R_stp,R_mtp,R_cwt);
    [F_bp,F_ev,F_stp,F_mtp,F_cwt] = extract_band_specific_heatmap(bandIdx,F_bp,F_ev,F_stp,F_mtp,F_cwt);
    [T_bp,T_ev,T_stp,T_mtp,T_cwt] = extract_band_specific_heatmap(bandIdx,T_bp,T_ev,T_stp,T_mtp,T_cwt);
    % Store Results in Structure
    HEATMAP = struct();
    HEATMAP.randseed = rndseed;
    HEATMAP.precision.bp = P_bp;    HEATMAP.recall.bp = R_bp;    HEATMAP.f1_score.bp = F_bp;    HEATMAP.concurrence.bp = T_bp;
    HEATMAP.precision.ev = P_ev;    HEATMAP.recall.ev = R_ev;    HEATMAP.f1_score.ev = F_ev;    HEATMAP.concurrence.ev = T_ev;
    HEATMAP.precision.stp = P_stp;  HEATMAP.recall.stp = R_stp;  HEATMAP.f1_score.stp = F_stp;  HEATMAP.concurrence.stp = T_stp;
    HEATMAP.precision.mtp = P_mtp;  HEATMAP.recall.mtp = R_mtp;  HEATMAP.f1_score.mtp = F_mtp;  HEATMAP.concurrence.mtp = T_mtp;
    HEATMAP.precision.cwt = P_cwt;  HEATMAP.recall.cwt = R_cwt;  HEATMAP.f1_score.cwt = F_cwt;  HEATMAP.concurrence.cwt = T_cwt;
    % Save Results
    if isfile(filepath_result)
        warning('Save attempt failed. A file with the specified name already exists.');
    else
        save(filepath_result,'HEATMAP');
    end
end

% Function #2: Select Heatmaps of Corresponding Frequency Band
function [varargout] = extract_band_specific_heatmap(bandIdx,varargin)
    varargout = cell(1,nargin-1);
    for n = 1:nargin-1
        varargout{n} = squeeze(varargin{n}(bandIdx,:,:));
    end
end
