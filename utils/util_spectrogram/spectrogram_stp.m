function [Spec_f, Spec_t, SpecMat, fft_interval] = spectrogram_stp(time, X, Fs, freq_range, dt, nOscCycle, verbose, interp_opt)
    %% Function 'spectrogram_stp'
    % DESCRIPTION
    % This function computes and plots a power spectrogram of an input data
    % smoothed by Hanning windows with frequency-dependent lengths using the FFT 
    % algorithm. The transform implemented here can be thought as an equivalent of
    % Constant-Q transform and is named as single-tapered STFT in the paper. The 
    % current version only supports a single-trial power spectrogram (i.e. works 
    % only for an input data of a single vector array).

    % USAGE
    % Full Input : spectrogram_stp(time, X, Fs, freq_range, dt, nOscCycle, verbose, interp_opt)
    % Example    : spectrogram_stp(time, data, 1024, unique([3:5 5:5:30 30:10:60]), 1/Fs, 6, true, true)

    % INPUT
    %    Variable       Data Type                    Description
    % 1. time           [1 x N vector]             : time array
    %                                                1) []
    %                                                2) vector array
    % 2. X              [1 x N vector]             : input data
    % 3. Fs             [integer Z > 0]            : sampling frequency rate
    % 4. freq_range     [1 x N vector]             : the array of frequencies of interest
    % 6. dt             [rational Q > 0]           : time step by which the window slides (seconds)
    %                                                default) 0.1
    % 7. nOscCycle      [integer Z > 0]            : number of oscillatory cycles to define the window size
    %                                                default) 6
    % 8. verbose        [logical]                  : option to plot the spectrogram and display summary
    %                                              : default) true
    % 9. interp_opt     [logical]                  : option to interpolate 2-D gridded spectrogram
    %                                              : default) true

    % OUTPUT
    %    Variable       Data Type                    Description
    % 1. Spec_f         [1 x N vector]             : vector of frequency bins for the spectrogram
    % 2. Spec_t         [1 x N vector]             : vector of times for the spectrogram
    % 3. SpecMat        [N x N matrix]             : matrix of calculated spectral power
    % 4. fft_interval   [N x 2 matrix]             : matrix of frequency ranges centered around each element of Spec_f

    % Written by SungJun Cho, March 24, 2021
    % Last modified on February 07, 2023
    %% Input Parameters
    % [1] Set Default Values
    if nargin < 8
        interp_opt = true;
    end
    if nargin < 7
        verbose = true;
    end
    if nargin < 6
        nOscCycle = 6;
    end
    if nargin < 5
        dt = 0.1;
    end
    if nargin < 4
        error('INPUT ERROR: Need to define a specific frequency range.');
    end
    if freq_range(1) == 0
        warning('Frequency-dependent window cannot be defined for DC frequency. The lowest frequency will be reset to 1 Hz.');
        freq_range(1) = 1;
    end
    if freq_range(end) > Fs/2
        warning('The highest frequency should be less than the Nyquist frequency and will be reset to the corresponding Nyquist frequency.');
        freq_range(end) = Fs/2;
    end
    % [2] Create Time Array
    if isempty(time)
        time = linspace(0, length(X)/Fs, length(X));
    end
    if length(time) ~= length(X)
        error('INPUT ERROR: The sizes of a time and data array should match.');
    end
    % [3] Check Input Type
    if ~isvector(time) || ~isvector(X)
        error('Your input data and/or time array should be be a vector array.');
    else
        if size(time,1) ~= 1
            time = time';
        end
        if size(X,1) ~= 1
            X = X';
        end
    end
    %% Spectrogram Parameters
    window_sizes = zeros(1,length(freq_range));
    dfs = zeros(1,length(freq_range));
    % [1] Define Window Sizes
    for w = 1:length(window_sizes)
        window_sizes(w) = round(nOscCycle*(Fs/freq_range(w)));
        if mod(window_sizes(w),2) ~= 0
            window_sizes(w) = window_sizes(w)+1;
            % ensure windows to be of even length to prevent non-integer indexing error
        end
        dfs(w) = Fs/window_sizes(w);
    end
    % [2] Define Time Arrays
    moving_step = floor(Fs*dt); % time step in sampling points for sliding windows
    time_new = linspace(time(1),time(end),length(time)/moving_step);
    %% Compute Spectrogram
    % [1] Preallocate Cell and Double Arrays
    SpecCell = cell(1,length(window_sizes));
    TimeCell = cell(1,length(window_sizes));
    freq_vec = zeros(1,length(window_sizes));
    IndCell = cell(1,length(window_sizes));
    fft_interval = zeros(length(window_sizes),2);
    
    for w = 1:length(window_sizes)
        % [2] Get sliding window indices
        idx_start = 1+floor(window_sizes(w)/2);
        idx_end = length(time) - floor(window_sizes(w)/2);
        idx_midpts = idx_start:moving_step:idx_end;

        idx_collection = zeros(length(idx_midpts),2);
        for idx = 1:length(idx_midpts)
            idx_collection(idx,1) = idx_midpts(idx) - window_sizes(w)*0.5;
            idx_collection(idx,2) = idx_midpts(idx) + window_sizes(w)*0.5 - 1;
        end

        idx_diff = idx_collection(:,2) - idx_collection(:,1);
        if idx_diff ~= (window_sizes(w)-1)
            error('ERROR: Index computation is incorrect.');
        end

        % [3] Compute Scaled Frequency Vector
        L = window_sizes(w);
        f_org = (Fs*(0:L/2)/L);
        f = round(f_org);
        if f ~= round(0:dfs(w):Fs/2)
            error('EEROR: Frequency vector computation has encountered an error.');
        end
        subtracted = abs(f_org - freq_range(w));
        select_freq = find(subtracted == min(subtracted));
        freq_vec(1,w) = f_org(select_freq);

        % [4] Pad Zeros to Equalize Different Time Vector Lengths
        time_temp = time(idx_midpts);
        t_vec = time_new;
        cut_start = find(time_new < time_temp(1), 1, 'last');
        cut_end = find(time_new > time_temp(end), 1, 'first');
        if ~isempty(cut_start)
            t_vec(1:cut_start) = 0;
        else
            cut_start = 0;
        end
        if ~isempty(cut_end)
            t_vec(cut_end:end) = 0;
        else
            cut_end = length(time_new)+1;
        end
        t_start = cut_start+1;
        t_end = cut_end-1;

        % [5]Temporal Smoothing with Frequency-Dependent Window Length
        Spec_temp = zeros(length(time_new),1);
        hanning = hann(window_sizes(w))';

        nfft = window_sizes(w);

        for tIdx = 1:length(t_start:t_end)
            x_hann = hanning .* X(idx_collection(tIdx,1):idx_collection(tIdx,2)); 
            x_fft = fft(x_hann,nfft)/Fs;
            abs_val = conj(x_fft).*x_fft;
            Spec_temp(tIdx+t_start-1,1) = abs_val(select_freq);
        end

        % [6] Store Results
        SpecCell{w} = Spec_temp;
        TimeCell{w} = t_vec;
        IndCell{w} = idx_collection;
        fft_interval(w,:) = [f_org(select_freq)-(dfs(w)/2) f_org(select_freq)+(dfs(w)/2)];
    end
    %% Output Results
    Spec_t = time_new;
    Spec_f = freq_vec;
    SpecMat = [SpecCell{:}];
    %% Plot Figure
    if ~verbose
        if length(Spec_f) ~= length(unique(Spec_f))
            error('Redundant frequency value occurred during the estimation of scaled frequency vector.');
        end
    else
        if length(Spec_f) == length(unique(Spec_f))
            figure();
            [~,obj] = contourf(Spec_t,Spec_f,SpecMat');
            obj.LineStyle = 'none';
            set(gca,'YTick',round(Spec_f,2));
            xlabel('Time (s)');
            ylabel('Frequency (Hz)');
            title('S-STFT Spectrogram (f-Dependent Window)')
            c = colorbar; colormap('jet');
            ylabel(c, 'Power (mV^2)','FontSize',12);
        else
            error('Redundant frequency value occurred during the estimation of scaled frequency vector.');
        end
        if interp_opt
            interp_method = 'spline'; % or alternatively: 'linear' or 'nearest'
            Spec_f_ip = Spec_f(1):Spec_f(end);
            Spec_t_ip = Spec_t;
            [Xq, Yq] = meshgrid(Spec_f_ip,Spec_t_ip);
            SpecInterp = interp2(Spec_f,Spec_t,SpecMat,Xq,Yq,interp_method);
            % visualized center frequency may deviate from the actual center
            % frequency after interpolation
            figure();
            imagesc(Spec_t_ip,Spec_f_ip,imgaussfilt(SpecInterp',2));
            % spectral power may decrease when Gaussian filtered
            axis xy;
            xlabel('Time (s)');
            ylabel('Frequency (Hz)');
            title('S-STFT Spectrogram (Frequency Interpolated)');
            c = colorbar; colormap('jet');
            ylabel(c, 'Power (mV^2)','FontSize',12);
        end
    end
    %% Display Spectrogram Properties
    if verbose
        display_summary(dfs, dt, window_sizes, Fs);
    end
end

%% Appendix: In-Script Functions
% Function #1: Displays a summary of S-STFT sepctrogram parameters and
% resoultions
function display_summary(dfs, dt, window_sizes, Fs)
    disp('*****************SUMMARY*****************');
    disp('Single-Tapered STFT Spectrogram Parameters:');
    disp(['    Frequency Bin / FFT binwdith (Hz): ' num2str(dfs)]);
    disp(['    Window Step / Time Bin (s): ' num2str(dt)]);
    disp(['    Window Length (s): ' num2str(window_sizes/Fs)]);
    disp('*****************************************');
end