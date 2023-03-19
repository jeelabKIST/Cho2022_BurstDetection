function [Spec_f, Spec_t, Spec] = spectrogram_stft(time, X, Fs, lf, hf, dt, window_size, verbose)
    %% Function 'spectrogram_stft'
    % DESCRIPTION
    % This function computes and plots a power spectrogram of an input
    % data based on the STFT algorithm, using a fixed Hanning window smoothing.
    % The current version only supports a single-trial power spectrogram 
    % (i.e. works only for an input data of a signle vector array).

    % USAGE
    % Full Input : spectrogram_stft(time, X, Fs, lf, hf, dt, window_size verbose)
    % Example    : spectrogram_stft(t_vec, data, 1024, 0, 100, 0.1, 1024, true)

    % INPUT
    %    Variable       Data Type                                    Description
    % 1. time           [1 x N vector]                             : time array
    %                                                                1) []
    %                                                                2) vector array
    % 2. X              [1 x N vector]                             : input data
    % 3. Fs             [integer Z > 0]                            : sampling frequency rate
    % 4. lf             [integer Z >= 0 & =< Nyq]                  : the lower threshold of frequency of interest
    % 5. hf             [integer Z >= 0 & =< Nyq]                  : the upper threshold of frequency of interest
    % 6. dt             [rational Q > 0]                           : time step by which the window slides (seconds)
    %                                                                default) 0.1
    % 7. window_size    [integer Z > 0]                            : the length of Hanning window
    %                                                                default) 2^(nextpow2(Fs))
    % 8. verbose        [logical]                                  : option to plot the spectrogram and display summary
    %                                                                default) true

    % OUTPUT
    %    Variable       Data Type                                    Description
    % 1. Spec_f         [1 x N vector]                             : vector of frequency bins for the spectrogram
    % 2. Spec_t         [1 x N vector]                             : vector of times for the spectrogram
    % 3. Spec           [length(Spec_t) x length(Spec_f) matrix]   : matrix of calculated spectral power

    % Written by SungJun Cho, March 24, 2021
    % Last modified on October 29, 2021
    %% Input Parameters
    % [1] Set Default Values
    if nargin < 8
        verbose = true;
    end
    if nargin < 7
        window_size = 2^(nextpow2(Fs));
    end
    if nargin < 6
        dt = 0.1;
    end
    if nargin < 4
        lf = 0;
        hf = Fs/2;
    end
    if (hf-lf) >= window_size
        error('Your frequency range should be smaller than the window size.');
    end
    % [2] Create Time Array
    if isempty(time)
        time = linspace(0, length(X)/Fs, length(X));
    end
    if length(time) ~= length(X)
        error('INPUT ERROR: The sizes of a time and data array should match.');
    end
    % [3] Optimize Computatioanl Speed
    if 2^nextpow2(window_size) ~= window_size
        warning('Consider having your window size in the power of 2 for efficient FFT computation.');
    end
    % [4] Adjust Window Size
    if mod(window_size,1) ~= 0
        window_size = round(window_size);
    end
    if mod(window_size,2) ~= 0
        window_size = window_size + 1;
        % ensure windows to be of even length to prevent non-integer indexing error
    end
    % [5] Check Input Type
    if ~isvector(time) || ~isvector(X)
        error('Your input data and/or time array should be in a vector format.');
    else
        if size(time,1) ~= 1
            time = time';
        end
        if size(X,1) ~= 1
            X = X';
        end
    end
    %% Calculate FFT-Hanning Spectrogram
    % [1] Get Sliding Window Indices
    moving_step = floor(Fs*dt); % moving step in sampling points
    idx_start = 1 + floor(window_size/2);
    idx_end = length(time) - floor(window_size/2);
    idx_midpts = idx_start:moving_step:idx_end;

    idx_collection = zeros(length(idx_midpts),2);
    for idx = 1:length(idx_midpts)
        idx_collection(idx,1) = idx_midpts(idx) - window_size*0.5;
        idx_collection(idx,2) = idx_midpts(idx) + window_size*0.5 - 1;
    end

    idx_diff = idx_collection(:,2) - idx_collection(:,1);
    if idx_diff ~= (window_size-1)
        error('ERROR: Index computation is incorrect.');
    end

    % [2] Define Frequency Resolution
    nfft = max(2^(nextpow2(window_size)),window_size);
    df = Fs/nfft;
    f = 0:df:Fs/2;

    select_freq = find(f >= lf & f <= hf);
    Spec_f = f(select_freq);
    Spec_t = time(idx_midpts);
    Spec = zeros(length(Spec_t),length(Spec_f));

    hanning = hann(window_size)';

    for tIdx = 1:length(Spec_t)
        x_hann = hanning .* X(idx_collection(tIdx,1):idx_collection(tIdx,2));
        x_fft = fft(x_hann,nfft)/Fs;
        % or alternatively: fft(x_hann)/length(x_hann); (only when window_size = Fs)
        abs_val = conj(x_fft).*x_fft;
        Spec(tIdx,:) = abs_val(select_freq);
        % or alternatively: Spec(tIdx,:) = abs(x_fft(select_freq));
    end

    % [3] Pad Zeros to Match Original Signal Size
    time_new = linspace(time(1),time(end),length(time)/moving_step);
    pad_len = (length(time_new) - length(Spec_t))/2;
    if floor(pad_len) == pad_len
        Spec = cat(1,zeros(pad_len,length(Spec_f)),Spec,zeros(pad_len,length(Spec_f)));
    else
        pad_len = floor(pad_len);
        Spec = cat(1,zeros(pad_len,length(Spec_f)),Spec,zeros(pad_len+1,length(Spec_f)));
    end
    Spec_t = time_new;
    %% Plot STFT Spectrogram
    if verbose
        figure();
        colormap('jet');
        imagesc(Spec_t, Spec_f, imgaussfilt(Spec',1));
        axis xy;
        xlabel('Time (s)');
        ylabel('Frequency (Hz)');
        c = colorbar;
        ylabel(c,'Power (mV^2)','FontSize',12);
        axis tight;
        title('Hanning Spectrogram with Fixed Window Length')
        % Display Spectrogram Properties
        display_summary(df, dt, window_size, Fs);
    end
end

%% Appendix: In-Script Functions
% Function #1: Displays a summary of STFT sepctrogram parameters and
% resoultions
function display_summary(df, dt, window_size, Fs)
    disp('*****************SUMMARY*****************');
    disp('Short-Time Fourier Transform Spectrogram Parameters:');
    disp(['    FFT Binwidth / Frequency Bin (Hz): ' num2str(df)]);
    disp(['    Window Step / Time Bin (s): ' num2str(dt)]);
    disp(['    Window Length (s): ' num2str(window_size/Fs)]);
    disp('*****************************************');
end