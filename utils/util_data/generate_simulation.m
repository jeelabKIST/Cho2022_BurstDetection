function [tvec, simulated_signal, Fs, data, btime, burst, bfreq, cyc] = generate_simulation(method, T, f, f_spans, nBurst, SNR_wn, rng_pks, rnd_seed, verbose)
    %% Function: 'generate_simulation'
    % DESCRIPTION
    % Generates a simulated signal containing aperiodic 1/f background noise,
    % white Gaussian noise, and bursts of three separate frequencies with predefined 
    % frequency spans. Bursts are rendered in random inter-burst intervals.

    % USAGE
    % Full Input : generate_simulation_rnd(sim_method, T, f, f_spans, n, SNR, range_pks, nSeed, sim_verbose)
    % Example    : generate_simulation_rnd('tukey', 60, 25, [1,2,4], 15, 5, [3,12], 97, true)

    % INPUT
    %    Variable       Data Type                    Description
    % 1. method         [string]                   : simulation method
    %                                                1) 'tukey' (default): uniform waveform
    %                                                2) 'gauss'          : wax-and-wane (bell-curved) waveform
    % 2. T              [number N]                 : total duration of the signal
    % 3. f              [number N]                 : lowest frequency for building the burst components
    % 4. f_spans        [array of number N]        : frequency spans (i.e.,bandwidths) for each burst frequency
    % 4. nBurst         [number N]                 : number of the bursts for each frequencies
    % 5. SNR_wn         [number N]                 : signal-to-noise ratio (dB) of the signal
    % 6. rng_pks        [array of number N]        : number of peaks (i.e., cycles) of the bursts
    % 7. rnd_seed       [number N]                 : seed number for the random number generator
    % 8. verbose        [logical]                  : indicator for generating an output plot
    %                                                1) false or 0 (default)
    %                                                2) true or 1

    % OUTPUT
    % INPUT
    %    Variable          Data Type                 Description
    % 1. tvec              [N x 1 vector]          : time vector array
    % 2. simulated_signal  [N x 1 vector]          : simulation of the signal
    % 3. Fs                [number N]              : sampling rate of the signal
    % 4. data              [N x N double]          : data consisting of background noise and two types of bursts
    % 5. btime             [N x N cell]            : set of time arrays for the bursts of different frequencies
    % 6. burst             [N x N cell]            : set of amplitude arrays for the bursts of different frequencies
    % 7. bfreq             [N x N double]          : set of frequencies characterizing each burst
    % 8. cyc               [N x N double]          : length (cycles) for the bursts of different frequencies

    % REFERENCES
    % [1] Quinn A.J., van Ede F., Brookes M.J., Heideman S.G., Nowak M., Seedat Z.A., Vidaurre D., Zich C., Nobre A.C., Woolrich M.W. (2019).
    %     Unpacking Transient Event Dynamics in  Electrophysiological Power Spectra. Brain Topography. https://doi.org/10.1007/s10548-019-00745-5.
    %     MIT License, Copyright (c) 2019 OHBA Analysis Group
    %     Partly adapted from "Quinn2019_BurstHMM/hmm_util_get_simulation.m"

    % Written by SungJun Cho, January 8, 2023
    % Last modified on February 04, 2023
    %% Configure Input Parameters
    % [1] Decide plotting availability
    if nargin < 9
        verbose = false;
    end
    % [2] Assign random seed number
    if nargin < 8
        rnd_seed = 97;
    end
    % [3] Check the method
    if nargin < 1 || isempty(method)
        method = 'tukey';
    end
    if isempty(validatestring(method,{'tukey','gauss'}))
        error('InputError: Simulation method must be either tukey or gauss.');
    end
    %% Generate Time-Series Simulation
    % [1] Control the random number generator with preset seed number
    rng(rnd_seed);
    % [2] Set simulation parameters
    Fs = 512; % sampling rate
    tvec = linspace(0,T,T*Fs)'; % time vector
    % [3] Create a noisy time course
    % The time series will approximately follow 1/f power scale by direct pole-placement
    a = -poly(.92); % the filter root; the denominator polynomial of a digital filter
    % The roots given by "poly" should be between 0 <= x <= 1 where 0 is white noise and
    % 1 is a spectrum with extreme slope
    pink_noise = 0.2*randn(size(tvec));
    x_noise = filter(1,a,pink_noise);
    % [4] Create synthetic oscillations
    chi = 0.15; % scaling constant
    steps = 4:5; % constants used to define higher frequencies
    shifts = -max(f_spans):max(f_spans); % frequency span
    nOsc = length(steps) + 1;
    osc = cell(1,length(shifts));
    for s = 1:length(shifts)
        [signals, freqs] = generate_multifreq_osc(f+shifts(s),tvec,Fs,chi,steps);
        osc{s} = signals;
        if shifts(s) == 0
            freq_list = freqs;
        end
    end
    % [5] Define burst occurrences and durations based on the center
    %     frequencies
    replacement = true;
    [start,duration] = define_on_off(Fs,freq_list,tvec,nBurst,rng_pks,replacement);
    % [6] Create burst-embedded time-series
    x = zeros(length(tvec), nOsc);
    [btime,burst] = deal(cell(nBurst,nOsc));
    [bfreq,cyc] = deal(zeros(nBurst,nOsc));
    switch method
        case 'tukey'
            r = 0.25; % cosine fraction; default = 0.5
                      % r = 0 -> a rectangular window
                      % r = 1 -> equivalent to a Hanning window
        case 'gauss'
            alpha = 1.2; % width factor; default = 2.5
                         % selected such that when the duration is about 3 cycles, all the
                         % peaks show sufficient magnitude of amplitude
    end
    for n = 1:nOsc
        for t = 1:length(start{n})
            span = f_spans(n);
            shift_idx = 1:length(shifts);
            shift_idx = shift_idx(abs(shifts) == span);
            k = datasample(shift_idx(1):shift_idx(end),1,'Replace',replacement);
            temp_sig = osc{k}(start{n}(t):start{n}(t)+duration{n}(t),n);
            if strcmp(method,'tukey')
                temp_sig = temp_sig.*tukeywin(length(temp_sig),r);
            elseif strcmp(method,'gauss')
                temp_sig = temp_sig.*gausswin(length(temp_sig),alpha);
            end
            x(start{n}(t):start{n}(t)+duration{n}(t),n) = temp_sig;
            % Store bursts
            btime{t,n} = tvec(start{n}(t):start{n}(t)+duration{n}(t));
            burst{t,n} = temp_sig;
            bfreq(t,n) = freq_list(n)+shifts(k);
            cyc(t,n) = round(duration{n}(t)/(Fs/bfreq(t,n)));
        end
    end
    % [7] Add noise and burst time-series to generate output signal
    simulated_signal = x_noise + sum(x,2); % superpose each signal component
    % [8] Add Gaussian white noise
    gwn = randn(size(simulated_signal)) * ((1/(f^chi))/db2mag(SNR_wn));
    % use the scaling factor for the lowest frequency as a variance
    simulated_signal = simulated_signal + gwn;
    % References: (1) https://www.mathworks.com/help/matlab/ref/randn.html
    %             (2) https://www.mathworks.com/help/matlab/math/random-numbers-with-specific-mean-and-variance.html
    %             (3) https://www.mathworks.com/matlabcentral/answers/472495-generating-random-variable-from-certain-standard-deviation-and-mean
    %             (4) https://en.wikipedia.org/wiki/Additive_white_Gaussian_noise
    % [9] Concatenate different parts of simulation
    data = cat(2,simulated_signal,x,x_noise,gwn);
    %% Plot Simulations
    if verbose
        % [1] Set visualization parameters
        nColumns = size(data,2);
        text_titles = cell(1,size(x,2));
        for i = 1:size(x,2)
            text_titles{i} = ['Burst #' num2str(i) ' (f = ' num2str(freq_list(i)) ' +/-' num2str(f_spans(i)) ' Hz)'];
        end
        text_titles = [{'Complete Data'}, text_titles, 'Aperiodic Noise', 'Gaussian Noise'];
        clr_array = [0,0,0; 123,213,245; 28,167,236; 31,47,152; 162,20,47; 217,83,25]/255;
        % [2] Visualize a summary of simulated signals
        figure();
        for i = 1:nColumns
            subplot(nColumns,1,i);
            plot(tvec, data(:,i),'Color',clr_array(i,:));
            xlabel('Time (s)');
            ylabel('Amplitude (a.u.)');
            title(text_titles{i},'FontWeight','bold');
            set(gca,'FontSize',12);
            set(gcf,'Color','w','Position',[83,37,1344,817]);
        end
    end
end

function [osc, freqs] = generate_multifreq_osc(f, t, Fs, chi, steps)
    %% Function: 'generate_multifreq_osc'
    % DESCRIPTION
    % Generates a set of synthetic sinusoidal oscillations, scaled by the
    % 1/f^chi power law.

    % USAGE
    % Full Input : generate_multifreq_osc(f, tvec, Fs, chi, steps)
    % Example    : generate_multifreq_osc(10, tvec, 512, 2, [2,3,4])

    % INPUT
    %    Variable       Data Type                       Description
    % 1. f              [number N]                    : lowest frequency of an oscillatory signal
    % 2. t              [array of integers Z >= 0]    : total duration of the signal
    % 3. Fs             [number N]                    : sampling frequency
    % 4. chi            [integer Z >= 0]              : scaling factor for scale-free power law
    % 5. steps          [array of number N]           : contants that define higher frequencies

    % OUTPUT
    % INPUT
    %    Variable          Data Type                 Description
    % 1. osc               [N x N vector]          : synthetic oscillations
    % 2. freqs             [1 x N vector]          : array of frequencies characterizing each oscillatory signal

    % Written by SungJun Cho, January 8, 2023
    % Last modified on January 8, 2023
    %% Generate Synthetic Sinusoidal Oscillations
    % [1] Compute oscillatory time-series
    nOsc = length(steps)+1;
    osc = zeros(length(t), nOsc);
    freqs = zeros(1,nOsc);
    for i = 1:length(steps)
        n = steps(i);
        f_prime = f+(Fs/(2^n));
        osc(:,i) = (1/(f_prime^chi))*sin(2*pi*f_prime*t);
        freqs(i) = f_prime;
    end
    osc(:,end) = (1/(f^chi)) * sin(2*pi*f*t);
    freqs(end) = f;
    % [2] Order outputs from the lowest to highest frequencies
    osc = fliplr(osc);
    freqs = sort(freqs);
end

function [start, duration] = define_on_off(Fs, freqs, tvec, nBurst, rng_pks, replacement)
    %% Function: 'define_on_off'
    % DESCRIPTION
    % Defines the starting time and duration of bursts of different frequencies

    % USAGE
    % Full Input : define_on_off(Fs, freqs, tvec, nBurst, rng_pks, replacement)
    % Example    : define_on_off(512, freq_list, tvec, 15, [3,12], true)

    % INPUT
    %    Variable       Data Type                       Description
    % 1. Fs             [number N]                    : sampling frequency
    % 2. freqs          [1 x N array]                 : list of frequencies used to generate bursts
    % 3. tvec           [N x 1 array]                 : time array (in seconds)
    % 4. nBurst         [number N]                    : number of the bursts for each frequencies
    % 5. rng_pks        [array of number N]           : number of peaks (i.e., cycles) of the bursts
    % 6. replacement    [boolean]                     : indicator for sampling with replacement

    % OUTPUT
    % INPUT
    %    Variable          Data Type                 Description
    % 1. start             [1 x N cell]            : cell array containing a set of burst onset times for each frequency
    % 2. duration          [1 x N cell]            : cell array containing a set of burst durations for each frequency

    % Written by SungJun Cho, January 8, 2023
    % Last modified on January 8, 2023
    %% Compute burst onsets and offsets
    % [1] Set parameters and preallocate outputs
    Nyq = Fs/2; % Nyquist frequency
    nFreqs = length(freqs);
    [start, duration] = deal(cell(1, nFreqs));
    % [2] Calucate starting time and duration of each burst
    for n = 1:nFreqs
        start{n} = sort(randsample(round(Nyq)+1:length(tvec)-(round(Nyq)+1)-(rng_pks(2)*(Fs/freqs(n))),nBurst)); % put ':Fs:' to enforce random burst intervals with at least 1s gap
        available_duration = round((rng_pks(1):rng_pks(2)).*(Fs/freqs(n)));
        duration{n} = randsample(available_duration,nBurst,replacement);
        if length(start{n}) ~= nBurst
            error('ComputationError: The number of generated bursts does not match the preset number.');
        end
        for t = 1:length(start{n})
            % Consider the case in which the predefined burst duration exceeds the signal length
            if t == length(start{n})
                if start{n}(t)+duration{n}(t) > length(tvec)
                    duration{n}(t) = length(tvec) - start{n}(t);
                end
                break
            end
            % Consider the case in which the next burst overwrites the previous burst
            % Reinforce the distance of 2 cycles as a gap between subsequent bursts in this case
            if start{n}(t)+duration{n}(t) > start{n}(t+1)
                start{n}(t+1) = start{n}(t+1) + (start{n}(t)+duration{n}(t)-start{n}(t+1))+round(2*Fs/freqs(n));
            end
        end
    end
    % Note: Nyquist frequency implemented here to avoid generating bursts at the area (i.e., at the start and the end) uncovered by 
    % the spectrogram computation (under the assumption that the window is less than 1s).
end
