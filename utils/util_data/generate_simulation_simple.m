function [tvec, simulated_signal, Fs, data, btime_f1, burst_f1, cyc1] = generate_simulation_simple(method, T, f1, f2, nBurst1, nBurst2, SNR_wn, rng_pks, rnd_seed, verbose)
    %% Function: 'generate_simulation_simple'
    % DESCRIPTION
    % Generates a simulated signal containing aperiodic 1/f background noise and
    % bursts of two separate frequencies. Bursts are rendered in random inter-burst
    % intervals.

    % USAGE
    % Full Input : generate_simulation_rnd(sim_method, T, f_low, f_high, n1, n2, SNR, range_pks, nSeed, sim_verbose)
    % Example    : generate_simulation_rnd('tukey', 60, 25, 40, 30, 30, 5, [3,12], 97, true)

    % INPUT
    %    Variable       Data Type                    Description
    % 1. method         [string]                   : simulation method
    %                                                1) 'tukey' (default): uniform waveform
    %                                                2) 'gauss'          : wax-and-wane (bell-curved) waveform
    % 2. T              [number N]                 : total duration of the signal
    % 3. f1             [number N]                 : frequency of the bursts
    % 4. f2             [number N]                 : frequency of the bursts (should be different from 'f1')
    % 5. nBurst1        [number N]                 : number of the bursts of 'f1' frequency
    % 6. nBurst2        [number N]                 : number of the bursts of 'f2' frequency
    % 7. SNR_wn         [number N]                 : signal-to-noise ratio (dB) of the signal
    % 8. rng_pks        [array of number N]        : number of peaks (i.e., cycles) of the bursts
    % 9. rnd_seed       [number N]                 : seed number for the random number generator
    % 10. verbose       [logical]                  : indicator for generating an output plot
    %                                                1) false or 0 (default)
    %                                                2) true or 1

    % OUTPUT
    % INPUT
    %    Variable          Data Type                 Description
    % 1. tvec              [N x 1 vector]          : time vector array
    % 2. simulated_signal  [N x 1 vector]          : simulation of the signal
    % 3. Fs                [number N]              : sampling rate of the signal
    % 4. data              [N x 3 matrix]          : data consisting of background noise and two types of bursts
    % 5. btime_f1          [N x 1 cell vector]     : set of time arrays for the bursts of 'f1' frequency
    % 6. burst_f1          [N x 1 cell vector]     : set of amplitude arrays for the bursts of 'f1' frequency
    % 7. cyc1              [N x 1 vector]          : length (cycles) for the bursts of 'f1' frequency

    % REFERENCES
    % [1] Quinn A.J., van Ede F., Brookes M.J., Heideman S.G., Nowak M., Seedat Z.A., Vidaurre D., Zich C., Nobre A.C., Woolrich M.W. (2019).
    %     Unpacking Transient Event Dynamics in  Electrophysiological Power Spectra. Brain Topography. https://doi.org/10.1007/s10548-019-00745-5.
    %     MIT License, Copyright (c) 2019 OHBA Analysis Group
    %     Partly adapted from "Quinn2019_BurstHMM/hmm_util_get_simulation.m"

    % Written by SungJun Cho, April 1, 2021
    % Last modified on October 29, 2021
    %% Configure Input Parameters
    % [1] Decide plotting availability
    if nargin < 10
        verbose = false;
    end
    % [2] Assign random seed number
    if nargin < 9
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

    % [2] Set Simulation Parameters
    Fs = 512; % sampling rate
    Nyq = Fs/2; % Nyquist frequency
    tvec = linspace(0,T,T*Fs)'; % time vector

    % [3] Create a noisy time course
    % The time series will approximately follow 1/f power scale by direct pole-placement
    a = -poly(.92); % the filter root; the denominator polynomial of a digital filter
    % The roots given by "poly" should be between 0 <= x <= 1 where 0 is white noise and
    % 1 is a spectrum with extreme slope.
    pink_noise = 0.2*randn(size(tvec));
    x_noise = filter(1,a,pink_noise);

    % [4] Define burst occurrences and durations
    replacement = true;

    start1 = sort(randsample(round(Nyq)+1:length(tvec)-(round(Nyq)+1)-(rng_pks(2)*(Fs/f1)),nBurst1)); % put ':Fs:' to enforce random burst intervals with at least 1s gap
    available_duration1 = round((rng_pks(1):rng_pks(2)).*(Fs/f1));
    duration1 = randsample(available_duration1,nBurst1,replacement); % with replacement

    start2 = sort(randsample(round(Nyq)+1:length(tvec)-round(Nyq)+1-(rng_pks(2)*(Fs/f2)),nBurst2));
    available_duration2 = round((rng_pks(1):rng_pks(2)).*(Fs/f2));
    duration2 = randsample(available_duration2,nBurst2,replacement); % with replacement

    % Note: Nyquist frequency implemented here to avoid generating bursts at the area (i.e., at the start and the end) uncovered by 
    % the spectrogram computation (under the assumption that the window is less than 1s).

    if length(start1) ~= nBurst1 || length(start2) ~= nBurst2
        error('ComputationError: The number of generated bursts does not match the preset number.');
    end

    for t = 1:length(start1)
        % Consider the case in which the predefined burst duration exceeds the signal length
        if t == length(start1)
            if start1(t)+duration1(t) > length(tvec)
                duration1(t) = length(tvec) - start1(t);
            end
            break
        end
        % Consider the case in which the next burst overwrites the previous burst
        % Reinforce the distance of 2 cycles as a gap between subsequent bursts in this case
        if start1(t)+duration1(t) > start1(t+1)
            start1(t+1) = start1(t+1) + (start1(t)+duration1(t)-start1(t+1))+round(2*Fs/f1);
        end
    end

    for t = 1:length(start2)
        % Consider the case in which the predefined burst duration exceeds the signal length
        if t == length(start2)
            if start2(t)+duration2(t) > length(tvec)
                duration2(t) = length(tvec) - start2(t);
            end
            break
        end
        % Consider the case in which the next burst overwrites the previous burst
        % Reinforce the distance of 2 cycles as a gap between subsequent bursts in this case
        if start2(t)+duration2(t) > start2(t+1)
            start2(t+1) = start2(t+1) + (start2(t)+duration2(t)-start2(t+1))+round(2*Fs/f2);
        end
    end

    % [5] Create burst-embedded time-series
    x_f1 = zeros(T*Fs,1);
    x_f2 = zeros(T*Fs,1);

    switch method
        case 'tukey'
            r = 0.25; % cosine fraction; default = 0.5
                      % r = 0 -> a rectangular window
                      % r = 1 -> equivalent to a Hanning window
            for t = 1:length(start1)
                % Add slow burst
                temp_sig = sin(2*pi*f1*tvec(start1(t):start1(t)+duration1(t)));
                temp_sig = temp_sig.*tukeywin(length(temp_sig),r);
                x_f1(start1(t):start1(t)+duration1(t)) = temp_sig;
            end
            for t = 1:length(start2)
                % Add fast burst
                temp_sig = sin(2*pi*f2*tvec(start2(t):start2(t)+duration2(t)));
                temp_sig = temp_sig.*tukeywin(length(temp_sig),r);
                x_f2(start2(t):start2(t)+duration2(t)) = temp_sig;
            end
            %wvtool(tukeywin(length(temp_sig),r)); % display Tukey (tapered cosine) window
        case 'gauss'
            alpha = 1.2; % width factor; default = 2.5
                         % selected such that when the duration is about 3 cycles, all the
                         % peaks show sufficient magnitude of amplitude
            for t = 1:length(start1)
                % Add slow burst
                temp_sig = sin(2*pi*f1*tvec(start1(t):start1(t)+duration1(t)));
                temp_sig = temp_sig.*gausswin(length(temp_sig),alpha);
                x_f1(start1(t):start1(t)+duration1(t)) = temp_sig;
            end
            for t = 1:length(start2)
                % Add fast burst
                temp_sig = sin(2*pi*f2*tvec(start2(t):start2(t)+duration2(t)));
                temp_sig = temp_sig.*gausswin(length(temp_sig),alpha);
                x_f2(start2(t):start2(t)+duration2(t)) = temp_sig;
            end
            %wvtool(gausswin(length(temp_sig),alpha)); % display Gaussian window
    end

    % [6] Add noise and burst time-series to generate output signal
    amp_f1 = 0.630;
    amp_f2 = 0.551;
    x_f1 = amp_f1.*x_f1;
    x_f2 = amp_f2.*x_f2;
    simulated_signal = x_noise + x_f1 + x_f2; % superpose each signal component

    % [7] Add Gaussian white noise
    gwn = randn(size(simulated_signal)) * (amp_f1/db2mag(SNR_wn));
    simulated_signal = simulated_signal + gwn;
    % References: (1) https://www.mathworks.com/help/matlab/ref/randn.html
    %             (2) https://www.mathworks.com/help/matlab/math/random-numbers-with-specific-mean-and-variance.html
    %             (3) https://www.mathworks.com/matlabcentral/answers/472495-generating-random-variable-from-certain-standard-deviation-and-mean
    %             (4) https://en.wikipedia.org/wiki/Additive_white_Gaussian_noise

    % [8] Concatenate different parts of simulation
    data = cat(2,simulated_signal,x_f1,x_f2,x_noise,gwn);
    %% Get Burst Information of Slow Bursts
    btime_f1 = cell(length(start1),1);
    burst_f1 = cell(length(start1),1);
    for t = 1:length(start1)
        btime_f1{t} = tvec(start1(t):start1(t)+duration1(t));
        burst_f1{t} = x_f1(start1(t):start1(t)+duration1(t));
    end
    cyc1 = round(duration1./(Fs/f1))'; % length of each burst in cycles
    %% Plot Simulations
    if verbose
        text_title = {'Complete Data',['Burst #1 (f = ' num2str(f1) 'Hz)'],['Burst #2 (f = ' num2str(f2) 'Hz)']};
        clr_array = [35,99,119; 124,173,62; 74,108,47]/255;
        figure();
        for i = 1:length(text_title)
            subplot(3,1,i);
            p = plot(tvec,data(:,i),'Color',clr_array(i,:));
            p.Color(4) = 0.7;
            xlabel('Time (s)');
            ylabel('Amplitude (a.u.)');
            title(text_title{i},'FontWeight','bold');
            set(gca,'FontSize',12);
        end
    end
end