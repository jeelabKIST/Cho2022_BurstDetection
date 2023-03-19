function [freq, time, powermat, coi] = spectrogram_cwt(time, X, Fs, wname, scale_opt, verbose)
    %% Function 'spectrogram_cwt'
    % DESCRIPTION
    % This function computes and plots a continuous wavelet spectrogram (i.e. a
    % magnitude spectrogram) of an input data. The current version only
    % supports the single-trial CWT spectrogram; that is, this code
    % works only for an input data of a single vector array as of now.

    % USAGE
    % Full Input : spectrogram_cwt(time, X, Fs, wname, scale_opt, verbose)
    % Example    : spectrogram_cwt([], data, 1024, 'Morse', 'expand', true)
    %            : spectrogram_cwt(time, data, 1024)

    % INPUT
    %    Variable     Data Type                 Description
    % 1. time         [1 x N vector]          : time array
    %                                           1) []
    %                                           2) vector array
    % 2. X            [1 x N vector]          : input data
    % 3. Fs           [integer Z > 0]         : sampling frequency rate
    % 4. wname        [string]                : type of analysis wavelet
    %                                         : 'Morse', 'amor' (default), 'bump'
    % 5. scale_opt    [string]                : option to define scales
    %                                           1) 'default' : the minimum and maximum scales will be
    %                                                          determine automatically by the function 'cwt'
    %                                           2) 'expand'  : uses 'cwtfreqbounds' to deterime the min and max 
    %                                                          wavelet bandpass freqencies for given signal length, 
    %                                                          sampling frequency, and wavelet
    % 6. verbose      [logical]               : display CWT spectrogram and its summary
    %                                           1) true or 1
    %                                           2) false or 0

    % OUTPUT
    %    Variable     Data Type                 Description
    % 1. freq         [1 x N vector]          : frequency vector of CWT
    % 2. time         [1 x N vector]          : time vector of CWT
    % 3. powermat     [N x N matrix]          : power of wavelet coefficients
    % 4. coi          [N x 1 vector]          : cone of influence which indicates the locations where edge effects occur

    % REFERENCE
    % [1] MATLAB Documentation. https://www.mathworks.com/help/wavelet/ref/cwt.html
    % [2] MATLAB Documentaiton. https://www.mathworks.com/help/wavelet/ref/cwtfreqbounds.html
    % [3] Ahmet Taspinar. (2018). A guide for using the Wavelet Transform in Machine Learning. Retrieved from 
    %     http://ataspinar.com/2018/12/21/a-guide-for-using-the-wavelet-transform-in-machine-learning/ on Dec 11, 2020.

    % Note: As indicated by MATLAB,the plot uses a logarithmic frequency axis (in powers of 10) since frequencies in the CWT are logarithmic.

    % Written by SungJun Cho, March 25, 2021
    % Last modified on October 29, 2021
    %% Input Parameters
    if nargin < 6
        verbose = true;
    end
    if nargin < 5
        scale_opt = 'default';
    end
    if nargin < 4
        wname = 'amor'; % analytic Morlet wavelet
    end
    % [1] Check Input Type
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
    % [2] Create Time Array
    if isempty(time)
        time = linspace(0, length(X)/Fs, length(X));
    end
    %% Calculate CWT Spectrogram
    switch scale_opt
        case 'default'
            [wtcoefs,freq,coi] = cwt(X,wname,Fs);
        case 'expand'
            [minf,maxf] = cwtfreqbounds(numel(X),Fs);
            fb = cwtfilterbank('SignalLength', numel(X), 'Wavelet', wname, 'SamplingFrequency', Fs,'FrequencyLimits',[minf,maxf]);
            [wtcoefs,freq,coi] = cwt(X,'FilterBank',fb);
            %figure(); freqz(fb);
    end
    powermat = abs(wtcoefs); % Or alternatively: abs(wtcoefs).^2;
    %                                            wtcoefs.*conj(wtcoefs);
    % wtcoefs: continuous wavelet transform of the input signal for specified scales and wavelet (returned as a complex matrix)
    %% Plot CWT Spectrogram
    % You can use the function 'surface' as suggested in MATLAB document, 
    % but note that it is slower than 'imagesc'.
    if verbose
        switch scale_opt
            case 'default'
                figure(); hold on;
                imagesc(time,freq,imgaussfilt(powermat,1));
                shading flat; colormap('jet');
                plot(time,coi,'--w','LineWidth',2);
                xlim([time(1) time(end)]);
                ylim([1 max(freq)]);
                ftemp = flip(freq);
                ax = gca;
                ax.YTick = round(ftemp(1:10:end),3);
                set(ax,'yscale','log');
                axis xy;
                xlabel('Time (s)');
                ylabel('Frequency (Hz)');
                title('Wavelet Spectrogram / Scalogram');
                c = colorbar;
                ylabel(c,'Coefficient Magnitude');
            case 'expand'
                figure(); hold on;
                fBin = 15;
                freq_log = logspace(log10(minf),log10(maxf),fBin);
                imagesc(time,freq,imgaussfilt(powermat,1));
                shading flat; colormap('jet');
                plot(time,coi,'--w','LineWidth',2);
                xlim([time(1) time(end)]);
                ylim([1 max(freq_log)]);
                ax = gca;
                set(ax,'yscale','log');
                ax.YTick = freq_log;
                axis xy;
                xlabel('Time (s)');
                ylabel('Frequency (Hz)');
                title('Wavelet Spectrogram / Scalogram');
                c = colorbar;
                ylabel(c,'Coefficient Magnitude');
        end
    end
    % Display spectrogram properties
    if verbose
        switch scale_opt
            case 'default'
                display_summary(wname,scale_opt,freq);
            case 'expand'
                display_summary(wname,scale_opt,freq,fb);
        end
    end
end

%% Appendix: In-Script Functions
% Function #1: Displays a summary of CWT sepctrogram parameters and
% resoultions
function display_summary(wname, verbose, freq, fb)
    if nargin < 4
        fb = [];
    end
    disp('*****************SUMMARY*****************');
    disp('CWT Spectrogram Parameters:');
    disp(['    Type of Wavelet: ' wname]);
    switch verbose
        case 'default'
            disp(['    Frequency Limits: ' num2str(min(freq)) 'Hz-' num2str(max(freq)) 'Hz']);
            disp(['    Number of Voices Per Octave: ' num2str(10)]);
        case 'expand'
            disp(['    Frequency Limits: ' num2str(fb.FrequencyLimits(1)) 'Hz-' num2str(fb.FrequencyLimits(2)) 'Hz']);
            disp(['    Number of Voices Per Octave: ' num2str(fb.VoicesPerOctave)]);
    end
    disp('*****************************************');
end