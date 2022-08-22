function [normSpec, muSpec, stdSpec] = normalize_spectrogram(type, method, window_scale, Spec_f, Spec_t, Spec, Fs, verbose)
    %% Function 'normalize_spectrogram'
    % DESCRIPTION
    % This function normalizes the input spectrogram on the basis of single frequency bands.

    % USAGE
    % Full Input : [normSpec, muSpec, stdSpec] = normalize_spectrogram(type, method, window_scale, Spec_f, Spec_t, Spec, Fs, verbose)
    % Example    : [normSpec,~,~] = normalize_spectrogram('FT', 'zscore', 'equal', Spec_f, Spec_t, Spec, Fs, 2);
    %            : [normSpec,~,~] = normalize_spectrogram('WT', 'minmax_scale', 'variable', freq, time, powermat', Fs);

    % INPUT
    %    Variable       Data Type                    Description
    % 1. type           [string]                   : type of spectrogram
    %                                                1) FT: Fourier transfrom based
    %                                                2) WT: Wavelet transfrom based
    % 2. method         [string]                   : time array
    %                                                1) minmax_scale
    %                                                2) minmax_mean
    %                                                3) zscore
    %                                                4) log (do not recommend for 'WT' as frequency is already in a logarithmic scale)
    % 3. window_scale   [string]                   : indicator of window size used
    %                                                1) equal (if Spec_f is uniformly distributed and fixed window was used for computation)
    %                                                2) variable (if Spec_f is not uniformly distributed and frequency-dependent window was used for computation)
    % 4. Spec_f         [1 x N vector]             : frequency vector of an input spectrogram
    % 5. Spec_t         [1 x N vector]             : time vector of an input spectrogram
    % 6. Spec           [N x N matrix]             : input power spectrogram
    % 7. Fs             [integer Z > 0]            : sampling frequency
    % 8. verbose        [{0, 1, 2}]                : type of the plot
    %                                                1) 0 - do not plot
    %                                                2) 1 - plot without smoothing (not recommended for 'FT'-'log')
    %                                                3) 2 - plot with smoothing

    % OUTPUT
    %    Variable       Data Type                    Description
    % 1. normSpec       [N x N matrix]             : the normalized spectrogram
    % 2. muSpec         [1 x N vector]             : mean values of the spectrogram at each frequency
    % 3. stdSpec        [1 x N matrix]             : standard deviation values of the spectrogram at each frequency

    % Written by SungJun Cho, March 26, 2021
    % Last modified on October 29, 2021
    %% Input Parameters
    if nargin < 8
        verbose = 0;
    end
    %% Initialize Empty Arrays
    muSpec = zeros(1,length(Spec_f));
    stdSpec = zeros(1,length(Spec_f));
    tempSpec = cell(1,length(Spec_f));
    cropIdx = zeros(2,length(Spec_f));
    normSpec = zeros(length(Spec_t),length(Spec_f));
    %% Mean and Standard Deviation of Spectrogram at Each Frequency Level
    for i = 1:length(Spec_f)
        sel_Spec = Spec(:,i)';
        sel_Spec(sel_Spec~=0) = 1;
        start = strfind([0 sel_Spec], [0 1]);
        stop = strfind([sel_Spec 0], [1 0]);
        if length(start) > 1 || length(stop) > 1
            start = start(1); stop = stop(end);
        end
        muSpec(1,i) = mean(Spec(start:stop,i)); % mean of a spectrogram at each frequency level
        stdSpec(1,i) = std(Spec(start:stop,i)); % standard deviation of a spectrogram at each frequency level
        tempSpec{1,i} = Spec(start:stop,i);
        cropIdx(:,i) = [start;stop];
    end
    %% Normalization
    switch method
        case 'minmax_scale'
            labelY = 'MinMax-Scale';
            for i = 1:length(Spec_f)
                temp_spec = tempSpec{1,i};
                start = cropIdx(1,i); stop = cropIdx(2,i);
                normSpec(start:stop,i) = (temp_spec - min(temp_spec)) ./ (max(temp_spec)-min(temp_spec));
            end
        case 'minmax_mean'
            labelY = 'MinMax-Mean';
            for i = 1:length(Spec_f)
                temp_spec = tempSpec{1,i};
                start = cropIdx(1,i); stop = cropIdx(2,i);
                normSpec(start:stop,i) = (temp_spec - muSpec(1,i)) ./ (max(temp_spec) - min(temp_spec));
            end
        case 'zscore'
            labelY = 'Z-Score';
            for i = 1:length(Spec_f)
                temp_spec = tempSpec{1,i};
                start = cropIdx(1,i); stop = cropIdx(2,i);
                normSpec(start:stop,i) = (temp_spec - muSpec(1,i)) ./ stdSpec(1,i);
            end
        case 'log'
            labelY = 'dB';
            for i = 1:length(Spec_f)
                temp_spec = tempSpec{1,i};
                start = cropIdx(1,i); stop = cropIdx(2,i);
                normSpec(start:stop,i) = pow2db(temp_spec);
            end
    end
    %% Plot the Normalized Spectrogram
    if verbose ~= 0
        switch type
            case 'FT'
                if strcmp(window_scale,'equal')
                    figure();
                    if verbose == 1
                        imagesc(Spec_t,Spec_f, normSpec');
                    elseif verbose == 2 % additional Gaussian filter
                        imagesc(Spec_t, Spec_f, imgaussfilt(normSpec',2));
                    end
                    axis xy; axis tight;
                    xlabel('Time (s)');
                    ylabel('Frequency (Hz)');
                    title('Normalized FT-Based Spectrogram');
                    c = colorbar; colormap(jet);
                    ylabel(c, ['Normalized Power (' labelY ')']);
                elseif strcmp(window_scale,'variable')
                    if verbose == 1
                        figure();
                        [~,obj] = contourf(Spec_t,Spec_f,normSpec');
                        obj.LineStyle = 'none';
                        set(gca,'YTick',round(Spec_f,2));
                        xlabel('Time (s)');
                        ylabel('Frequency (Hz)');
                        title('Normalized FT-Based Spectrogram');
                        c = colorbar; colormap(jet);
                        ylabel(c, ['Normalized Power (' labelY ')']);
                        % Interpolation
                        interp_method = 'spline';
                        Spec_f_ip = Spec_f(1):Spec_f(end);
                        Spec_t_ip = Spec_t;
                        [Xq, Yq] = meshgrid(Spec_f_ip,Spec_t_ip);
                        SpecInterp = interp2(Spec_f,Spec_t,normSpec,Xq,Yq,interp_method);
                        figure();
                        imagesc(Spec_t_ip,Spec_f_ip,SpecInterp'); axis xy;
                        xlabel('Time (s)');
                        ylabel('Frequency (Hz)');
                        title('Interpolated Normalized FT-Based Spectrogram');
                        c = colorbar; colormap(jet);
                        ylabel(c, ['Normalized Power (' labelY ')']);
                    elseif verbose == 2 % additional Gaussian filter
                        figure();
                        [~,obj] = contourf(Spec_t,Spec_f,imgaussfilt(normSpec',2));
                        obj.LineStyle = 'none';
                        set(gca,'YTick',round(Spec_f,2));
                        xlabel('Time (s)');
                        ylabel('Frequency (Hz)');
                        title('Normalized FT-Based Spectrogram');
                        c = colorbar; colormap(jet);
                        ylabel(c, ['Normalized Power (' labelY ')']);
                        % Interpolation
                        interp_method = 'spline';
                        Spec_f_ip = Spec_f(1):Spec_f(end);
                        Spec_t_ip = Spec_t;
                        [Xq, Yq] = meshgrid(Spec_f_ip,Spec_t_ip);
                        SpecInterp = interp2(Spec_f,Spec_t,normSpec,Xq,Yq,interp_method);
                        figure();
                        imagesc(Spec_t_ip,Spec_f_ip,imgaussfilt(SpecInterp',2)); axis xy;
                        xlabel('Time (s)');
                        ylabel('Frequency (Hz)');
                        title('Interpolated Normalized FT-Based Spectrogram');
                        c = colorbar; colormap(jet);
                        ylabel(c, ['Normalized Power (' labelY ')']);
                    end
                end
            case 'WT'
                if ~strcmp(window_scale,'variable')
                    warning('The frequency scale of wavelet scalogram is not equally distributed. The variable window scale reset to "variable". \n');
                    window_scale = 'variable';
                end
                if strcmp(window_scale,'variable')
                    figure();
                    if verbose == 1
                        imagesc(Spec_t, Spec_f, normSpec');
                    elseif verbose == 2
                        imagesc(Spec_t, Spec_f, imgaussfilt(normSpec',2));
                    end
                    ylim([1 max(Spec_f)+1]);
                    shading flat;
                    ax = gca;
                    [minf,maxf] = cwtfreqbounds(numel(Spec_t),Fs);
                    if round(maxf) == round(max(Spec_f))
                        fBin = 15;
                        freq_log = logspace(log10(minf),log10(maxf),fBin);
                        ax.YTick = freq_log;
                    else
                        ftemp = flip(Spec_f);
                        ax.YTick = round(ftemp(1:10:end),3);
                    end
                    set(ax,'yscale','log');
                    axis xy;
                    xlabel('Time (s)');
                    ylabel('Frequency (Hz)');
                    title('Normalized Spectrogram');
                    c = colorbar; colormap(jet);
                    ylabel(c, ['Normalized Power (', labelY ')']);
                end
        end
    end
end