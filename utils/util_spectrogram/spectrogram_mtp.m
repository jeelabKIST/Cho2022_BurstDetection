function [Spec_f, Spec_t, SpecMat, fft_interval] = spectrogram_mtp(time, X, Fs, nw, ntp, freq_range, dt, nOscCycle, interp_opt, plot_psd, check_opt, plot_tp, verbose)
    %% Function 'spectrogram_mtp'
    % DESCRIPTION
    % This function computes and plots a multitaper spectrogram of an input
    % data. The current version only supports the single-trial multitaper 
    % spectrogram (i.e. works only for an input data of a single vector array).

    % USAGE
    % Full Input   : spectrogram_mtp(time, X, Fs, nw, ntp, freq_range, dt, nOscCycle, interp_opt, plot_psd, check_opt, plot_tp, verbose);
    % Example      : spectrogram_mtp(time, data, 1024, 2.5, 4, 1/Fs, 6, true, true, 'method1', true, true);

    % INPUT
    %    Variable     Data Type                   Description
    % 1. time         [1 x N vector]            : time array
    %                                             1) []
    %                                             2) vector array
    % 2. X            [1 x N vector]            : input data
    % 3. Fs           [integer Z > 0]           : sampling frequency rate
    % 4. nw           [integer Z > 0]           : the time-half bandwidth for the tapers
    %                                           : 2.5, 3, 3.5, 4 recommended
    % 5. ntp          [integer Z > 0]           : the number of tapers
    % 6. freq_range   [1 x N vector]            : the array of frequencies of interest
    % 7. dt           [rational Q > 0]          : time step by which the window slides (seconds)
    %                                           : default) 0.1
    % 8. nOscCycle    [integer Z > 0]           : number of oscillatory cycles to define the window sizes
    %                                             1) 0 for the window size fixed to one second
    %                                             2) > 0 for the frequency-dependent window sizes
    % 9. interp_opt   [logical]                 : option to interpolate 2-D gridded spectrogram
    %                                             Note) counts only when nOscCycle > 0
    %                                             1) true or 1
    %                                             2) false or 0
    % 10. plot_psd    [boolean / logical]       : display multitaper power spectral density estimates
    %                                             1) true or 1
    %                                             2) false or 0
    % 11. check_opt   [string]                  : check optimal number of tapers based on user's selected time-half bandwidth
    %                                             1) []
    %                                             2) 'method1'
    %                                             3) 'method2'
    % 12. plot_tp     [boolean / logical]       : display DPSS tapers
    %                                             1) true or 1
    %                                             2) false or 0
    % 13. verbose     [boolean / logical]       : option to plot the computed spectrogram and display summary
    %                                             1) true or 1
    %                                             2) false or 0

    % OUTPUT
    %    Variable         Data Type               Description
    % 1. Spec_f           [1 x N vector]        : vector of frequency bins for the spectrogram
    % 2. Spec_t           [1 x N vector]        : vector of times for the spectrogram
    % 3. SpecMat          [N x N matrix]        : matrix of calculated spectral power
    % 4. fft_interval     [N x 2 matrix]        : matrix of frequency ranges centered around each element of Spec_f

    % REFERENCE
    % [1] Prerau M.J., Brown R.E., Bianchi M.T., Ellenbogen J.M., Purdon P.L. (2016). Sleep Neurophysiological Dynamics
    %     Through the Lens of Multitaper Spectral Analysis. Physiology, 32(1): 60-92. DOI: 10.1152/physiol.00062.2015.
    %     Copyright 2019 Michael J. Prerau, Ph.D., Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License (adapted the idea of parallel computing)
    % [2] Percival D.B., Walden A.T. (1993). Spectral Analysis for Physical Applications. UK: Cambridge University Press.
    % [3] Chronux Analysis Software. http://chronux.org/chronuxFiles/Documentation/chronux/spectral_analysis/helper/dpsschk.html
    % [4] Chronux Analysis Software. http://chronux.org/chronuxFiles/Documentation/chronux/spectral_analysis/continuous/mtspectrumc.html
    % [5] MATLAB Documentation. https://www.mathworks.com/help/signal/ref/dpss.html

    % Written by SungJun Cho, March 25, 2021
    % Last modified on February 07, 2023
    %% Input Parameter Configuration
    if nargin < 13
        verbose = true;
    end
    if nargin < 10
        plot_psd = true;
        check_opt = 'method1';
        plot_tp = true;
    end
    if nargin < 9
        interp_opt = true;
    end
    if nargin < 8
        nOscCycle = 0;
    end
    if nargin < 7
        dt = 0.1;
    end
    if nargin < 6
        error('INPUT ERROR: Must provide adequate inputs for the function.');
    end
    if nOscCycle ~= 0
        plot_tp = false;
        if freq_range(1) == 0
            warning('Frequency-dependent window cannot be defined for DC frequency. The lowest frequency will be reset to 1 Hz.');
            freq_range(1) = 1;
        end
    else
        interp_opt = [];
    end
    if freq_range(end) > Fs/2
        warning('The highest frequency should be defined smaller than the Nyquist frequency and will be reset to the corresponding Nyquist frequency.');
        freq_range(end) = Fs/2;
    end
    % [1] Create Time Array
    if isempty(time)
        time = linspace(0, length(X)/Fs, length(X));
    end
    if length(time) ~= length(X)
        error('INPUT ERROR: The sizes of a time and data array should match.');
    end
    % [2] Check Input Type
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
    % [3] Plot Multitaper Power Spectral Density Estimate
    dpss_params = {nw, ntp};
    if plot_psd
        figure();
        pmtm(X, dpss_params, length(X), Fs, 'DropLastTaper', false);
    end
    % [4] Define Time Arrays
    moving_step = floor(Fs*dt); % time step in sampling points for sliding windows
    time_new = linspace(time(1),time(end),length(time)/moving_step);
    %% Compute Fixed Window Multitaper Spectrogram
    %fprintf('Computing a multitaper spectrogram ... '); tic;

    if nOscCycle == 0
        window_size = Fs; % 1s fixed window
        % [1] Create Tapers for Fixed Window
        [dpss_win, eigvals] = dpss(window_size, dpss_params{1}, dpss_params{2});
        % [2] Check for the Optimal Number of the Tapers
        if ~isempty(check_opt)
            switch check_opt
                case 'method1'
                    if floor(2*nw-1) < ntp
                        opt = floor(2*nw-1);
                        warning(['The optimal number of tapers for the chosen time-half bandwidth is: ', num2str(opt)]);
                        ntp = opt;
                        fprintf('The number of tapers is set to its optimal value. \n');
                    end
                case 'method2'
                    figure();
                    stem(1:length(eigvals), eigvals, 'filled');
                    hold on;
                    plot(1:length(eigvals), 0.99*ones(length(eigvals),1));
                    ylim([0 1.2]);
                    title('Proportion of Energy in [-w,w] of k-th Slepian Sequence');
                    opt = find(eigvals>0.99,1,'last');
                    if ntp ~= opt
                        warning(['The optimal number of tapers for the chosen time-half bandwidth is: ', num2str(opt)]);
                        ntp = opt;
                        fprintf('The number of tapers is set to its optimal value. \n');
                    end
            end
        end
        dpss_win = dpss_win * sqrt(Fs); % for the normalization of the tapers (cf. [3])
        % [3] Plot DPSS Tapers
        if plot_tp
            figure(); hold on;
            for i = 1:ntp
                txt = ['win #', num2str(i)];
                plot(dpss_win(:,i), 'DisplayName', txt);
            end
            xlabel('Samples');
            title(['DPSS, M=', num2str(size(dpss_win,1)), ', NW=', num2str(nw)]);
            xlim([0 size(dpss_win,1)]);
            legend show;
        end
        % [4] Get Sliding Window Indices
        idx_start = 1+floor(window_size/2);
        idx_end = length(time) - floor(window_size/2);
        idx_midpts = idx_start:moving_step:idx_end;

        idx_collection = zeros(length(idx_midpts),2);
        for idx = 1:length(idx_midpts)
            idx_collection(idx,1) = idx_midpts(idx) - window_size*0.5;
            idx_collection(idx,2) = idx_midpts(idx) + window_size*0.5 - 1;
        end

        idx_diff = idx_collection(:,2) - idx_collection(:,1);
        if idx_diff ~= (window_size - 1)
            error('ERROR: Index computations for sliding windows are incorrect.');
        end
        % [5] Compute Scaled Frequency Vector
        nfft = max(2^(nextpow2(window_size)),window_size);
        df = Fs/nfft;
        f = 0:df:Fs/2;
        select_freq = find(f >= freq_range(1) & f <= freq_range(end));
        freq_vec = f(select_freq);
        % [6] Compute Time Vector
        time_temp = time(idx_midpts);
        cut_start = find(time_new < time_temp(1), 1, 'last');
        cut_end = find(time_new > time_temp(end), 1, 'first');    
        if isempty(cut_start)
            cut_start = 0;
        end
        if isempty(cut_end)
            cut_end = length(time_new) + 1;
        end
        t_start = cut_start + 1;
        t_end = cut_end - 1;
        % [7] Temporal Smoothing With Fixed Window Length
        Spec_temp = zeros(length(time_new),length(freq_vec),ntp);
        data_segments = zeros(length(time_new),window_size);
        for tIdx = 1:length(t_start:t_end)
            data_segments(tIdx,:) = X(idx_collection(tIdx,1):idx_collection(tIdx,2));
        end
        for time_idx = 1:length(t_start:t_end)
            x_seg = data_segments(time_idx,:);
            x_ipt = repmat(x_seg',1,ntp).*dpss_win;
            x_mtp = fft(x_ipt,nfft)/Fs;
            abs_val = conj(x_mtp).*x_mtp;
            Spec_temp(time_idx+t_start-1,:,:) = abs_val(select_freq,:);
            % or alternatively: abs(x_mtp(select_freq,:));
        end
        Spec_taper_mean = mean(Spec_temp,3);
        % [8] Store Results
        fft_interval = zeros(length(freq_vec),2);
        for j = 1:length(freq_vec)
            fft_interval(j,:) = [freq_vec(j)-(df/2) freq_vec(j)+(df/2)];
        end
        % [9] Output Results
        Spec_t = time_new;
        Spec_f = freq_vec;
        SpecMat = Spec_taper_mean;
    end
    %% Compute Frequency-Dependent Window Multitaper Spectrogram
    if nOscCycle > 0
        % [1] Define Window Sizes
        window_sizes = zeros(length(freq_range),1);
        dfs = zeros(length(freq_range),1);
        for w = 1:size(window_sizes,1)
            window_sizes(w) = round(nOscCycle*(Fs/freq_range(w)));
            if mod(window_sizes(w),2) ~= 0
                window_sizes(w) = window_sizes(w) + 1;
                % ensure windows to be of even length to prevent non-integer
                % indexing error
            end
            dfs(w,1) = Fs/window_sizes(w);
        end
        % [2] Preallocate Cell and Double Arrays for the Final Outputs
        SpecCell = cell(1,size(window_sizes,1));
        TimeCell = cell(1,size(window_sizes,1));
        IndCell = cell(1,size(window_sizes,1));
        freq_vec = zeros(1,size(window_sizes,1));
        fft_interval = zeros(size(window_sizes,1),2);
        for w = 1:size(window_sizes,1)
            % [3] Construct DPSS Tapers
            window_size = window_sizes(w);
            [dpss_win, eigvals] = dpss(window_size, nw, ntp);
            % [4] Check for the Optimal Number of the Tapers
            if ~isempty(check_opt)
                switch check_opt
                    case 'method1'
                        if floor(2*nw-1) < ntp
                            opt = floor(2*nw-1);
                            warning(['The optimal number of tapers for the chosen time-half bandwidth is: ', num2str(opt)]);
                            ntp = opt;
                            fprintf('The number of tapers is set to its optimal value. \n');
                        end
                    case 'method2'
                        figure();
                        stem(1:length(eigvals), eigvals, 'filled');
                        hold on;
                        plot(1:length(eigvals), 0.99*ones(length(eigvals),1));
                        ylim([0 1.2]);
                        title('Proportion of Energy in [-w,w] of k-th Slepian Sequence');
                        opt = find(eigvals>0.99,1,'last');
                        if ntp ~= opt
                            warning(['The optimal number of tapers for the chosen time-half bandwidth is: ', num2str(opt)]);
                            ntp = opt;
                            fprintf('The number of tapers is set to its optimal value. \n');
                        end
                end
            end
            dpss_win = dpss_win * sqrt(Fs); % for the normalization of the tapers (cf. [3])
            % [5] Get Sliding Window Indices
            idx_start = 1 + floor(window_size/2);
            idx_end = length(time) - floor(window_size/2);
            idx_midpts = idx_start:moving_step:idx_end;
            idx_collection = zeros(length(idx_midpts),2);
            for idx = 1:length(idx_midpts)
                idx_collection(idx,1) = idx_midpts(idx) - window_size*0.5;
                idx_collection(idx,2) = idx_midpts(idx) + window_size*0.5 - 1;
            end
            idx_diff = idx_collection(:,2) - idx_collection(:,1);
            if idx_diff ~= (window_size - 1)
                error('ERROR: Index computations for sliding windows are incorrect.');
            end
            % [6] Compute Scaled Frequency Vector
            f_org = (Fs*(0:window_size/2)/window_size);
            f = round(f_org);
            if f ~= round(0:dfs(w):Fs/2)
                error('ERROR: Frequency vector computation has encountered an error.');
            end
            subtracted = abs(f_org - freq_range(w));
            select_freq = find(subtracted == min(subtracted));
            freq_vec(1,w) = f_org(select_freq);
            % [7] Pad Zeros to Equalize Different Time Vector Lengths
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
                cut_end = length(time_new) + 1;
            end
            t_start = cut_start + 1;
            t_end = cut_end - 1;
            % [8] Temporal Smoothing with Frequency-Dependent Window Length
            nfft = window_size;
            Spec_temp = zeros(length(time_new),ntp);
            for tIdx = 1:length(t_start:t_end)
                x_seg = X(idx_collection(tIdx,1):idx_collection(tIdx,2));
                x_mtp = repmat(x_seg',1,ntp).*dpss_win;
                x_fft = fft(x_mtp,nfft)/Fs;
                abs_val = conj(x_fft).*x_fft;
                Spec_temp(tIdx+t_start-1,:) = abs_val(select_freq,:);
                % or alternatively: abs(x_fft(select_freq,:));
            end
            Spec_taper_mean = mean(Spec_temp,2);
            % [9] Store Results
            SpecCell{w} = Spec_taper_mean;
            TimeCell{w} = t_vec;
            IndCell{w} = idx_collection;
            fft_interval(w,:) = [f_org(select_freq)-(dfs(w)/2) f_org(select_freq)+(dfs(w)/2)];
        end
        % [10] Output Results
        Spec_t = time_new;
        Spec_f = freq_vec;
        SpecMat = [SpecCell{:}];
    end
    %% Plot Figure
    if ~verbose
        if length(Spec_f) ~= length(unique(Spec_f))
            error('Redundant frequency value occurred during estimation of scaled frequency vector.');
        end
    else
        if length(Spec_f) == length(unique(Spec_f))
            if nOscCycle > 0
                figure();
                [~,obj] = contourf(Spec_t,Spec_f,SpecMat');
                obj.LineStyle = 'none';
                set(gca,'YTick',round(Spec_f,2));
                xlabel('Time (s)');
                ylabel('Frequency (Hz)');
                title('Multitaper Spectrogram (f-Dependent Window)');
                colormap('jet'); c = colorbar;
                ylabel(c, 'Power (mV^2)','FontSize', 12);
                if interp_opt
                    interp_method = 'spline'; % or alternatively: 'linear'
                    Spec_f_ip = Spec_f(1):Spec_f(end);
                    Spec_t_ip = Spec_t;
                    [Xq,Yq] = meshgrid(Spec_f_ip,Spec_t_ip);
                    SpecInterp = interp2(Spec_f,Spec_t,SpecMat,Xq,Yq,interp_method);
                    % visualized center frequency may deviate from the actual center
                    % frequency after interpolation
                    figure();
                    imagesc(Spec_t_ip,Spec_f_ip,imgaussfilt(SpecInterp',2));
                    % spectral power may decrease when Gaussian filtered
                    axis xy;
                    xlabel('Time (s)');
                    ylabel('Frequency (Hz)');
                    title('Multitaper Spectrogram (Frequency Interpolated)');
                    c = colorbar; colormap('jet');
                    ylabel(c, 'Power (mV^2)','FontSize',12);
                end
            else
                figure();
                imagesc(Spec_t,Spec_f,SpecMat');
                axis xy; axis tight;
                xlabel('Time (s)');
                ylabel('Frequency (Hz)');
                title('Multitaper Spectrogram (Fixed Window)');
                colormap('jet'); c = colorbar;
                ylabel(c, 'Power (mV^2)','FontSize', 12);
            end
        else
            error('Redundant frequency value occurred during estimation of scaled frequency vector.');
        end
    end

    %elapsed_time = toc;
    %fprintf('Complete! The process has taken %.3f seconds. \n', elapsed_time);
    %% Display Spectrogram Properties
    if verbose
        if nOscCycle == 0
            display_summary(df,dt,window_size,Fs,nw,ntp);
        else
            display_summary(dfs',dt,window_sizes',Fs,nw,ntp);
        end
    end
end

%% Appendix: In-Script Functions
% Function #1: Displays a summary of MTP sepctrogram parameters and
% resoultions
function display_summary(df, dt, window_size, Fs, nw, ntp)
    disp('*****************SUMMARY*****************');
    disp('Multitaper Spectrogram Parameters:');
    disp(['    Spectral Resolution (Hz): ' num2str(2*nw./(window_size./Fs))]);
    disp(['    Frequency Bin (FFT binwdith) (Hz): ' num2str(df)]);
    disp(['    Window Step / Time Bin: ' num2str(dt) 's']);
    disp(['    Window Length (s): ' num2str(window_size/Fs)]);
    disp(['    Time Half-Bandwidth Product: ' num2str(nw)]);
    disp(['    Number of Tapers: ' num2str(ntp)]);
    disp('*****************************************');
    % The spectral resolution corresponds to the bandwidth of the main lobe in
    % the spectral estimate, which specifies the minimum distance between peaks
    % that can be resolved [1].
end