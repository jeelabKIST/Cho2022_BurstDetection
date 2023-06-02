function [bscoreMat, auc] = compute_ecdf_auc(trialID, metric_name, scoreMat, btimeMat, dtMat, tvec, signal, Fs, lo_f, hi_f, verbose)
    %% Function 'compute_ecdf_auc'
    % DESCRIPTION
    % Plots empirical cumulative distribution functions (eCDFs) and
    % calculates their AUC
    
    % INPUT
    %    Variable       Data Type           Description
    % 1. trialID        [integer Z > 0]   : number index of a trial
    % 2. metric_name    [string]          : name of the chosen performance metric
    % 3. scoreMat       [1 x N cell]      : set of heatmaps of a chosen performance metric for each algorithm
    % 4. btimeMat       [1 x N cell]      : set of time arrays of detected bursts for each algorithm
    % 5. dtMat          [1 x N cell]      : set of a sampling rate or sliding window time step for each algorithm
    % 6. tvec           [1 x N double]    : time array for the input signal
    % 7. signal         [1 x N double]    : input signal
    % 8. Fs             [number N]        : sampling frequency
    % 9. lo_f           [number N]        : low cutoff frequency of the input signal
    % 10. hi_f          [number N]        : high cutoff frequency of the input signal
    % 11. verbose       [boolean]         : whether to plot burst occurrences and eCDFs
    
    % OUTPUT
    %    Variable       Data Type           Description
    % 1. bscoreMat     [1 x N cell]       : performance metric values mapped from burst detections of each algorithm
    % 2. auc           [1 x N double]     : AUC values for each algorithm
    
    % Written by SungJun Cho, August 15, 2021
    % Last Modified on April 30, 2023
    %% Set Input Parameters
    if nargin < 11
        verbose = true;
    end
    nMethod = length(btimeMat); % number of algorithms
    okabe_ito = [[0,114,178];
                 [213,94,0];
                 [0,158,115];
                 [204,121,167];
                 [240,228,66]]./255; % color map
    %% Visualize Burst Time Points
    if verbose
        plot_burst_raster_plot(btimeMat);
    end
    %% Compute Burst Cycle Lengths
    foi = (lo_f + hi_f) / 2;
    wrapper = @(x,y) sec2cyc(x,foi,y);
    blenMat = cellfun(wrapper, btimeMat, dtMat, 'UniformOutput', false); % length (in cycles) of all bursts
    % Since bursts longer than 12 cycles or shorter than 3 cycles can be 
    % detected, we bound maximum length of bursts to 12 cycles and its 
    % minimum to 3 cycles in order to compare them with decision matrices.
    for n = 1:nMethod
        for i = 1:length(blenMat{n})
            if blenMat{n}(i) > 12
                blenMat{n}(i) = 12;
            end
            if blenMat{n}(i) < 3
                blenMat{n}(i) = 3;
            end
        end
    end
    %% Compute Burst Signal-to-Noise Ratio
    snrMat = get_snr(btimeMat,tvec,signal,Fs,lo_f,hi_f);
    % Since bursts with SNR higher than 10 dB or lower than -10 dB can be 
    % detected, we bound maximum SNR of bursts to 10 dB and minimum to -10 
    % dB in order to compare them with decision matrices.
    for n = 1:nMethod
        for i = 1:length(snrMat{n})
            if snrMat{n}(i) > 10
                snrMat{n}(i) = 10;
            end
            if snrMat{n}(i) < -10
                snrMat{n}(i) = -10;
            end
        end
    end
    %% Get Metric Scores for Each Burst
    cycList = 3:12;
    snrList = -10:2:10;
    bscoreMat = cell(size(btimeMat)); % metric scores for each burst
    for n = 1:nMethod
        score = scoreMat{n};
        blen = blenMat{n};
        bsnr = snrMat{n};
        if length(blen) < 2 % if no bursts or only one burst is detected
            bscoreMat{n} = NaN;
            continue;
        end
        if length(blen) == length(bsnr)
            bscore = zeros(size(blen));
            for i = 1:length(bscore)
                cycIdx = find(blen(i) == cycList);
                if mod(bsnr(i),2) == 0
                    snrIdx = find(bsnr(i) == snrList);
                    bscore(i) = score(snrIdx,cycIdx);
                else
                    snrIdx_lo = find(bsnr(i)-1 == snrList);
                    snrIdx_hi = find(bsnr(i)+1 == snrList);
                    bscore(i) = mean([score(snrIdx_lo,cycIdx),score(snrIdx_hi,cycIdx)]);
                end
            end
        end
        bscoreMat{n} = bscore;
    end
    %% Compute and Plot Empirical Cumulative Distribution Function and Its AUC
    [cdf_x, cdf_y] = deal(cell(1,nMethod));
    auc = zeros(1,nMethod);
    figure(); hold on;
    xmin = 0;
    xmax = 0;
    h_set = zeros(1,nMethod);
    for n = 1:nMethod
        if isnan(bscoreMat{n})
            auc(n) = NaN;
            continue;
        end
        [h, ~] = cdfplot(sort(bscoreMat{n}));
        h.Color = okabe_ito(n,:);
        h.Color(4) = 0.9;
        h.LineWidth = 4.5;
        h.XData(1) = 0;
        h.XData(end) = max(cat(2,bscoreMat{:}),[],'all');
        if xmax < h.XData(end)
            xmax = h.XData(end);
        end
        cdf_x{n} = h.XData;
        cdf_y{n} = h.YData;
        auc(n) = trapz(cdf_x{n},cdf_y{n});
        h_set(n) = h;
    end
    if sum(isnan([bscoreMat{:}])) == 0
        baseline = plot(linspace(xmin, xmax, 100),linspace(0, 1, 100),'Color','k','LineWidth',3.5,'LineStyle','--');
        baseline.Color(4) = 0.5;
        uistack(baseline,'bottom');
        h.Parent.XLabel.String = metric_name;
        h.Parent.YLabel.String = 'Cumulative Probability';
        title({['Example: Trial #' num2str(trialID)], [metric_name ' eCDF']});
        yticks(linspace(0,1,6));
        leg_labels = {'BP','ENV','S-STFT','MTP','CWT'};
        for n = 1:nMethod
            leg_labels{n} = [leg_labels{n} ' (AUC = ' num2str(auc(n),'%.3f') ')'];
        end
        leg = legend(h_set,leg_labels,'Location','best');
        leg.ItemTokenSize = [15,18];
        grid off;
        set(leg,'FontSize',24,'FontWeight','bold');
        set(gca,'FontSize',35,'FontWeight','bold','LineWidth',5,'YLim',[-0.01,1.01],'TickDir','Out','TickLength',[0.015, 0.025], ...
            'XMinorTick','off','YMinorTick','off','Position',[0.191107069608701,0.151302987340527,0.718640222810071,0.712761615620442]);
        set(gcf,'Color','w','Position',[703,199,831,743]);
    else
        warning('Cannot plot ECDF. At least one algorithm did not detect any bursts.');
    end
    if ~verbose
        close(gcf);
    end
end

%% Appendix: In-Script Functions
% Function #1: Plots time intervals of the detected bursts
function plot_burst_raster_plot(btimeMat)
    nMethod = length(btimeMat);
    figure(); hold on;
    for n = 1:nMethod
        btime = btimeMat{n};
        for b = 1:length(btime)
            scatter(btime{b}, ones(size(btime{b}))*n, 'k', 'filled');
        end
    end
    ylim([0 nMethod+1]);
    yticks(1:nMethod); yticklabels({'BP','ENV','STP','MTP','CWT'});
    xlabel('Time (s)');
    ylabel('Algorithm');
    title('Raster Plots');
    set(gca,'FontSize',15,'TickDir','out','LineWidth',2);
    set(gcf,'Color','w');
end

% Function #2: Converts burst lengths from samples to cycles
function [num_cycles] = sec2cyc(btime,foi,dt)
    num_cycles = cellfun(@(b) round(length(b)*(foi*dt)), btime);
end

% Function #3: Compute SNR (dB) of each burst
function [snrMat] = get_snr(btimeMat,tvec,signal,Fs,lo_f,hi_f)
    % [1] Get bandpass-filtered signals with and without bursts
    filt_sig = eegfilt(signal,Fs,lo_f,hi_f);
    % [2] Measure SNR for each burst
    snrMat = cell(size(btimeMat));
    for i = 1:length(btimeMat)
        btime = btimeMat{i};
        bsig = cell(size(btime));
        qnoise = ones(1,length(tvec));
        for n = 1:length(btime)
            % (A) For time domain methods
            if i < 3
                bIdx = ismember(tvec,btime{n});
                bsig{n} = filt_sig(bIdx);
                qnoise(bIdx) = 0;
            % (B) For time-frequency domain methods
            else
                binterval = [btime{n}(1) btime{n}(end)];
                bIdx = find(tvec >= binterval(1) & tvec <= binterval(end));
                bsig{n} = filt_sig(bIdx);
                qnoise(bIdx) = 0;
            end
        end
        amp_sig = cell2mat(cellfun(@(x) rms(x), bsig, 'UniformOutput', false));
        amp_noise = rms(filt_sig(logical(qnoise)));
        snrMat{i} = round(10*log10((amp_sig/amp_noise).^2));
    end
end
