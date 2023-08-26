function [stime, sburst, binSpec] = detect_burst_spectrogram(Spec_f, Spec_t, Spec, dt, thrA_factor, verbose)
    %% Function: 'detect_burst_spectrogram'
    % DESCRIPTION
    % Detects bursts in the frequency domain by imposing a power and duration
    % threshold on a power spectrogram.

    % USAGE
    % Full Input : detect_burst_spectrogram(Spec_f,Spec_t,Spec,dt,verbose)
    % Example    : detect_burst_spectrogram(Spec_f,Spec_t,Spec,1/Fs,true)

    % INPUT
    %    Variable     Data Type                 Description
    % 1. Spec_f       [1 x N vector]          : frequency vector of the spectrogram
    % 2. Spec_t       [1 x N vector]          : time vector of the spectrogram
    % 3. Spec         [N x N matrix]          : power spectrogram
    % 4. dt           [rational Q > 0]        : time step used for spectrogram computation
    % 5. thrA_factor  [number N]              : power threshold factor (i.e., multiplier of standard deviation)
    % 6. verbose      [logical]               : indicator for plotting a quantized spectrogram
    %                                           1) false or 0 (default)
    %                                           2) true or 1
    %                                              NOT recommended for the spectrograms with high temporal resolution (i.e., small dt)

    % NOTE: Only the bursts of frequencies lower than or equal to 1/dt Hz can be detected using this algorithm, because the spectrogram is sampled by the
    % period dt. Therefore, to detect the bursts in higher frequency bands, a spectrogram with windows moving in smaller steps should be used. However, this 
    % should not be confused with the sampling period of the actual data or the temporal resolution of a spectrogram given by the window length.

    % OUTPUT
    %    Variable     Data Type                 Description
    % 1. stime        [1 x N cell vector]     : set of time arrays of the detected bursts
    % 2. sburst       [1 x N cell vector]     : set of power arrays of the detected bursts
    % 3. binSpec      [N x N binary matrix]   : spectrogram quantized to 0 and 1 according to the thresholds

    % Written by SungJun Cho, April 7, 2021
    % Modified on October 29, 2021
    %% Set Input Parameters
    if nargin < 6
        verbose = false;
    end
    %% Apply Power Threshold
    binSpec = zeros(size(Spec)); % preallocate binary spectrogram matrix
    % [1] Quantize Spectrogram in Each Frequency Band
    for j = 1:size(Spec,2) % for each frequency in a column
        % Consider the windowing effect of the spectrogram
        % Since the beginning and end of the spectrogram are cut off (i.e.,
        % zero-padded), including those values in the threshold calculation
        % will mitigate the threshold value.
        sel_Spec = Spec(:,j)';
        sel_Spec(sel_Spec~=0) = 1;
        start = strfind([0 sel_Spec], [0 1]);
        stop = strfind([sel_Spec 0], [1 0]);
        if length(start) > 1 || length(stop) > 1
            start = start(1); stop = stop(end);
        end
        thr = mean(Spec(start:stop,j)) + thrA_factor*std(Spec(start:stop,j)); % define power threshold
        for i = 1:size(Spec,1) % for each time point in a row
            if Spec(i,j) > thr
                binSpec(i,j) = 1;
            else
                binSpec(i,j) = 0;
            end
        end
    end
    %% Extract Separate Spectral Bursts
    % [1] Preallocate Arrays of Burst Duration and Power Values
    stime = cell(length(Spec_f),size(binSpec,1));
    sburst = cell(length(Spec_f),size(binSpec,1));
    for k = 1:size(binSpec,2) % for each frequency in a column
        burstend = zeros(1,size(binSpec,1));
        temp_true = find(binSpec(:,k) == 1);
        if ~isempty(temp_true)
            for t = 1:length(temp_true)
                if t == length(temp_true)
                    break
                end
                % If indices above the threshold exceeds at least 1 cycle difference,
                % then regard them as separate bursts.
                if temp_true(t+1) - temp_true(t) > ((1/dt)/Spec_f(k))
                % Or alternatively: if Spec_t(temp_true(t+1)) - Spec_t(temp_true(t)) > 1/Spec_f(k)
                %                 : if temp_true(t+1) - temp_true(t) > Fs/Spec_f(k) [NOTE - only when dt=1/Fs]
                    burstend(1,t) = temp_true(t);
                end
            end
            % Consider the case where there is only one burst
            if sum(burstend) == 0
                burstend(1,length(temp_true)+1) = temp_true(end);
            end
            % Include the last detected burst if not already included
            if ~ismember(burstend,temp_true(end))
                burstend(1,length(temp_true)+1) = temp_true(end);
            end
            burstend = nonzeros(burstend)';
            for t = 1:length(burstend)
                if t == 1
                    sburst{k,t} = Spec(temp_true(t):burstend(t),k);
                    stime{k,t} = Spec_t(temp_true(t):burstend(t));
                elseif t > 1
                    start = find(temp_true == burstend(t-1),1,'last') + 1;
                    sburst{k,t} = Spec(temp_true(start):burstend(t),k);
                    stime{k,t} = Spec_t(temp_true(start):burstend(t));
                end
            end
        end
    end
    stime = stime(:,any(~cellfun('isempty',stime)));
    sburst = sburst(:,any(~cellfun('isempty',sburst)));
    if size(stime) ~= size(sburst)
        error('ComputationError: The size of "stime" and "sburst" should be identical.');
    end
    %% Apply Duration Threshold
    % [1] Eliminate Bursts with Duration Less Than 3 Oscillatory Cycles
    elimIdx = zeros(size(sburst));
    for k = 1:size(sburst,1) % for each frequency band
        criteria_dur = 3*((1/dt)/Spec_f(k)); % duration of 3 cycles for specified frequency
        % Or alternatively: 3*(Fs/Spec_f(k)) [Note - only when dt=1/Fs]
        for s = 1:size(sburst,2)
            if length(sburst{k,s}) < criteria_dur
                elimIdx(k,s) = s;
            end
        end
        eliminate = nonzeros(elimIdx(k,:))';
        for eIdx = eliminate
            sburst{k,eIdx} = [];
            stime{k,eIdx} = [];
        end
        catch_empty = cellfun('isempty',sburst(k,:));
        % 'isempty' is faster than @isempty since the former is performed inside the Mex
        % function while the latter uses mexCallMATLAB for each cell element
        if ~isequal(catch_empty,sort(catch_empty)) % if there is any empty cell indexed prior to the bursts
            temp_sburst = sburst(k,:);
            temp_stime = stime(k,:);
            temp_sburst = temp_sburst(~cellfun('isempty',temp_sburst));
            temp_stime = temp_stime(~cellfun('isempty',temp_stime));
            temp_sburst_upd = cell(1,size(sburst,2));
            temp_sburst_upd(1:length(temp_sburst)) = temp_sburst;
            temp_sburst_upd(length(temp_sburst)+1:end) = cell(1,size(sburst,2)-length(temp_sburst));
            temp_stime_upd = cell(1,size(stime,2));
            temp_stime_upd(1:length(temp_stime)) = temp_stime;
            temp_stime_upd(length(temp_stime)+1:end) = cell(1,size(stime,2)-length(temp_stime));
            sburst(k,:) = temp_sburst_upd;
            stime(k,:) = temp_stime_upd;
        end
        % Note: more concise programming is available when using Python or ignoring
        % MATLAB preallocation warning
    end
    sburst = sburst(:,any(~cellfun('isempty',sburst)));
    stime = stime(:,any(~cellfun('isempty',stime)));
    %% Update Binary Spectrogram Based on Bursts Defined by the Treshold
    for k = 1:size(stime,1)
        entire_stime = unique(sort([stime{k,:}]));
        temp_tIdx = ismember(Spec_t,entire_stime);
        temp_tIdx_not = ~ismember(Spec_t,entire_stime);
        binSpec(temp_tIdx,k) = 1;
        binSpec(temp_tIdx_not,k) = 0;
    end
    %% Visualize Quantized Spectrogram
    if verbose
        if numel(binSpec) > 1e6
            warning("binSpec too large to visualize. Manually adjust the code to visualize nevertheless.")
        else
            % [1] Adjust Spectrogram Dimensions
            if Spec_f(end) - Spec_f(1) < 0
                binSpec = fliplr(binSpec);
                Spec_f = flipud(Spec_f);
            end
            % [2] Plot Burst Events
            plot_raster(Spec_t, binSpec, Spec_f)
            title('Quantized Spectrogram');
        end
    end
    %% Auxiliary Function
    function plot_raster(time, data, trials)
        % Adopted from Alwin Gieselmann & Felix Schneider, Newcastle University (Feb 2020)
        % [1] Define Number of Trials
        nTrial = size(data,2);
        % [2] Visualize Raster Plot
        figure(); hold on;
        for n = 1:nTrial
            events = logical(data(:,n));
            t_events = time(events);
            t_events = repmat(t_events,3,1);
            y_events = nan(size(t_events));
            if ~isempty(y_events)
                y_events(1,:) = n-1;
                y_events(2,:) = n;
            end
            plot(t_events,y_events,'Color','k');
        end
        % [3] Adjust Axis and Figure Settings
        ticks = 1:nTrial;
        yticks(ticks-0.5); % center the ticks
        yticklabels(trials);
        xlabel('Time (s)');
        ylabel('Frequency (Hz)');
        box on;
        set(gca,'TickLength',[0 0],'FontSize',12);
        set(gcf,'Color','w');
    end
end