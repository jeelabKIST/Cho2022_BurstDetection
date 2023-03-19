function [ANNOT] = read_annotation_file(sheet_name, nTrial, tvec, util_path, data_path)
    %% Function 'read_annotation_file'
    % DESCRIPTION
    % Reads a .xlsx file containing burst annotations and outputs onsets
    % and offsets of the bursts detected by every evaluators
    
    % INPUT
    %    Variable       Data Type           Description
    % 1. sheet_name     [string]          : name of the excel sheet to read
    % 2. nTrial         [integer Z > 0]   : number of total trials available from the data
    % 3. tvec           [1 x N double]    : time vector array of the data
    % 4. util_path      [string]          : path to relevant functions
    % 5. data_path      [string]          : path to burst annotations
    
    % OUTPUT
    %    Variable       Data Type           Description
    % 1. ANNOT          [struct]          : burst onsets and offsets of each trial
    
    % Written by SungJun Cho, February 12, 2023
    % Last Modified on February 17, 2023
    %% Configure Library Paths
    addpath(util_path);
    addpath(data_path);
    %% Load Annotation File
    data = readmatrix('annot/annotations.xlsx','Sheet',sheet_name);
    nResearcher = length(unique(data(:,1)));
    %% Separate Data by Researcher
    data_annot = cell(1,nResearcher);
    idx_annot = cell(1,nResearcher);
    for n = 1:nResearcher
        data_annot{n} = data(data(:,1)==n,[2,6,7]);
        idx_annot{n} = unique(data_annot{n}(:,1));
    end
    if idx_annot{1} ~= mean([idx_annot{:}],2)
        error("The number of trials annotated by each researcher does not match.");
    end
    idx_annot = mean([idx_annot{:}],2);
    %% Get Bursts Annotated by All Researchers
    final_data_annot = cell(1,length(idx_annot));
    for i = 1:length(idx_annot)
        trial_annot = cellfun(@(x) x(x(:,1)==idx_annot(i),:),data_annot,'UniformOutput',false); % annotations of all researchers in each trial
        % Get Binarized Burst Detection Time Courses
        binary_dtc = cell(1,nResearcher);
        for n = 1:nResearcher
            binary_dtc{n} = get_detection_time_course(tvec,trial_annot{n}(:,2:end));
        end
        binary_dtc = sum(cat(1,binary_dtc{:}),1);
        % Get Center Time of Each Detection
        threshold_count = nResearcher;
        center_time = get_center_time(tvec,binary_dtc,threshold_count);
        % Get Bursts that are Detected by Every Evaluators
        trial_annot_selected = cell(size(trial_annot));
        for n = 1:nResearcher
            trial_annot_selected{n} = select_burst_detections(trial_annot{n},center_time);
        end
        final_trial_annot = mean(cat(3,trial_annot_selected{:}),3);
        final_data_annot{i} = final_trial_annot;
    end
    final_data_annot = cat(1,final_data_annot{:});
    %% Organize Data by Trials
    ANNOT = struct();
    for n = 1:nTrial
        if ismember(n,final_data_annot(:,1)) % if there exists human annotations for a trial
            data_name = ['trial' num2str(n)];
            data_idx = final_data_annot(:,1) == n;
            ANNOT.(data_name) = final_data_annot(data_idx,end-1:end);
        else
            continue;
        end
    end
end

%% Appendix: In-Script Functions
% Function #1: Extracts a binary detection time course
function [dtc] = get_detection_time_course(tvec, data)
    dtc = zeros(size(tvec));
    nBurst = size(data,1);
    for n = 1:nBurst
        onset = data(n,1);
        offset = data(n,2);
        onidx = find(tvec >= onset,1,'First');
        offidx = find(tvec <= offset,1,'Last');
        dtc(onidx:offidx) = 1;
    end   
end

% Function #2: Extracts central time points of the bursts detected by every
% evaluators
function [center_time] = get_center_time(tvec, dtc, threshold_count)
    dtc(dtc < threshold_count) = 0;
    dtc(dtc >= threshold_count) = 1;
    hp = helper;
    [idx_start, idx_end] = hp.find_binary_idx(dtc);
    center_time = zeros(1,length(idx_start));
    for i = 1:length(center_time)
        center_time(i) = mean([tvec(idx_start(i)), tvec(idx_end(i))]);
    end
end

% Function #3: Select bursts that occurred during specific time points
function [selected_data] = select_burst_detections(data, center_time)
    nBurst = size(data,1);
    keepidx = zeros(1,nBurst);
    for ct = center_time
        for n = 1:nBurst
            onset = data(n,end-1);
            offset = data(n,end);
            if (ct >= onset) && (ct <= offset)
                keepidx(n) = 1;
            end
        end
    end
    selected_data = data(logical(keepidx),:);
end