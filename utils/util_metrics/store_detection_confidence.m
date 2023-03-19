function [detection_confidence] = store_detection_confidence(filepath, save_opt)
    %% Function 'store_detection_confidence'
    % DESCRIPTION
    % Computes detection confidence heatmaps for each algorithm using its 
    % F1-score and temporal concurrence values and stores the results as a
    % '.mat' file
    
    % INPUT
    %    Variable                Data Type         Description
    % 1. filepath                [string]        : path to a file containing heatmaps of F1-scores
    %                                              and temporal concurrence
    % 2. save_opt                [boolean]       : whether to save the detection confidence heatmaps or not
    
    % OUTPUT
    %    Variable                Data Type         Description
    % 1. detection_confidence    [1 x N cell]    : set of heatmaps containing detection confidence values 
    %                                              for each algorithm
    
    % Written by SungJun Cho, February 02, 2023
    % Last Modified on February 02, 2023
    %% Set Parameters
    if nargin < 2
        save_opt = true;
    end
    method_names = {'bp', 'ev', 'stp', 'mtp', 'cwt'};
    nMethod = length(method_names);
    %% Load Data
    F1_SCORE = load(filepath).HEATMAP.f1_score;
    TEMPORAL_CONCURRENCE = load(filepath).HEATMAP.concurrence;
    %% Store Detection Confidence Heatmaps
    detection_confidence = cell(1,nMethod);
    for n = 1:nMethod
        F_mat = F1_SCORE.(method_names{n});
        T_mat = TEMPORAL_CONCURRENCE.(method_names{n});
        DC = (F_mat.*exp(T_mat))./exp(1);
        detection_confidence{n} = DC;
    end
    if save_opt
        savepath = replace(filepath,'HM','DC');
        save(savepath,'detection_confidence');
        savename = split(savepath,'/');
        fprintf(['Progress: File "' savename{end} '" saved. \n']);
    end
end