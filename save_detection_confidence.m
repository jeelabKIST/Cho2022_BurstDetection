%% Configure Library Paths
util_path = genpath('/Users/scho/Neuroscience_KIST/Cho2022_BurstDetection/utils');
data_path = genpath('/Users/scho/Neuroscience_KIST/Cho2022_BurstDetection/data');
addpath(util_path);
addpath(data_path);
%% Store Detection Confidence Heatmaps
base_path = '/Users/scho/Neuroscience_KIST/Cho2022_BurstDetection/data/simulation_data/HM_*.mat';
band_names = {'theta','beta','gamma'};
for i = 1:length(band_names)
    file_path = replace(base_path,'*',band_names{i});
    store_detection_confidence(file_path);
end
