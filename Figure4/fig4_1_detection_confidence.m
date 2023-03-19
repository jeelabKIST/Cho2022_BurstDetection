%% Configure Library Path
util_path = genpath('/Users/scho/Neuroscience_KIST/Cho2022_BurstDetection/utils');
data_path = genpath('/Users/scho/Neuroscience_KIST/Cho2022_BurstDetection/data');
addpath(util_path);
addpath(data_path);
%% Load Data
HEATMAP = load('HM_beta.mat').HEATMAP;
Fs = 512; f = 25;
listCycle = (((3:12).*(Fs/f))./Fs)*1000; % convert length of cycles to milliseconds
listNoise = -10:2:10;
nMethod = 5;
%% Compute Detection Confidence
dc = compute_detection_confidence(HEATMAP);
dc = cellfun(@(mat) mat.*100, dc, 'UniformOutput', false); % convert to percentage scale
%% Set Visualization Parameters
% [1] Set Axes Parameters
vmin = 0;
vmax = 100;
f1_cutoff = 0.75;
tc_cutoff = 0.75;
dc_cutoff = (f1_cutoff * exp(tc_cutoff) / exp(1)) * 100;
% [2] Create Custom Colormap: Red2Blue
r2w_len = round(((dc_cutoff-vmin)/(vmax-vmin))*1000);
w2b_len = round(((vmax-dc_cutoff)/(vmax-vmin))*1000);
red2white = interp1(1:2,[159,59,50;255,255,255]/255,linspace(1,2,r2w_len),'linear');
white2blue = interp1(1:2,[255,255,255;30,106,156]/255,linspace(1,2,w2b_len),'linear');
red2blue = [red2white;white2blue];
%% Plot Detection Confidence Heatmaps For Each Algorithm
fig_pos = [1571,30,1230,976];
ax_pos = struct('pos_type','Position','pos_coord',[0.1967,0.1771,0.5769,0.7245]);
titles = {'Bandpass Filtering', 'Envelope-Based', 'Single-Tapered', 'Multitaper', 'Wavelet'};
cbar_opts = {'off','off','on','off','on'};
ylbl_opts = {'on','off','off','on','off'};
for n = 1:nMethod
    axis_params = struct('fnt_sz',52,'txt_sz',45,'fig_pos',fig_pos,'ax_pos',ax_pos,'annot_fmt','%.0f',...
        'annot_rng',[68,35],'vint',20,'cbar_opt',cbar_opts{n},'cbar_lbl','Detection Confidence (%)','ylbl_opt',ylbl_opts{n});
    figure();
    plot_heatmap(listCycle,listNoise,dc{n},red2blue,vmin,vmax,titles{n},axis_params);
end
%% Appendix: In-Script Functions
% Function #1: Compute Detection Confidence
function [dc] = compute_detection_confidence(DATA)
    hp = helper;
    [f1_bp,f1_ev,f1_stp,f1_mtp,f1_cwt] = hp.unpack_data(DATA.f1_score);
    [tc_bp,tc_ev,tc_stp,tc_mtp,tc_cwt] = hp.unpack_data(DATA.concurrence);
    f1 = {f1_bp,f1_ev,f1_stp,f1_mtp,f1_cwt};
    tc = {tc_bp,tc_ev,tc_stp,tc_mtp,tc_cwt};
    dc = cell(1,5);
    for i = 1:length(dc)
        dc{i} = (f1{i} .* exp(tc{i})) / exp(1);
    end
end