%% Configure Library Paths
util_path = genpath('/Users/jeelab/Desktop/Cho2022_BurstDetection/utils');
data_path = genpath('/Users/jeelab/Desktop/Cho2022_BurstDetection/data');
addpath(util_path);
addpath(data_path);
%% Load Data
% [1] Organize Data
HEATMAP = load('HM_beta.mat').HEATMAP;
hp = helper;
[T_bp,T_ev,T_stp,T_mtp,T_cwt] = hp.unpack_data(HEATMAP.concurrence);
[T_bp,T_ev,T_stp,T_mtp,T_cwt] = hp.convert_to_percent(T_bp,T_ev,T_stp,T_mtp,T_cwt);
T_set = {T_bp,T_ev,T_stp,T_mtp,T_cwt};
% [2] Set Parameters
Fs = 512; % sampling rate
f1 = 25;  % frequency of interest
listCycle = (((3:12).*(Fs/f1))./Fs)*1000; % convert length of cycles to milliseconds
listNoise = -10:2:10;
nMethod = length(T_set);
%% Create Custom Colormap (Red2Blue)
red2white = interp1(1:2,[159,59,50;255,255,255]/255,linspace(1,2,650),'linear');
white2blue = interp1(1:2,[255,255,255;50,115,156]/255,linspace(1,2,350),'linear');
red2blue = [red2white;white2blue];
%% Plot Heatmaps For Each Individual Algorithm
% [1] Set Visualization Parameters
vmin = 40; vmax = 90;
fig_pos = [352,1,1230,976];
ax_pos = struct('pos_type','Position','pos_coord',[0.1967,0.1771,0.5769,0.7245]);
axis_params1 = struct('fnt_sz',52,'txt_sz',45,'fig_pos',fig_pos,'ax_pos',ax_pos,'annot_fmt','%.0f',...
    'annot_rng',[80,60],'vint',10,'cbar_opt','off','cbar_lbl','Temporal Concurrence (%)','ylbl_opt','on');
axis_params2 = struct('fnt_sz',52,'txt_sz',45,'fig_pos',fig_pos,'ax_pos',ax_pos,'annot_fmt','%.0f',...
    'annot_rng',[80,60],'vint',10,'cbar_opt','on','cbar_lbl','Temporal Concurrence (%)','ylbl_opt','off');
axis_params3 = axis_params1; axis_params3.ylbl_opt = 'off';
% [2] Construct Heatmaps
figure(); plot_heatmap(listCycle,listNoise,T_bp,red2blue,vmin,vmax,'Bandpass Filtering',axis_params1);
figure(); plot_heatmap(listCycle,listNoise,T_ev,red2blue,vmin,vmax,'Envelope-Based',axis_params3);
figure(); plot_heatmap(listCycle,listNoise,T_stp,red2blue,vmin,vmax,'Single-Tapered',axis_params2);
figure(); plot_heatmap(listCycle,listNoise,T_mtp,red2blue,vmin,vmax,'Multitaper',axis_params1);
figure(); plot_heatmap(listCycle,listNoise,T_cwt,red2blue,vmin,vmax,'Wavelet',axis_params2);