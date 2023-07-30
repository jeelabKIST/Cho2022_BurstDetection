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
%% Set Visualization Parameters
% [1] Create Custom Colormaps
% (A) Okabe-Ito (colorblind-friendly)
okabe_ito = [[0,114,178];
             [213,94,0];
             [0,158,115];
             [204,121,167];
             [240,228,66]]./255;
okabe_ito_flip = flipud(okabe_ito);
okabe_ito_extend = cat(1,flipud(okabe_ito),[182,188,191]./255,[0,0,0]);
% (B) Blue2Red
blue2white = interp1(1:2,[50,115,156;255,255,255]/255,linspace(1,2,200),'linear');
white2red = interp1(1:2,[255,255,255;159,59,50]/255,linspace(1,2,800),'linear');
blue2red = [blue2white;white2red];
% [2] Set Axes Parameters
vmin = 1;
vmin_diff = 0;
vmax_diff = 100;
fig_pos = [1571,30,1286,1001];
ax_pos = struct('pos_type','InnerPosition','pos_coord',[0.2054,0.1201,0.6284,0.7619]);
axis_params1 = struct('fnt_sz',42,'txt_sz',36,'fig_pos',fig_pos,'ax_pos',ax_pos,...
    'annot_fmt','%.0f','cbar_opt','off','xlbl_opt','off','ylbl_opt','on');
axis_params2 = struct('fnt_sz',42,'txt_sz',36,'fig_pos',fig_pos,'ax_pos',ax_pos,'annot_fmt','%.0f','annot_rng',[55, 10],...
    'cbar_opt','on','cbar_loc','eastoutside','cbar_lbl','\Delta Percentage','vint',20,'xlbl_opt','off','ylbl_opt','on');
%% Visualize Precision Heatmaps
vmax = 7;
convert_opt = true;
round_opt = 'integer';
% [1] Compute Heatmaps by Performances
[matP_best,matC_best,matP_poor,matC_poor] = get_best_and_poor_heatmaps(listNoise,listCycle,HEATMAP.precision,convert_opt,round_opt);
% [2] Plot Heatmaps (Best Methods)
figure();
plot_optimal_heatmap(listCycle,listNoise,matC_best,matP_best,okabe_ito_extend,vmin,vmax,axis_params1);
title({'Precision (%) by'; 'Best Performance'});
% [3] Plot Heatmaps (Poor Methods)
figure();
plot_optimal_heatmap(listCycle,listNoise,matC_poor,matP_poor,okabe_ito_extend,vmin,vmax,axis_params1);
title({'Precision (%) by'; 'Worst Performance'});
% [4] Plot Difference in Performances
figure();
plot_difference_heatmap(listCycle,listNoise,matP_best,matP_poor,blue2red,vmin_diff,vmax_diff,axis_params2);
title({'Performance';'Difference in Precision'},'Position',[-23.8749,100.7569,0.5]);
%% Visualize Recall Heatmaps
vmax = 7;
convert_opt = true;
round_opt = 'integer';
% [1] Compute Heatmaps by Performances
[matR_best,matC_best,matR_poor,matC_poor] = get_best_and_poor_heatmaps(listNoise,listCycle,HEATMAP.recall,convert_opt,round_opt);
% [2] Plot Heatmaps (Best Methods)
figure();
plot_optimal_heatmap(listCycle,listNoise,matC_best,matR_best,okabe_ito_extend,vmin,vmax,axis_params1);
title({'Recall (%) by'; 'Best Performance'});
% [3] Plot Heatmaps (Poor Methods)
figure();
plot_optimal_heatmap(listCycle,listNoise,matC_poor,matR_poor,okabe_ito_extend,vmin,vmax,axis_params1);
title({'Recall (%) by'; 'Worst Performance'});
% [4] Plot Difference in Performances
figure();
plot_difference_heatmap(listCycle,listNoise,matR_best,matR_poor,blue2red,vmin_diff,vmax_diff,axis_params2);
title({'Performance';'Difference in Recall'},'Position',[-23.8749,100.7569,0.5]);
%% Visualize F1-Score Heatmaps
vmax = 7;
convert_opt = true;
round_opt = 'integer';
% [1] Compute Heatmaps by Performances
[matF_best,matC_best,matF_poor,matC_poor] = get_best_and_poor_heatmaps(listNoise,listCycle,HEATMAP.f1_score,convert_opt,round_opt);
% [2] Plot Heatmaps (Best Methods)
figure();
plot_optimal_heatmap(listCycle,listNoise,matC_best,matF_best,okabe_ito_extend,vmin,vmax,axis_params1);
title({'F1 Score (%) by'; 'Best Performance'});
% [3] Plot Heatmaps (Poor Methods)
figure();
plot_optimal_heatmap(listCycle,listNoise,matC_poor,matF_poor,okabe_ito_extend,vmin,vmax,axis_params1);
title({'F1 Score (%) by'; 'Worst Performance'});
% [4] Plot Difference in Performances
figure();
plot_difference_heatmap(listCycle,listNoise,matF_best,matF_poor,blue2red,vmin_diff,vmax_diff,axis_params2);
title({'Performance';'Difference in F1 Score'},'Position',[-23.8749,100.7569,0.5]);
%% Create Manual Legend Box
% [1] Legends for Individual Algorithms
leg_fig = figure(); hold on;
set(gcf,'Color','w');
set(gca,'DefaultLineLineWidth',10.0,'Color','none');
plot(rand(1,2),rand(1,2),'Color',okabe_ito(1,:));
plot(rand(1,2),rand(1,2),'Color',okabe_ito(2,:));
plot(rand(1,2),rand(1,2),'Color',okabe_ito(3,:));
plot(rand(1,2),rand(1,2),'Color',okabe_ito(4,:));
plot(rand(1,2),rand(1,2),'Color',okabe_ito(5,:));
lgnd = legend({'BP','ENV','S-STFT','MTP','CWT'},'FontSize',34,'Color','w','NumColumns',5);
lgnd.EdgeColor = 'w';
lgnd.ItemTokenSize = [50,8];
makeLegendToImage(leg_fig,lgnd,'line');
% [2] Legends for Multiple Algorithms
leg_fig = figure(); hold on;
set(gcf,'Color','w');
set(gca,'DefaultLineLineWidth',10.0,'Color','none');
plot(rand(1,2),rand(1,2),'Color',okabe_ito_extend(end-1,:));
plot(rand(1,2),rand(1,2),'Color',okabe_ito_extend(end,:));
lgnd = legend({'Equal for Multiple Algorithms','Equal for All Algorithms'},'FontSize',34,'Color','w','NumColumns',2);
lgnd.EdgeColor = 'w';
lgnd.ItemTokenSize = [50,8];
makeLegendToImage(leg_fig,lgnd,'line');
