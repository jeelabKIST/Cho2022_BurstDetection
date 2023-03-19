%% Configure Library Paths
util_path = genpath('/Users/scho/Neuroscience_KIST/Cho2022_BurstDetection/utils');
data_path = genpath('/Users/scho/Neuroscience_KIST/Cho2022_BurstDetection/data');
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
f = 25;  % frequency of interest (central frequency of a frequency band)
listCycle = (((3:12).*(Fs/f))./Fs)*1000; % convert length of cycles to milliseconds
listNoise = -10:2:10;
nMethod = length(T_set);
%% Create Custom Colormaps
% [1] Okabe-Ito (colorblind-friendly)
okabe_ito = [[0,114,178];
             [213,94,0];
             [0,158,115];
             [204,121,167];
             [240,228,66]]./255;
okabe_ito_flip = flipud(okabe_ito);
okabe_ito_extend = cat(1,flipud(okabe_ito),[182,188,191]./255,[0,0,0]); % add grey and black
% [2] Blue2Red
blue2white = interp1(1:2,[50,115,156;255,255,255]/255,linspace(1,2,200),'linear');
white2red = interp1(1:2,[255,255,255;159,59,50]/255,linspace(1,2,800),'linear');
blue2red = [blue2white;white2red];
%% Schematics of Temporal Concurrence
% [1] Set Relevant Parameters
r = 0.25; % cosine fraction
f = 25; % frequency of interest
Fs = 512; % sampling rate
dt = 1/Fs; % sampling period
t = 1:dt:3; % time vector
x = zeros(1,length(t));
location = 1.5;
duration = round(8*(Fs/f));
% [2] Build Tukey Bursts
tstart = find(t == location);
burst = sin(2*pi*f*t(tstart:tstart+duration));
tukey_burst = burst.*tukeywin(length(burst),r)';
x(tstart:tstart+duration) = tukey_burst;
% [3] Build Square Pulses
sp = zeros(4,length(t));
idx = 1:length(t);
idx = idx(x~=0);
tstart = idx(1);
tend = idx(end);
move = round(3*(Fs/f));
shft = round(1.5*(Fs/f));
sp(1,tstart:tend) = 1.2;
sp(2,tstart+shft:tend-shft) = 1.2;
sp(3,tstart+move:tend+move) = 1.2;
sp(4,tstart-move:tend+move) = 1.2;
% [4] Set Visualization Parameters
fnt_sz = 30;
ovlp_clr = '#fcc792'; % light orange
orange = '#d15b00';   % orange
blue = '#026ba3';     % dark blue
grey = '#4d4d4d';     % dark grey
% [5] Plot Schematics
figure();
axs = cell(1,size(sp,1));
for i = 1:length(axs)
    axs{i} = subplot(1,4,i); hold on;
    % Distance of the burst
    s_vector = x~=0;
    [s_start,s_end] = hp.find_binary_idx(s_vector);
    % Distance of the detection
    d_vector = sp(i,:)~=0;
    [d_start,d_end] = hp.find_binary_idx(d_vector);
    hmin = min([t(s_start), t(d_start)]);
    hmax = max([t(s_end), t(d_end)]);
    if i >= 3
        d_end = s_end;
    end
    if i == 4
        d_start = s_start;
    end
    % Plot Overlap
    rectangle('Position',[t(d_start),-3.2,t(d_end)-t(d_start),4.7],'EdgeColor',ovlp_clr,'FaceColor',ovlp_clr);
    % Plot bursts
    plot(t,x,'LineWidth',3,'Color',grey);
    line([hmin hmax],[1.2 1.2],'Color',blue,'LineWidth',4);
    line([hmax hmax],[1.0 1.4],'Color',blue,'LineWidth',4);
    line([hmin hmin],[1.0 1.4],'Color',blue,'LineWidth',4);
    % Plot Square Pulses
    plot(t,sp(i,:)-3,'LineWidth',3,'Color',orange);
    xlim([t(1)+0.3 t(tend+(tstart-1))-0.3]);
    ylim([-5 2]);
    % Text Box
    mid_pt = hmin + (hmax-hmin)/2;
    text(mid_pt,2,'Union','Color',blue,'LineWidth',5,'FontSize',fnt_sz,'FontWeight','bold','HorizontalAlignment','center');
    mid_pt = t(d_start) + (t(d_end)-t(d_start))/2;
    text(mid_pt,-3.6,'Overlap','Color',orange,'LineWidth',5,'FontSize',fnt_sz,'FontWeight','bold','HorizontalAlignment','center');
end
% Plot Scale Bar
xlim_r = t(d_end+(d_start-1));
line([xlim_r-0.6,xlim_r-0.3],[-4.5,-4.5],'Color','k','LineWidth',3);
text(xlim_r-0.56,-4.92,'300ms','Color','k','LineWidth',5,'FontSize',fnt_sz,'FontWeight','bold');
% Figure Settings
set([axs{:}],'TickLength',[0,0],'XColor','none','YColor','none','Color','none');
set(gcf,'Color','w','Position',[6,411,1915,428]);
%% Visualize Temporal Concurrence Performances
% [1] Set Axes Parameters
vmin = 1; vmax = 5;
convert_opt = true;
round_opt = 'integer';
vmin_diff = 0;
vmax_diff = 100;
annot_fmt = '%.0f';
fig_pos = [1571,30,1230,976];
ax_pos = struct('pos_type','Position','pos_coord',[0.2545,0.1576,0.5665,0.7113]);
axis_params1 = struct('fnt_sz',52,'txt_sz',45,'fig_pos',fig_pos,'ax_pos',ax_pos,...
    'annot_fmt','%.0f','cbar_opt','off','xlbl_opt','on','ylbl_opt','on');
ax_pos = struct('pos_type','Position','pos_coord',[0.29,0.3486,0.4277,0.5428]);
axis_params2 = struct('fnt_sz',40,'txt_sz',34,'fig_pos',fig_pos,'ax_pos',ax_pos,'annot_fmt','%.0f','annot_rng',[55, 10],...
    'cbar_opt','on','cbar_loc','southoutside','cbar_lbl','\Delta Percentage','vint',20,'xlbl_opt','on','ylbl_opt','off');
% [2] Compute Heatmaps by Performances
[matT_best,matC_best,matT_poor,matC_poor] = get_best_and_poor_heatmaps(listNoise,listCycle,HEATMAP.concurrence,convert_opt,round_opt);
% [3] Plot Heatmaps (Best Algorithms)
figure();
plot_optimal_heatmap(listCycle,listNoise,matC_best,matT_best,okabe_ito_flip,vmin,vmax,axis_params1);
title({'Temporal Concurrence'; 'by Best Performance'},'Position',[298.2407,10.9731,0]);
% [4] Plot Heatmaps (Worst Algorithms)
axis_params.ylbl_opt = 'off';
figure();
plot_optimal_heatmap(listCycle,listNoise,matC_poor,matT_poor,okabe_ito_flip,vmin,vmax,axis_params1);
title({'Temporal Concurrence'; 'by Worst Performance'},'Position',[298.2407,10.9731,0]);
% [5] Plot Performance Differences
figure();
plot_difference_heatmap(listCycle,listNoise,matT_best,matT_poor,blue2red,vmin_diff,vmax_diff,axis_params2);
title({'Performance Difference'; 'in Temporal Concurrence'},'Position',[50.3803,41.5063,0.5]);
%% Create Manual Legend Box - 1
leg_fig = figure(); hold on;
set(gcf,'Color','w');
set(gca,'DefaultLineLineWidth',10.0,'Color','none');
plot(rand(1,2),rand(1,2),'Color',okabe_ito(1,:));
plot(rand(1,2),rand(1,2),'Color',okabe_ito(2,:));
plot(rand(1,2),rand(1,2),'Color',okabe_ito(3,:));
plot(rand(1,2),rand(1,2),'Color',okabe_ito(4,:));
plot(rand(1,2),rand(1,2),'Color',okabe_ito(5,:));
lgnd = legend({'BP','ENV','S-STFT','MTP','CWT'},'FontSize',34,'Color','w','NumColumns',5);
lgnd.EdgeColor = 'k';
lgnd.LineWidth = 3.5;
lgnd.ItemTokenSize = [50,8];
makeLegendToImage(leg_fig,lgnd,'line');
%% Plot Temporal Concurrence over Duration of Bursts for Low & High Noises
% [1] Set Input Data and Parameters
durIdx = find(ismember(listCycle,[160,240,360,480]));
% [2] Concurrence over Duration under Low Noise
vmin = 60; vmax = 90;
snrIdx = find(ismember(listNoise,[2,4,6]));
tcLN = zeros(length(snrIdx)*nMethod,length(listCycle));
k = 1;
for i = 1:nMethod
    for s = snrIdx
        tcLN(k,:) = T_set{i}(s,:);
        k = k + 1;
    end
end
tcLN = tcLN(:,durIdx);
tcLN_avg = mean(tcLN,1);
tcLN_std = std(tcLN,1);
pval_LN = kruskalwallis(tcLN,listCycle(durIdx),'off');
plot_noise_bar(tcLN,tcLN_avg,tcLN_std,durIdx,snrIdx,vmin,vmax,okabe_ito,true);
text(0.5,87,'2 – 6dB','FontSize',28,'HorizontalAlignment','left');
% [3] Concurrence over Duration under High Noise
vmin = 40; vmax = 80;
snrIdx = find(ismember(listNoise,[-2,-4,-6]));
tcHN = zeros(length(snrIdx)*nMethod,length(listCycle));
k = 1;
for i = 1:nMethod
    for s = snrIdx
        tcHN(k,:) = T_set{i}(s,:);
        k = k + 1;
    end
end
tcHN = tcHN(:,durIdx);
tcHN_avg = mean(tcHN,1);
tcHN_std = std(tcHN,1);
pval_HN = kruskalwallis(tcHN,listCycle(durIdx),'off');
plot_noise_bar(tcHN,tcHN_avg,tcHN_std,durIdx,snrIdx,vmin,vmax,okabe_ito,false);
text(4.7,79,'-6 – -2dB','FontSize',28,'HorizontalAlignment','right');
%% Plot Temporal Concurrence over SNR for Short & Long Burst Durations
% [1] Set Input Data and Parameters
snrIdx = find(ismember(listNoise,listNoise(1:2:end)));
% [2] Concurrence over Noise for Short Bursts
vmin = 50; vmax = 90;
durIdx = find(ismember(listCycle,[120,160,200]));
tcSD = zeros(length(durIdx)*nMethod,length(listNoise));
k = 1;
for i = 1:nMethod
    for d = durIdx
        tcSD(k,:) = T_set{i}(:,d);
        k = k + 1;
    end
end
tcSD = tcSD(:,snrIdx);
tcSD_avg = mean(tcSD,1);
tcSD_std = std(tcSD,1);
pval_SD = kruskalwallis(tcSD,listNoise(snrIdx),'off');
plot_duration_bar(tcSD,tcSD_avg,tcSD_std,durIdx,snrIdx,vmin,vmax,okabe_ito,true);
text(0.5,87,'120 – 200ms','FontSize',28,'HorizontalAlignment','left');
% [3] Concurrence over Noise for Longer Bursts
vmin = 30; vmax = 90;
durIdx = find(ismember(listCycle,[400,440,480]));
tcLD = zeros(length(durIdx)*nMethod,length(listNoise));
k = 1;
for i = 1:nMethod
    for d = durIdx
        tcLD(k,:) = T_set{i}(:,d);
        k = k + 1;
    end
end
tcLD = tcLD(:,snrIdx);
tcLD_avg = mean(tcLD,1);
tcLD_std = std(tcLD,1);
pval_LD = kruskalwallis(tcLD,listNoise(snrIdx),'off');
plot_duration_bar(tcLD,tcLD_avg,tcLD_std,durIdx,snrIdx,vmin,vmax,okabe_ito,false);
text(0.5,89,'400 – 480ms','FontSize',28,'HorizontalAlignment','left');
%% Create Manual Legend Box - 2
mkrList = {'o','d','s','*','x'};
leg_fig = figure(); hold on;
set(gcf,'Color','w');
set(gca,'Color','none');
p1 = plot(rand(1,2),rand(1,2),'Color',okabe_ito(1,:),'Marker',mkrList{1},'MarkerFaceColor',okabe_ito(1,:));
p2 = plot(rand(1,2),rand(1,2),'Color',okabe_ito(2,:),'Marker',mkrList{2},'MarkerFaceColor',okabe_ito(2,:));
p3 = plot(rand(1,2),rand(1,2),'Color',okabe_ito(3,:),'Marker',mkrList{3},'MarkerFaceColor',okabe_ito(3,:));
p4 = plot(rand(1,2),rand(1,2),'Color',okabe_ito(4,:),'Marker',mkrList{4},'MarkerFaceColor',okabe_ito(4,:));
p5 = plot(rand(1,2),rand(1,2),'Color',okabe_ito(5,:),'Marker',mkrList{5},'MarkerFaceColor',okabe_ito(5,:));
set([p1,p2,p3,p4,p5],'LineWidth',3,'MarkerSize',20);
lgnd = legend('BP','ENV','S-STFT','MTP','CWT','FontSize',34,'Color','w','NumColumns',5);
lgnd.LineWidth = 4;
lgnd.ItemTokenSize = [50,8];
makeLegendToImage(leg_fig,lgnd,'line');
%% Appendix: In-Script Functions
% Function #1: Plot the Effect of Noise Level on Concurrence across Burst Durations
function plot_noise_bar(concurrence,concurrence_avg,concurrence_std,durIdx,snrIdx,vmin,vmax,cmap,lbl_opt)
    % Set Visualization Parameters
    mkrList = {'o','d','s','*','x'};
    clrList = cmap;
    fig_pos = [2154,69,885,420];
    % Plot Bar Graph
    figure(); hold on;
    b = bar(concurrence_avg);
    x_axis = 1:length(durIdx);
    p = cell(size(x_axis));
    s = cell(size(x_axis));
    k = 1;
    for i = 1:size(concurrence,1)
        p{i} = plot(x_axis,concurrence(i,:),'Color',clrList(k,:));
        p{i}.Color(4) = 0.8;
        s{i} = scatter(x_axis,concurrence(i,:),250,clrList(k,:),'filled','Marker',mkrList{k},'LineWidth',3,'MarkerFaceColor',clrList(k,:));
        if k > 3
            s{i}.MarkerEdgeColor = clrList(k,:);
        end
        s{i}.MarkerFaceAlpha = 0.7;
        s{i}.MarkerEdgeAlpha = 0.7;
        if mod(i,length(snrIdx)) == 0
            k = k + 1;
        end
    end
    er = errorbar(x_axis,concurrence_avg,concurrence_std,concurrence_std);
    lg = legend('Location','southeastoutside'); lg.String = lg.String{1:4};
    xlim([0.2 4.7]);
    ylim([vmin vmax]);
    xticks(x_axis); yticks(vmin:10:vmax);
    % Figure Settings
    set(b,'EdgeColor','none','BaseValue',vmin+0.003,'ShowBaseLine','off','FaceColor','k','FaceAlpha',0.4);
    set([p{:}],'LineWidth',3);
    set(er,'Color','k','LineStyle','none','LineWidth',3);
    set(lg,'Visible','off');
    set(gca,'Box','off','TickDir','out','TickLength',[0.03,0.03],'LineWidth',3,'FontSize',28,...
        'FontWeight','bold', 'XTickLabelRotation',0,...
        'Position',[0.333857627118644,0.369306429156933,0.336198870056497,0.46174667943211]);
    set(gcf,'Color','w','OuterPosition',fig_pos);
    if ~lbl_opt
        xticklabels({'160','240','360','480'});
        xlabel('Duration (ms)');
    else
        xticklabels(categorical(NaN(1,length(x_axis))))
    end
end

% Function #2: Plot the Effect of Burst Duration on Concurrence across Noise Levels
function plot_duration_bar(concurrence,concurrence_avg,concurrence_std,durIdx,snrIdx,vmin,vmax,cmap,lbl_opt)
    % Set Visualization Parameters
    mkrList = {'o','d','s','*','x'};
    clrList = cmap;
    fig_pos = [2154,69,885,420];
    % Plot Bar Graph
    figure(); hold on;
    b = bar(concurrence_avg);
    x_axis = 1:length(snrIdx);
    p = cell(size(x_axis));
    s = cell(size(x_axis));
    k = 1;
    for i = 1:size(concurrence,1)
        p{i} = plot(x_axis,concurrence(i,:),'Color',clrList(k,:));
        p{i}.Color(4) = 0.8;
        s{i} = scatter(x_axis,concurrence(i,:),250,clrList(k,:),'filled','Marker',mkrList{k},'LineWidth',3,'MarkerFaceColor',clrList(k,:));
        if k > 3
            s{i}.MarkerEdgeColor = clrList(k,:);
        end
        s{i}.MarkerFaceAlpha = 0.7;
        s{i}.MarkerEdgeAlpha = 0.7;
        if mod(i,length(durIdx)) == 0
            k = k + 1;
        end
    end
    er = errorbar(x_axis,concurrence_avg,concurrence_std,concurrence_std);
    lg = legend('Location','southeastoutside'); lg.String = lg.String{1:4};
    xlim([0.2 6.7]);
    ylim([vmin vmax]);
    xticks(x_axis); yticks(vmin:20:vmax);
    % Figure Settings
    set(b,'EdgeColor','none','BaseValue',vmin+0.003,'ShowBaseLine','off','FaceColor','k','FaceAlpha',0.4);
    set([p{:}],'LineWidth',3);
    set(er,'Color','k','LineStyle','none','LineWidth',3);
    set(lg,'Visible','off');
    set(gca,'Box','off','TickDir','out','TickLength',[0.03,0.03],'LineWidth',3,'FontSize',28,...
        'FontWeight','bold','XTickLabelRotation',0,... 
        'Position',[0.333857627118644,0.369306429156933,0.336198870056497,0.46174667943211]);
    set(gcf,'Color','w','OuterPosition',fig_pos);
    if ~lbl_opt
        xticklabels({'-10','-6','-2','2','6','10'});
        xlabel('SNR (dB)');
    else
        xticklabels(categorical(NaN(1,length(x_axis))));
    end
end
