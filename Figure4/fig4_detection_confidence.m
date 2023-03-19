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
vmin_diff = 0;
vmax_diff = vmax - vmin;
f1_cutoff = 0.75;
tc_cutoff = 0.75;
dc_cutoff = f1_cutoff * exp(tc_cutoff) / exp(1) * 100;
% [2] Create Custom Colormaps
% (A) Okabe-Ito (colorblind-friendly)
okabe_ito = [[0,114,178];
             [213,94,0];
             [0,158,115];
             [204,121,167];
             [240,228,66]]./255;
okabe_ito_flip = flipud(okabe_ito);
okabe_ito_extend = cat(1,flipud(okabe_ito),[182,188,191]./255,[0,0,0]);
% (B) Blue2Red
b2w_len = 200;
w2r_len = 800;
blue2white = interp1(1:2,[50,115,156;255,255,255]/255,linspace(1,2,b2w_len),'linear');
white2red = interp1(1:2,[255,255,255;159,59,50]/255,linspace(1,2,w2r_len),'linear');
blue2red = [blue2white;white2red];
%% Visualize Detection Confidence Schematics
% [1] Set Visualization Parameters
f1_min = 0; f1_max = 1;
[tc_min,tc_max] = find_metric_minmax(HEATMAP.concurrence);
cmap = parula; cmap = cmap(1:235,:);
fig_pos = [1571,30,1230,976];
ax_pos = [0.2138,0.1537,0.5455,0.6926];
% [2] Plot Schematics
figure();
x_axis = linspace(f1_min,f1_max,11);
y_axis = linspace(tc_min,tc_max,11);
[X,Y] = meshgrid(x_axis,y_axis);
Z = (X.*exp(Y))/exp(1);
colormap(cmap);
imagesc(x_axis,y_axis,Z);
xticks(x_axis);
yticks(y_axis);
ln{1} = xline(f1_cutoff);
ln{2} = yline(tc_cutoff);
axis xy;
val = cellfun(@(x) num2str(x,'%.2f'),num2cell(Z),'UniformOutput',false); % convert values to strings
for i = 1:length(x_axis)
    for j = 1:length(y_axis)
        text(X(i,j),Y(i,j),val(i,j),'HorizontalAlignment','Center','FontSize',27,'Color','w')
    end
end
pause(0.8);
xlabel('F1-Score');
ylabel('Temporal Concurrence');
title({'Schematics of';'Detection Confidence Scores'});
% (A) Axis Settings
ax = gca;
ax.XTickLabelRotation = 0;
ax.XRuler.Axle.LineStyle = 'none';
ax.YRuler.Axle.LineStyle = 'none';
ax.XAxis.MajorTickChild.LineWidth = 4;
ax.YAxis.MajorTickChild.LineWidth = 4;
% (B) Figure Settings
set([ln{:}],'Color','#f66e45','LineWidth',7,'LineStyle','--','Alpha',1.0);
set(ax,'TickDir','out','Box','off','FontSize',34,'FontWeight','bold','LineWidth',2,'Color','none','Position',ax_pos);
set(ax,'XTickLabel',num2str(get(ax,'xtick')','%.1f'),'YTickLabel',num2str(get(ax,'ytick')','%.2f'));
set(gcf,'Color','w','Position',fig_pos);
%% Plot Manual Horizontal Colorbar
lbl_name = 'Detection Confidence';
dc_min = min(Z,[],'all');
dc_max = max(Z,[],'all');
figure();
plot_manual_colorbar(x_axis,y_axis,dc_cutoff,cmap,dc_min,dc_max,lbl_name);
%% Plot Heatmaps by Best and Poor Performances
% [1] Set Visualization Parameters
vmin = 1; vmax = 7; 
convert_opt = false;
round_opt = 'decimal';
fig_pos = [1571,30,1230,976];
ax_pos = struct('pos_type','Position','pos_coord',[0.2138,0.1537,0.5455,0.6926]);
axis_params = struct('fnt_sz',52,'txt_sz',42,'fig_pos',fig_pos,'ax_pos',ax_pos,...
    'annot_fmt','%.0f','cbar_opt','off','xlbl_opt','on','ylbl_opt','on');
% [2] Compute Heatmaps by Performances
[matD_best,matC_best,matD_poor,matC_poor] = get_best_and_poor_heatmaps(listNoise,listCycle,dc,convert_opt,round_opt);
% [3] Plot Heatmaps (Best Methods)
figure();
plot_optimal_heatmap(listCycle,listNoise,matC_best,matD_best,okabe_ito_extend,vmin,vmax,axis_params);
title({'Detection Confidence by'; 'Best Performance'});
% [4] Plot Heatmaps (Poor Methods)
axis_params.ylbl_opt = 'off';
figure();
plot_optimal_heatmap(listCycle,listNoise,matC_poor,matD_poor,okabe_ito_extend,vmin,vmax,axis_params);
title({'Detection Confidence by'; 'Worst Performance'});
% [5] Plot Difference in Performances
axis_params.vint = 20;
axis_params.annot_rng = [31, 10];
axis_params.cbar_opt = 'on';
axis_params.cbar_loc = 'eastoutside';
axis_params.cbar_lbl = '\Delta Confidence (%)';
figure();
plot_difference_heatmap(listCycle,listNoise,matD_best,matD_poor,blue2red,vmin_diff,vmax_diff,axis_params);
title({'Performance Difference';'in Detection Confidence'},'Position',[-22.06244283914566,102.1153846153846,0.5]);
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
%% Plot Decision Matrix
% [1] Compute Decision Matrix
[decMat,clrMat,list_uqid,list_algo] = extract_decision_matrix(matD_best,dc);
% [2] Set Visualization Parameters
okabe_ito_expand = cat(1,okabe_ito,[182,188,191]./255,[0,0,0]);
fig_pos = [1571,30,1230,976];
ax_pos = [0.2138,0.1537,0.5455,0.6926];
% [3] Visualize Decision Matrix
figure();
imagesc(listCycle,listNoise,clrMat);
xticks(listCycle);
yticks(listNoise);
xvar_lbl = categorical(cellfun(@(x) mod(x,listCycle(1)) == 0, num2cell(listCycle)) .* listCycle);
xvar_lbl(xvar_lbl == categorical(0)) = ' ';
xticklabels(xvar_lbl);
colormap(okabe_ito_expand); clim([1,7]);
axis xy;
x = repmat(listCycle,length(listNoise),1);
y = repmat(listNoise,length(listCycle),1)';
val = cellfun(@(x) num2str(x,'%.0f'), num2cell(decMat), 'UniformOutput', false); % convert values to strings
for i = 1:length(listNoise)
    for j = 1:length(listCycle)
        if clrMat(i,j) == 1
            clr = 'k';
        else
            clr = 'w';
        end
        text(x(i,j),y(i,j),val(i,j),'HorizontalAlignment','Center','FontSize',42,'Color',clr)
    end
end
pause(0.8);
xlabel('Duration (ms)');
ylabel('SNR (dB)');
title({'Detection Confidence'; 'Decision Matrix'});
ax = gca;
ax.XTickLabelRotation = 0;
ax.XRuler.Axle.LineStyle = 'none';
ax.YRuler.Axle.LineStyle = 'none';
ax.XAxis.MajorTickChild.LineWidth = 4;
ax.YAxis.MajorTickChild.LineWidth = 4;
set(ax,'TickDir','out','Box','off','FontSize',52,'FontWeight','bold','LineWidth',2,'Color','none','Position',ax_pos);
set(gcf,'Color','w','Position',fig_pos);
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

% Function #2: Find minimum and maximum metric values
function [val_min,val_max] = find_metric_minmax(data_metric)
    if isstruct(data_metric)
        data_metric = struct2cell(data_metric);
    end
    values = [data_metric{:}];
    val_min = min(values,[],'all');
    val_max = max(values,[],'all');
end

% Function #3: Plot Color Scale Bar Manually
function plot_manual_colorbar(xvar,yvar,dc_cutoff,custom_cmap,vmin,vmax,lbl_name)
    fnt_sz = 34;
    vmin = round(vmin,2);
    vmax = round(vmax,2);
    imagesc(xvar,yvar,ones(length(xvar),length(yvar))*dc_cutoff);
    colormap(custom_cmap); cb = colorbar('southoutside'); clim([vmin,vmax]);
    pause(0.5);
    
    ax = gca;
    ax.XRuler.Axle.LineStyle = 'none';
    ax.YRuler.Axle.LineStyle = 'none';
    set(ax,'Color','none','Position',[0.2138,0.2537,0.5455,0.6926],'Visible','off');
    set(cb,'Visible','on','XColor','none','YColor','none');
    fig_pos = [1571,30,1230,976];
    set(gcf,'Color','w','Position',fig_pos);
    
    dm_pos = cb.Position + [-0.0014 0 0 0];
    dm_ax = axes('Position',dm_pos,'Color','none','TickDir','out','Layer','bottom','YAxisLocation','right','FontSize',fnt_sz,'FontWeight','bold');
    dm_ax.YTick = [];
    dm_ax.XTick = linspace(vmin,vmax,6);
    dm_ax.XTickLabel = cellfun(@(x) num2str(x,'%.2f'),num2cell(linspace(vmin,vmax,6)),'UniformOutput',false);
    dm_ax.XLim = [vmin vmax];
    dm_ax.YColor = 'w';
    dm_ax.TickLength = [0.015,0];
    dm_ax.LineWidth = 4;
    dm_ax.Box = 'off';    
    xlabel(dm_ax,lbl_name,'FontSize',fnt_sz,'FontWeight','bold');
    uistack(dm_ax, 'bottom');
    uistack(ax, 'top');
end

% Function #4: Construct Decision Matrix Based on the Best Performing Algorithms
function [decMat,clrMat,list_uqid,list_algo] = extract_decision_matrix(matD_best,dc)
    % [1] Assign binary vectors to each set of algorithms that has a best
    %     decision confidence value
    algorithms = {'BP','ENV','STP','MTP','CWT'};
    algorithms_id = 1:length(algorithms);
    nMethod = length(dc);
    binary_decMat = cell(size(matD_best));
    for row = 1:size(matD_best,1)
        for col = 1:size(matD_best,2)
            binary_id = zeros(1,nMethod);
            for n = 1:nMethod
                dc_mat = round(dc{n},1); % consider rounded values
                if matD_best(row,col) == dc_mat(row,col)
                    binary_id(n) = 1;
                end
            end
            binary_decMat{row,col} = binary_id;
        end
    end
    % [2] Create matrix of unique number IDs and corresponding algorithm
    %     sets based on the acquired binary matrix
    uqid_decMat = cellfun(@(x) bin2dec(strjoin(string(x),"")), binary_decMat,'UniformOutput',false);
    algo_decMat = cellfun(@(x) algorithms(logical(x)),binary_decMat,'UniformOutput',false);
    % [3] Reassign reader-friendly numbers to unique number IDs
    list_uqid = unique(cell2mat(uqid_decMat));
    list_algo = cell(1,length(list_uqid));
    for i = 1:length(list_uqid)
        [map_row, map_col] = find(cell2mat(uqid_decMat) == list_uqid(i));
        list_algo{i} = algo_decMat{map_row(1),map_col(1)};
    end
    [num_algo,idx_ord] = sort(cellfun(@(x) length(x), list_algo));
    list_uqid = list_uqid(idx_ord);
    list_algo = list_algo(idx_ord);
    idx_cut = find(num_algo <= 1,1,'last');
    list_uqid = zeros(size(list_uqid));
    k = nMethod + 1;
    for i = 1:length(num_algo)
        if i <= idx_cut
            list_uqid(i) = algorithms_id(strcmpi(list_algo{i},algorithms));
        else
            list_uqid(i) = k;
            k = k + 1;
        end
    end
    if length(list_uqid) ~= length(list_algo)
        error('The number of unique number IDs and algorithm sets should be indentical.');
    end
    % [4] Construct the decision matrix
    decMat = zeros(size(matD_best));
    for row = 1:size(decMat,1)
        for col = 1:size(decMat,2)
            test = algo_decMat{row,col};
            target = cellfun(@(x) [x{:}],list_algo,'UniformOutput',false);
            identifier = strcmpi([test{:}],target);
            if sum(identifier) == 1
                decMat(row,col) = list_uqid(identifier);
            end
        end
    end
    % [5] Construct the color matrix
    clrMat = zeros(size(matD_best));
    for row = 1:size(clrMat,1)
        for col = 1:size(clrMat,2)
            num_algo = length(algo_decMat{row,col});
            if num_algo > 1
                if num_algo == 5
                    clrMat(row,col) = 7; % every algorithms
                else
                    clrMat(row,col) = 6; % multiple algorithms
                end
            else
                clrMat(row,col) = find(strcmpi(algo_decMat{row,col},algorithms)==1); % single algorithms
            end 
        end
    end
end