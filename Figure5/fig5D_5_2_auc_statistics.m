%% Set Parameters
% The metrics available when reading 'AUC.mat' are: (1) precision, (2)
% recall, (3) f1_score, (4) concurrence, and (5) detection_confidence.
% To replicate the figures, use:
%   * (5)     for Figure 5D
%   * (1)-(4) for Figure 5-2
metric = 'detection_confidence';
metric_name = 'Detection Confidence';
algorithms = {'BP','ENV','S-STFT','MTP','CWT'};
nMethod = length(algorithms);
fprintf(['\n* Selected Metric: ' upper(metric) '\n']);
%% Load Data
% [1] Get Area Under Curve (AUC)
auc = load('AUC_stage2.mat').AUC.(metric);
% [2] Select Trials with Human Annotations
nTotTrial = size(auc,1);
sel_trials = sort(cat(2,1:8:nTotTrial,5:8:nTotTrial)); % 1st trials of Day 1 & Day 2 for 8 mice
auc = auc(sel_trials,:);
% [3] Compute Inverse AUC and Its Statistics
inverse_auc = 1./auc;
auc_mean = mean(inverse_auc,1);
auc_std = std(inverse_auc,1);
nTrial = size(inverse_auc,1);
%% Visualize Burst Detection Confidence AUC
% [1] Set Visualization Parameters
okabe_ito = [[0,114,178];
             [213,94,0];
             [0,158,115];
             [204,121,167];
             [240,228,66]]./255;
grey = [97,97,97]./255;
vmin = round(min(inverse_auc, [], 'all'),1)-0.1;
vmax = round(max(inverse_auc, [], 'all'),1)+0.1;
if mod(vmin,0.2) ~= 0
    vmin = vmin - 0.1;
end
if mod(vmax,0.2) ~= 0
    vmax = vmax + 0.1;
end
alpha = 0.4;
% [2] Plot Bar Graphs: AUC across Trials
figure(); hold on;
x_axis = 1:nMethod;
b = bar(x_axis, auc_mean);
b.FaceColor = 'flat';
for n = 1:nMethod
    b.CData(n,:) = okabe_ito(n,:);
end
s = cell(size(x_axis));
for n = 1:nMethod
    rnd = unifrnd(-0.2,0.2,1,nTrial);
    s{n} = scatter(ones(1,nTrial)*n+rnd,inverse_auc(:,n),80,grey,'filled');
    s{n}.MarkerFaceAlpha = 0.65;
    s{n}.LineWidth = 1.5;
end
er = cell(size(x_axis));
for n = 1:nMethod
    er{n} = errorbar(x_axis(n),auc_mean(n),auc_std(n),auc_std(n));
    er{n}.Color = okabe_ito(n,:);
    if n == nMethod
        er{n}.Color = [209,206,23]./255; % darker yellow for visibility
    end
end
xticks(x_axis);
xticklabels(algorithms);
yticks(linspace(vmin,vmax,5));
xlim([0.2,nMethod+0.8]);
ylim([vmin,vmax]);
xlabel('Algorithms');
ylabel({'Inverse Area Under Curve'; '(1/AUC)'});
% (A) Figure Settings
set(b,'EdgeColor','none','BaseValue',vmin+0.004,'ShowBaseLine','off','FaceAlpha',alpha);
set([er{:}],'CapSize',25,'LineStyle','none','LineWidth',10);
% set(er.Bar, 'ColorType', 'truecoloralpha', 'ColorData', [er.Line.ColorData(1:3); 255*alpha]);
% set(er.Cap, 'EdgeColorType', 'truecoloralpha', 'EdgeColorData', [er.Cap.EdgeColorData(1:3); 255*alpha]);
set(gca,'FontSize',35,'LineWidth',5,'TickDir','out','TickLength',[0.015, 0.025],'Box','off','FontWeight','bold',...
    'Position',[0.221648616102467,0.146554506483428,0.73623345369296,0.778445493516572]);
set(gcf,'Color','w','Position',[703,199,831,743]);
%% Statistical Test
% [1] One-Way ANOVA
group = algorithms;
[p,tbl,stats] = anova1(inverse_auc,group,'off');
% [2] Multiple Comparisons Test: Tukey's HSD
[c,m_est,~,gnames] = multcompare(stats,'CType','tukey-kramer','Display','off');
mc_tbl = array2table(c,'VariableNames',{'Group','Control Group','Lower Limit','Difference','Upper Limit','p-value'});
mc_tbl.('Group') = gnames(mc_tbl.('Group'));
mc_tbl.('Control Group') = gnames(mc_tbl.('Control Group'));
% [3] Print Statistical Test Summary
fprintf('* One-Way ANOVA with HSD (multiple comparisons test): \n');
disp(mc_tbl);