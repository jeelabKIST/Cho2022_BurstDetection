%% Configure Library Path
util_path = genpath('/Users/scho/Neuroscience_KIST/Cho2022_BurstDetection/utils');
data_path = genpath('/Users/scho/Neuroscience_KIST/Cho2022_BurstDetection/data');
addpath(util_path);
addpath(data_path);
%% Load Data
DC_ALPHA = load('DC_alpha.mat').detection_confidence;
DC_BETA = load('DC_beta.mat').detection_confidence;
DC_GAMMA = load('DC_gamma.mat').detection_confidence;
[DC_ALPHA,DC_BETA,DC_GAMMA] = convert2percentage(DC_ALPHA,DC_BETA,DC_GAMMA); % convert to percentage scale
Fs = 512;
listNoise = -10:2:10;
nMethod = 5;
%% Plot Decision Matrix
% [1] Set Visualization Parameters
convert_opt = false;
round_opt = 'decimal';
okabe_ito_expand = [[0,114,178];
             [213,94,0];
             [0,158,115];
             [204,121,167];
             [240,228,66];
             [182,188,191];
             [0,0,0]]./255;
fig_pos = [1571,30,1230,976];
ax_pos = [0.2538,0.1837,0.5455,0.6926];
% [2] Visualize Decision Matrix
DC = {DC_ALPHA, DC_BETA, DC_GAMMA};
cnt_f = [9,25,41];
band_names = {'Alpha','Beta','Low Gamma'};
for n = 1:length(DC)
    % (A) Compute Decision Matrix
    listCycle = round(((3:12).*Fs/cnt_f(n)./Fs)*1000); % convert length of cycles to milliseconds
    [matD_best,~,~,~] = get_best_and_poor_heatmaps(listNoise,listCycle,DC{n},convert_opt,round_opt);
    [decMat,clrMat,list_uqid,list_algo] = extract_decision_matrix(matD_best,DC{n});
    if sum(list_uqid > nMethod) ~= 0
        for id_idx = find((cellfun(@(x) length(x), list_algo) > 1) == 1)
            fprintf(['Multi-algorithm ID - ' num2str(list_uqid(id_idx)) ': ' sprintf('%s ', string(list_algo{id_idx})) '\n']);
        end
    end
    % (B) Construct Figure
    figure();
    imagesc(listCycle,listNoise,clrMat);
    xticks(listCycle);
    yticks(listNoise);
    xvar_lbl = categorical(listCycle);
    xvar_lbl(setdiff(1:length(listCycle),1:3:length(listCycle))) = ' ';
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
    pause(1);
    xlabel('Duration (ms)');
    ylabel('SNR (dB)');
    title(band_names{n});
    ax = gca;
    ax.XTickLabelRotation = 0;
    ax.XRuler.Axle.LineStyle = 'none';
    ax.YRuler.Axle.LineStyle = 'none';
    ax.XAxis.MajorTickChild.LineWidth = 4;
    ax.YAxis.MajorTickChild.LineWidth = 4;
    set(ax,'TickDir','out','Box','off','FontSize',52,'FontWeight','bold','LineWidth',2,'Color','none','Position',ax_pos);
    set(gcf,'Color','w','Position',fig_pos);
end
%% Appendix: In-Script Functions
% Function #1: Convert Matrices in Cell Array to Percentage Scale
function [varargout] = convert2percentage(varargin)
    varargout = cell(1,nargin);
    for n = 1:nargin
        varargout{n} = cellfun(@(mat) mat.*100,varargin{n},'UniformOutput',false);
    end
end

% Function #2: Construct Decision Matrix Based on the Best Performing Algorithms
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
