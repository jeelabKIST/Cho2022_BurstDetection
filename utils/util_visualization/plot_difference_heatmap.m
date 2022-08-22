function plot_difference_heatmap(xvar,yvar,matV_best,matV_poor,custom_cmap,vmin,vmax,axis_params)
    %% Function 'plot_difference_heatmap'
    % DESCRIPTION
    % This function plots difference between heatmaps of best and worst
    % performing algorithms.
    
    % INPUT
    %    Variable      Data Type            Description
    % 1. xvar          [1 x N vector]     : values of x-axis
    % 2. yvar          [1 x N vector]     : values of y-axis
    % 3. matV_best     [N x N matrix]     : input hetamap of best performing algorithms
    % 4. matV_poor     [N x N matrix]     : input heatmap of worst performing algorithms
    % 5. custom_cmap   [N x N matrix]     : custom color map
    % 6. vmin          [integer Z]        : minimum value for colorbar
    % 7. vmax          [integer Z]        : maximum value for colorbar
    % 8. axis_params   [struct]           : axis parameters for visualization
    
    % Written by SungJun Cho, October 11, 2021
    % Last modified on October 29, 2021
    %% Set Axis Parameters
    fnt_sz = axis_params.fnt_sz;
    txt_sz = axis_params.txt_sz;
    fig_pos = axis_params.fig_pos;
    ax_pos = axis_params.ax_pos;
    annot_fmt = axis_params.annot_fmt;
    annot_rng = axis_params.annot_rng;
    cbar_opt = axis_params.cbar_opt;
    cbar_loc = axis_params.cbar_loc;
    cbar_lbl = axis_params.cbar_lbl;
    vint = axis_params.vint;
    xlbl_opt = axis_params.xlbl_opt;
    ylbl_opt = axis_params.ylbl_opt;
    %% Visualize Heatmap
    % [1] Plot Heatmap
    matDiff = matV_best - matV_poor;
    imagesc(xvar,yvar,matDiff);
    xticks(xvar);
    yticks(yvar);
    if strcmp(xlbl_opt, 'on')
        xvar_lbl = categorical(cellfun(@(x) mod(x,xvar(1)) == 0, num2cell(xvar)) .* xvar);
        xvar_lbl(xvar_lbl == categorical(0)) = ' ';
        xticklabels(xvar_lbl);
    end
    if strcmp(ylbl_opt,'on')
        ylabel('SNR (dB)');
    else
        yticklabels(categorical(NaN(1,length(yvar))));
    end
    colormap(custom_cmap); cb = colorbar(cbar_loc); caxis([vmin,vmax]);
    axis xy;
    x = repmat(xvar,length(yvar),1);
    y = repmat(yvar,length(xvar),1)';
    val = num2cell(matDiff);
    val = cellfun(@(k) num2str(k,annot_fmt), val, 'UniformOutput', false); % convert values to strings
    for i = 1:length(yvar)
        for j = 1:length(xvar)
            if str2double(val{i,j}) < annot_rng(1) && str2double(val{i,j}) > annot_rng(2)
                clr = 'k';
            else
                clr = 'w';
            end
            text(x(i,j),y(i,j),val(i,j),'HorizontalAlignment','Center','FontSize',txt_sz,'Color',clr);
        end
    end
    pause(0.5);
    xlabel('Duration (ms)');
    ax = gca;
    ax.XRuler.Axle.LineStyle = 'none';
    ax.YRuler.Axle.LineStyle = 'none';
    set(ax,'TickDir','out','Box','off','FontSize',fnt_sz,'FontWeight','bold','LineWidth',2,'Color','none',ax_pos.pos_type,ax_pos.pos_coord);
    ax.XAxis.MajorTickChild.LineWidth = 4;
    ax.YAxis.MajorTickChild.LineWidth = 4;
    set(cb,'Visible',cbar_opt,'XColor','none','YColor','none');
    set(gcf,'Color','w','Position',fig_pos);
    % [2] Customize Colorbar
    if strcmp(cbar_opt,'on')
        dm_pos = cb.Position + [-0.0014 0 0 0];
        dm_ax = axes('Position',dm_pos,'Color','none','TickDir','out','Layer','bottom','YAxisLocation','right','FontSize',fnt_sz,'FontWeight','bold');
        if strcmp(cbar_loc, 'southoutside')
            dm_ax.YTick = [];
            dm_ax.XTick = vmin:vint:vmax;
            dm_ax.XLim = [vmin vmax];
            dm_ax.YColor = 'w';
            xlabel(dm_ax,cbar_lbl,'FontSize',fnt_sz,'FontWeight','bold');
        else
            dm_ax.XTick = [];
            dm_ax.YTick = vmin:vint:vmax;
            dm_ax.YLim = [vmin vmax];
            dm_ax.XColor = 'w';
            ylabel(dm_ax,cbar_lbl,'FontSize',fnt_sz,'FontWeight','bold');
        end
        dm_ax.TickLength = [0.015,0];
        dm_ax.LineWidth = 4;
        dm_ax.Box = 'off';    
        uistack(dm_ax, 'bottom');
        uistack(ax, 'top');
    end
end