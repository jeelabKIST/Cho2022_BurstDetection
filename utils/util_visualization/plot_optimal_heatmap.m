function plot_optimal_heatmap(xvar, yvar, matC, matV, custom_cmap, vmin, vmax, axis_params)
    %% Function 'plot_optimal_heatmap'
    % DESCRIPTION
    % This function plots heatmaps of the best and worst performing algorithms.
    
    % INPUT
    %    Variable      Data Type            Description
    % 1. xvar          [1 x N vector]     : values of x-axis
    % 2. yvar          [1 x N vector]     : values of y-axis
    % 3. matC          [N x N matrix]     : matrix consisting of color labels
    % 4. matV          [N x N matrix]     : input heatmap
    % 5. custom_cmap   [N x N matrix]     : custom color map
    % 6. vmin          [integer Z]        : minimum value for colorbar
    % 7. vmax          [integer Z]        : maximum value for colorbar
    % 8. axis_params   [struct]           : axis parameters for visualization
    
    % Written by SungJun Cho, October 11, 2021
    % Last modified on January 29, 2023
    %% Set Axis Parameters
    fnt_sz = axis_params.fnt_sz;
    txt_sz = axis_params.txt_sz;
    fig_pos = axis_params.fig_pos;
    ax_pos = axis_params.ax_pos;
    annot_fmt = axis_params.annot_fmt;
    cbar_opt = axis_params.cbar_opt;
    xlbl_opt = axis_params.xlbl_opt;
    ylbl_opt = axis_params.ylbl_opt;
    %% Visualize Heatmap
    imagesc(xvar,yvar,matC);
    xticks(xvar);
    yticks(yvar);
    if strcmp(xlbl_opt,'on')
        xvar_lbl = categorical(cellfun(@(x) mod(x,xvar(1)) == 0, num2cell(xvar)) .* xvar);
        xvar_lbl(xvar_lbl == categorical(0)) = ' ';
        xticklabels(xvar_lbl);
    end
    if strcmp(ylbl_opt,'on')
        ylabel('SNR (dB)');
    else
        yticklabels(categorical(NaN(1,length(yvar))));
    end
    colormap(custom_cmap); cb = colorbar; caxis([vmin,vmax]);
    axis xy;
    x = repmat(xvar,length(yvar),1);
    y = repmat(yvar,length(xvar),1)';
    val = num2cell(matV); % extract values into cells
    val = cellfun(@(x) num2str(x,annot_fmt), val, 'UniformOutput', false); % convert values to strings
    for i = 1:length(yvar)
        for j = 1:length(xvar)
            if matC(i,j) == 1
                clr = 'k';
            else
                clr = 'w';
            end
            text(x(i,j),y(i,j),val(i,j),'HorizontalAlignment','Center','FontSize',txt_sz,'Color',clr);
        end
    end
    pause(0.8);
    xlabel('Duration (ms)');
    ylabel(cb,'Algorithms');
    ax = gca;
    ax.XTickLabelRotation = 0;
    ax.XRuler.Axle.LineStyle = 'none';
    ax.YRuler.Axle.LineStyle = 'none';
    set(ax,'TickDir','out','Box','off','FontSize',fnt_sz,'FontWeight','bold','LineWidth',2,'Color','none',ax_pos.pos_type,ax_pos.pos_coord);
    ax.XAxis.MajorTickChild.LineWidth = 4;
    ax.YAxis.MajorTickChild.LineWidth = 4;
    set(cb,'Visible',cbar_opt,'XColor','none','YColor','none');
    set(gcf,'Color','w','Position',fig_pos);
end