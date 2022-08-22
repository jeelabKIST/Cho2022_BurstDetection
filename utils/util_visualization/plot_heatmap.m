function plot_heatmap(xvar,yvar,data,custom_cmap,vmin,vmax,title_name,axis_params)
    %% Function 'plot_heatmap'
    % DESCRIPTION
    % This function plots a heatmap of a specific metric type for each algorithm.
    
    % INPUT
    %    Variable      Data Type            Description
    % 1. xvar          [1 x N vector]     : values of x-axis
    % 2. yvar          [1 x N vector]     : values of y-axis
    % 3. data          [N x N matrix]     : input hetamap
    % 4. custom_cmap   [N x N matrix]     : custom color map
    % 5. vmin          [integer Z]        : minimum value for colorbar
    % 6. vmax          [integer Z]        : maximum value for colorbar
    % 7. title_name    [string]           : title name of the plot
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
    vint = axis_params.vint;
    cbar_opt = axis_params.cbar_opt;
    cbar_lbl = axis_params.cbar_lbl;
    ylbl_opt = axis_params.ylbl_opt;
    %% Visualize Heatmap
    % [1] Plot Heatmap
    imagesc(xvar,yvar,data);
    xvar_lbl = categorical(cellfun(@(x) mod(x,xvar(1)) == 0, num2cell(xvar)) .* xvar);
    xvar_lbl(xvar_lbl == categorical(0)) = ' ';
    xticks(xvar); xticklabels(xvar_lbl);
    yticks(yvar);
    if strcmp(ylbl_opt,'on')
        ylabel('SNR (dB)');
    else
        yticklabels(categorical(NaN(1,length(yvar))));
    end
    colormap(custom_cmap); cb = colorbar; caxis([vmin,vmax]);
    axis xy;
    x = repmat(xvar,length(yvar),1);
    y = repmat(yvar,length(xvar),1)';
    val = num2cell(data); % extract values into cells
    val = cellfun(@(x) num2str(x,annot_fmt), val, 'UniformOutput', false); % convert values to strings
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
    pause(0.8);
    xlabel('Duration (ms)'); title(title_name);
    ax = gca;
    ax.XRuler.Axle.LineStyle = 'none';
    ax.YRuler.Axle.LineStyle = 'none';
    set(ax,'TickDir','out','Box','off','FontSize',fnt_sz,'FontWeight','bold','LineWidth',4,'Color','none',ax_pos.pos_type,ax_pos.pos_coord);
    set(cb,'Visible',cbar_opt,'XColor','none','YColor','none');
    set(gcf,'Color','w','Position',fig_pos);
    % [2] Customize Colorbar
    if strcmp(cbar_opt,'on')
        dm_pos = cb.Position + [-0.0014 0 0 0];
        dm_ax = axes('Position',dm_pos,'Color','none','TickDir','out','Layer','bottom','YAxisLocation','right','FontSize',fnt_sz,'FontWeight','bold');
        dm_ax.XTick = [];
        dm_ax.YTick = vmin:vint:vmax;
        dm_ax.YLim = [vmin vmax];
        dm_ax.XColor = 'w';
        dm_ax.TickLength = [0.015,0];
        dm_ax.LineWidth = 4;
        dm_ax.Box = 'off';    
        ylabel(dm_ax,cbar_lbl,'FontSize',fnt_sz,'FontWeight','bold');
        uistack(dm_ax, 'bottom');
        uistack(ax, 'top');
    end
end

% Helpful References:
% https://stackoverflow.com/questions/25838789/remove-only-axis-lines-without-affecting-ticks-and-tick-labels
% https://www.mathworks.com/matlabcentral/answers/352117-how-can-i-access-the-axle-property-in-an-axis-ruler