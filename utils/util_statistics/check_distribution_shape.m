function check_distribution_shape(dataset, standardize, verbose)
    %% Function: 'check_distribution_shape'
    % DESCRIPTION
    % This function checks pair-wise if distributions come from the same continuous
    % distribution using the Kolmogorov-Smirnov test. Distributions can also be 
    % plotted for a qualitative/visual comparison.
    
    % INPUT
    %    Variable    Data Type            Description
    % 1. dataset     [1 x N cells]      : input data with each cell containing 1D vector array
    % 2. standardize [boolean]          : whether to standardize the data
    %                                     default) true
    % 3. verbose     [boolean]          : whether to plot the data distributions
    %                                     default) true
    
    % Written by SungJun Cho, May 19, 2023
    %% Inspect Distribution Shapes
    if nargin < 3
        verbose = true;
    end
    if nargin < 2
        standardize = true;
    end
    % [1] Preprocess the Data
    nData = length(dataset);
    if standardize
        dataset = cellfun(@(x) normalize(x),dataset,'UniformOutput',false);
    end
    % [2] Two-sample Kolmogorov-Smirnov Test
    permutations = nchoosek(1:nData, 2);
    for i = 1:size(permutations,1)
        idx1 = permutations(i,1);
        idx2 = permutations(i,2);
        [h,p] = kstest2(dataset{idx1},dataset{idx2});
        if h == 0
            message = ['Null hypothesis accepted (p=' num2str(p) '). ' ...
                       'Distributions are from the same continuous distribution. \n'];
        else
            message = ['Null hypothesis accepted (p=' num2str(p) '). ' ...
                       'Distributions are from different continuous distributions. \n'];
        end
        fprintf(['Group ' num2str(idx1) ' vs ' num2str(idx2) ': ' message]);
    end
    % [3] Plot Distributions
    if verbose
        figure(); hold on;
        lg_lbl = cell(1,nData);
        cmap = colormap(lines(nData));
        for n = 1:nData
            [f,xi] = ksdensity(dataset{n});
            plot(xi,f,'Color',cmap(n,:),'LineWidth',2);
            lg_lbl{n} = ['Group ' num2str(n)];
        end
        xlabel('Data');
        ylabel('Probability Density');
        legend(lg_lbl);
        set(gca,'Fontsize',14,'LineWidth',2,'TickDir','out');
        set(gcf,'Color','white')
    end
end
