function [pval] = check_equal_variance(dataset, test_type)
    %% Function: 'check_equal_variance'
    % DESCRIPTION
    % This function checks whether the distributions of the input data have
    % equal variances using multiple-samples tests for equal variances.
    
    % INPUT
    %    Variable    Data Type            Description
    % 1. dataset     [1 x N cells]      : input data with each cell containing 1D vector array
    % 2. test_type   [string]           : type of hypothesis test
    %                                     For the available options, refer to the function `vartestn`.
    
    % OUTPUT
    %    Variable    Data Type            Description
    % 1. pval        [float]            : p-values of the test
    
    % Written by SungJun Cho, May 19, 2023
    %% Test Equal Variance Assumption
    % [1] Preprocess the Data
    nData = length(dataset);
    if ~iscolumn(dataset{1})
        dataset = cellfun(@(x) x',dataset,'UniformOutput',false);
    end
    data = vertcat(dataset{:});
    % [2] Set Grouping Variable
    group = zeros(size(data));
    i = 1;
    j = 0;
    for n = 1:nData
        j = j + length(dataset{n});
        group(i:j) = n;
        i = i + length(dataset{n});
    end
    % [3] Conduct Multiple-Sample Test for Equal Variances
    pval = vartestn(data,group,'TestType',test_type,'Display','off');
    if pval >= 0.05
        fprintf(['Null hypothesis accepted (p=' num2str(pval) '). ' ...
                 'Distributions have equal variances. \n']);
    else
        fprintf(['Null hypothesis rejected (p=' num2str(pval) '). ' ...
                 'Distributions have unequal variances. \n']);
    end
end