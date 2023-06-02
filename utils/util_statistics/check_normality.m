function [pvals] = check_normality(varargin)
    %% Function: 'check_normality'
    % DESCRIPTION
    % This function checks whether the input vector has a normal
    % distribution using the Kolmogorov-Smirnov (KS) test.
    
    % INPUT
    %    Variable    Data Type            Description
    % 1. varargin    [1 x N doubles]    : data vector array
    
    % OUTPUT
    %    Variable    Data Type            Description
    % 1. pvals       [1 x N double]     : p-values of the KS test for each input data
    
    % Written by SungJun Cho, May 19, 2023
    %% Test Normality Assumption
    pvals = zeros(1,nargin);
    for i = 1:nargin
        % [1] Normalize the Data
        data_nm = normalize(varargin{i});
        % [2] Compre eCDF to the standard normal CDF
        [h,p] = kstest(data_nm);
        if h == 0
            message = ['Null hypothesis accepted (p=' num2str(p) '). ' ...
                'Distribution can be considered normal. \n'];
        elseif h == 1
            message = ['Null hypothesis rejected (p=' num2str(p) '). ' ...
                'Distributions cannot be considered normal. \n'];
        end
        fprintf(['Group #' num2str(i) ': ' message]);
        pvals(i) = p;
    end
end