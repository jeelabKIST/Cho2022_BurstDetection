% Written by SungJun Cho, October 25, 2021
% Last modified on August 8, 2022

function h = helper
    %% Function 'helper'
    % DESCRIPTION
    % Returns handles to local functions and initiates all the functions
    % below
    h.unpack_data = @unpack_data;
    h.pack_data = @pack_data;
    h.convert_to_percent = @convert_to_percent;
    h.dec2int = @dec2int;
    h.round_dec = @round_dec;
    h.find_binary_idx = @find_binary_idx;
end

function [mat_bp,mat_ev,mat_stp,mat_mtp,mat_cwt] = unpack_data(dataset)
    %% Function 'unpack_data'
    % DESCRIPTION
    % Extracts a dataset of struct or cell type to double
    
    % INPUT
    %    Variable    Data Type           Description
    % 1. dataset     [struct / cell]   : dataset consisting of heatmaps of each algorithm
    
    % OUTPUT
    %    Variable    Data Type           Description
    % 1. mat_bp      [N x N matrix]    : heatmap of BP algorithm
    % 2. mat_ev      [N x N matrix]    : heatmap of ENV algorithm
    % 3. mat_stp     [N x N matrix]    : heatmap of S-STFT algorithm
    % 4. mat_mtp     [N x N matrix]    : heatmap of MTP algorithm
    % 5. mat_cwt     [N x N matrix]    : heatmap of CWT algorithm
    %% Dataset Type: Structure
    if isa(dataset,'struct')
        mat_bp = dataset.bp;
        mat_ev = dataset.ev;
        mat_stp = dataset.stp;
        mat_mtp = dataset.mtp;
        mat_cwt = dataset.cwt;
    end
    %% Dataset Type: Cell
    if isa(dataset,'cell')
        mat_bp = dataset{1};
        mat_ev = dataset{2};
        mat_stp = dataset{3};
        mat_mtp = dataset{4};
        mat_cwt = dataset{5};
    end
end

function [package] = pack_data(varargin)
    %% Function 'pack_data'
    % DESCRIPTION
    % Packs the input data into a single cell array
    
    % OUTPUT
    %    Variable    Data Type           Description
    % 1. package     [1 x N cell]      : cell array containing input data
    package = cell(1,nargin);
    for i = 1:nargin
        package{i} = varargin{i};
    end
end

function [varargout] = convert_to_percent(varargin)
    %% Function 'convert_to_percent'
    % DESCRIPTION
    % Converts each element of input vectors and matrices into percentage scale
    varargout = cell(1,nargin);
    for i = 1:nargin
        varargout{i} = varargin{i}.*100;
    end
end

function [varargout] = dec2int(varargin)
    %% Function 'dec2int'
    % DESCRIPTION
    % Converts each element of input vectors and matrices into nearest integers
    varargout = cell(1,nargin);
    for i = 1:nargin
        varargout{i} = round(varargin{i});
    end
end

function [varargout] = round_dec(varargin)
    %% Function 'round_dec'
    % DESCRIPTION
    % Converts each element of input vectors and matrices into first digit decimals    
    varargout = cell(1,nargin);
    for i = 1:nargin
        varargout{i} = round(varargin{i},1);
    end
end

function [idx_start,idx_end] = find_binary_idx(test_vector)
    %% Function 'find_binary_idx'
    % DESCRIPTION
    % Finds indices where 1 (True) starts and ends from the given array
    
    % INPUT
    %    Variable      Data Type              Description
    % 1. test_vector   [double / logical]   : vector array of binary values (0 or 1)
    
    % OUTPUT
    %    Variable      Data Type              Description
    % 1. idx_start     [1 x N matrix]       : vector array of indices where 1 starts
    % 2. idx_end       [1 x N matrix]       : vector array of indices where 1 ends
    idx_start = strfind([0 test_vector],[0 1]);
    idx_end = strfind([test_vector 0], [1 0]);
end