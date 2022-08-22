function makeLegendToImage(figHandle, legHandle, figType)
    %% Function 'makeLegendToImage'
    % DESCRIPTION
    % Creates a manual legend for specific figures
    
    % INPUT
    %    Variable       Data Type          Description
    % 1. figHandle      [object]         : figure handle
    % 2. legHandle      [object]         : legend handle
    % 3. figType        [string]         : plot types for the legend
    
    % Source: https://stackoverflow.com/questions/18117664/how-can-i-show-only-the-legend-in-matlab
    % Last modified by SungJun Cho, August 09, 2022
    %   1. Added `figType` argument
    %   2. Deleted saving functionality
    %% Generate Figure Legend
    % [1] Make All Contents in the Figure Invisible
    allHandles = findall(figHandle, 'type', figType);
    for i = 1:length(allHandles)
        allHandles(i).XData = NaN; %ignore warnings
        allHandles(i).YData = NaN;
    end
    % [2] Make Axis Invisible
    axis off
    % [3] Relocate Legend to Lower Left Corner of the Figure Window
    legHandle.Units = 'pixels';
    boxLineWidth = legHandle.LineWidth;
    % [4] Ensure Legend Box Position is Adequate for Saving
    legHandle.Position = [6 * boxLineWidth, 6 * boxLineWidth, ...
        legHandle.Position(3), legHandle.Position(4)];
    legLocPixels = legHandle.Position;
    % [5] Adjut Figure Window to Fit Legend
    figHandle.Units = 'pixels';
    figHandle.InnerPosition = [1, 1, legLocPixels(3) + 12 * boxLineWidth, ...
        legLocPixels(4) + 12 * boxLineWidth];
end