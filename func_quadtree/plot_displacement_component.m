function plot_displacement_component(Img, elementsFEM, coordinatesFEMWorldDef, disp_data, um2px, transparency, plotParams)
% PLOT_DISPLACEMENT_COMPONENT: A helper function to plot a single displacement field.
% This function encapsulates the plotting logic to avoid code duplication.
%
% INPUTS:
%   Img                   - Path to the background image.
%   elementsFEM           - FE mesh elements.
%   coordinatesFEMWorldDef- FE mesh coordinates (in pixels).
%   disp_data             - Displacement data to be plotted.
%   um2px                 - Conversion factor from micrometers to pixels.
%   transparency          - Transparency of the displacement overlay.
%   plotParams            - A struct with plotting parameters, containing:
%                             .caxis: Color axis limits (e.g., [-10, 10] or 'auto').
%                             .colormap: Colormap to use (e.g., 'turbo' or 'jet').

% Create a new figure for the plot
figure; 
ax1 = axes;

% Display the background image
try
    % Use imshow for standard image formats
    imshow(flipud(imread(Img)), 'InitialMagnification', 'fit');
catch
    % Fallback to surf for other data types, which is a pragmatic way
    % to handle different kinds of image data.
    surf(flipud(imread(Img)), 'EdgeColor', 'none', 'LineStyle', 'none');
end
axis on; axis equal; axis tight; box on; 
set(gca, 'fontSize', 18);
view(2);
set(gca, 'ydir', 'normal');
hold on;

% Create a second axes for the displacement overlay
ax2 = axes;
h2 = show([], elementsFEM(:,1:4), coordinatesFEMWorldDef / um2px, disp_data, 'NoEdgeColor');

set(gca, 'fontSize', 18);
set(gca, 'ydir', 'reverse');
view(2);
box on;
axis equal;
axis tight;

% Apply plotting parameters
alpha(h2, transparency);
colormap(ax2, plotParams.colormap);
caxis(ax2, plotParams.caxis);

% Link the two axes so they zoom and pan together
linkaxes([ax1, ax2]);

% Configure axis visibility and colormaps
ax2.Visible = 'off'; 
ax2.XTick = []; 
ax2.YTick = [];
colormap(ax1, 'gray'); % Background is always gray

% Set final layout and position
set([ax1, ax2], 'Position', [.17 .11 .685 .815]);
ax1.Visible = 'on';

% Convert axis labels from pixels to physical units
xticklabels(ax1, num2cell(round(um2px * ax1.XTick * 10) / 10, length(ax1.XTick))');
yticklabels(ax1, num2cell(round(um2px * ax1.YTick * 10) / 10, length(ax1.YTick))');

% Add a colorbar for the displacement data
colorbar('Position', [.17 + 0.685 + 0.012 .11 + .128 .03 .557]);

hold off;

end