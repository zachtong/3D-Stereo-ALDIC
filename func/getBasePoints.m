function Base_Points_2D = getBasePoints(imageLeft,maskLeft)
%GETBASEPOINTS Interactive selection of 3 base points to define coordinate system
%
% INPUTS:
%   imageLeft - Input image for point selection
%
% OUTPUTS:
%   Base_Points_2D - 3x2 matrix containing [O; X; Y] point coordinates
%
% DESCRIPTION:
%   This function allows user to interactively select 3 points on an image
%   to define a new coordinate system: Origin (O), X-axis direction, Y-axis direction

%% Input validation
if nargin < 1
    error('getBasePoints: Input image is required');
end

if isempty(imageLeft)
    error('getBasePoints: Input image is empty');
end

%% Prepare image for display
% Handle different image formats
if size(imageLeft, 3) > 1
    imageLeft = rgb2gray(imageLeft);  % Convert to grayscale if RGB
end

% Handle different data types
if any(imageLeft(:) > 255)
    % 16-bit or higher
    displayImg = uint16(imageLeft);
else
    % 8-bit
    displayImg = uint8(imageLeft);
end

%% Interactive point selection
figure('Name', 'Base Points Selection', 'NumberTitle', 'off');
imshow(displayImg);  % Remove transpose - display image correctly
axis on; grid on;

% --- mask overlay (new code) ---
if nargin > 1 && ~isempty(maskLeft)
    hold on;
    
    % find mask boundary
    B = bwboundaries(maskLeft, 'noholes');
    
    % create a semi-transparent overlay for the mask interior
    h = imshow(cat(3, zeros(size(maskLeft)), ones(size(maskLeft)), zeros(size(maskLeft)))); % green color
    set(h, 'AlphaData', 0.3 * maskLeft); % set transparency for the mask region
    
    % draw a thick, opaque boundary
    for k = 1:length(B)
        boundary = B{k};
        plot(boundary(:,2), boundary(:,1), 'r-', 'linewidth', 3); % red boundary
    end
end
% --- end of new code ---












title({'Click 3 points to define coordinate system:', ...
       '1st: Origin (O)', '2nd: X-axis direction', '3rd: Y-axis direction'}, ...
       'FontWeight', 'normal', 'FontSize', 14);

Base_Points_2D = zeros(3, 2);
pointLabels = {'Origin (O)', 'X-axis direction', 'Y-axis direction'};
colors = {'ro', 'bo', 'go'};

for i = 1:3
    % Update title for current point
    title({sprintf('Click point %d: %s', i, pointLabels{i}), ...
           'Right-click to redo previous point'}, ...
           'FontWeight', 'normal', 'FontSize', 14);
    
    % Get point coordinates
    [x, y, button] = ginput(1);
    
    % Handle right-click (redo previous point)
    if button == 3 && i > 1
        % Remove previous point marker
        children = get(gca, 'Children');
        if length(children) > 1  % More than just the image
            delete(children(1));  % Remove last added marker
        end
        i = i - 1;  % Go back one step
        continue;
    end
    
    % Validate point is within image bounds
    if x < 1 || x > size(imageLeft, 2) || y < 1 || y > size(imageLeft, 1)
        warning('Point is outside image bounds. Please click within the image.');
        i = i - 1;  % Retry current point
        continue;
    end
    
    % Store coordinates
    Base_Points_2D(i, :) = [x, y];
    
    % Add visual marker
    hold on;
    plot(x, y, colors{i}, 'MarkerSize', 8, 'LineWidth', 2);
    text(x+10, y-10, pointLabels{i}, 'Color', colors{i}(1), ...
         'FontSize', 10, 'FontWeight', 'bold', 'BackgroundColor', 'white');
    
    % Display coordinates
    fprintf('Point %d (%s): (%.3f, %.3f)\n', i, pointLabels{i}, x, y);
end

%% Validation and final display
% Check if points are collinear
vec1 = Base_Points_2D(2,:) - Base_Points_2D(1,:);
vec2 = Base_Points_2D(3,:) - Base_Points_2D(1,:);
crossProduct = vec1(1)*vec2(2) - vec1(2)*vec2(1);

if abs(crossProduct) < 1e-6
    warning('getBasePoints: Selected points are nearly collinear. This may cause numerical issues.');
end

% Draw coordinate system visualization
hold on;
% X-axis (red arrow)
quiver(Base_Points_2D(1,1), Base_Points_2D(1,2), ...
       vec1(1), vec1(2), 0, 'r', 'LineWidth', 2, 'MaxHeadSize', 0.3);
% Y-axis (green arrow)  
quiver(Base_Points_2D(1,1), Base_Points_2D(1,2), ...
       vec2(1), vec2(2), 0, 'g', 'LineWidth', 2, 'MaxHeadSize', 0.3);

title('Selected coordinate system (Red: X-axis, Green: Y-axis)', ...
      'FontWeight', 'normal', 'FontSize', 12);

% Final confirmation
fprintf('\n=== Selected Base Points ===\n');
fprintf('Origin (O): (%.3f, %.3f)\n', Base_Points_2D(1,1), Base_Points_2D(1,2));
fprintf('X-axis:     (%.3f, %.3f)\n', Base_Points_2D(2,1), Base_Points_2D(2,2));
fprintf('Y-axis:     (%.3f, %.3f)\n', Base_Points_2D(3,1), Base_Points_2D(3,2));
fprintf('============================\n\n');

% Optional: Ask for confirmation
response = questdlg('Are you satisfied with the selected points?', ...
                   'Confirm Selection', 'Yes', 'No', 'Yes');
if strcmp(response, 'No')
    close;
    Base_Points_2D = getBasePoints(imageLeft,maskLeft);  % Recursive call to restart
    return;
end

end