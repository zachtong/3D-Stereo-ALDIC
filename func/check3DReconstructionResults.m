function check3DReconstructionResults(reprojectionErrors, FinalResult, RD_L, frameNum)
% check3DReconstructionResults: Checks 3D reconstruction results and displacements.
%
%   check3DReconstructionResults(reprojectionErrors, FinalResult, RD_L, frameNum)
%
%   Inputs:
%     reprojectionErrors - Cell array containing reprojection errors.
%     FinalResult        - Structure containing 3D coordinates and displacements.
%     RD_L               - Structure containing FEMesh information (for Plotdisp_show3D_STAQ).
%     frameNum           - (Optional) The frame number to view. Defaults to 2.

if nargin < 4
    frameNum = 2; % Default to viewing results for the second frame
end

close all; % Close all existing figure windows

% Check and display reprojection errors for the specified frame
if frameNum <= size(reprojectionErrors, 1)
    disp(['--- Frame ', num2str(frameNum), ' Results ---']);
    disp(['3D_reconstruction_error_mean = ', num2str(mean(reprojectionErrors{frameNum,1}(:)))]);
    disp(['3D_reconstruction_error_std = ', num2str(std(reprojectionErrors{frameNum,1}(:)))]);
else
    warning('Specified frameNum (%d) is out of bounds for reprojectionErrors. Skipping error display.', frameNum);
end

figure;
sgtitle(['3D Reconstruction Results for Frame ', num2str(frameNum)]); % Super title for the figure

% Reprojection error histogram for Frame 1
subplot(3,2,1);
% Always display error for Frame 1 for comparison
if ~isempty(reprojectionErrors{1,1})
    histogram(reprojectionErrors{1,1}, 100);
    title('Reprojection Error Frame 1');
else
    title('Reprojection Error Frame 1 (Data Missing)');
end

% Reprojection error histogram for the specified frame
subplot(3,2,2);
if frameNum <= size(reprojectionErrors, 1) && ~isempty(reprojectionErrors{frameNum,1})
    histogram(reprojectionErrors{frameNum,1}, 100);
    title(['Reprojection Error Frame ', num2str(frameNum)]);
else
    title(['Reprojection Error Frame ', num2str(frameNum), ' (Data Missing)']);
end


% 3D Coordinates scatter plot for Frame 1
subplot(3,2,3);
if ~isempty(FinalResult.Coordinates{1,1})
    scatter3(FinalResult.Coordinates{1,1}, FinalResult.Coordinates{1,2}, FinalResult.Coordinates{1,3}, 'filled');
    title('Coordinates at Frame 1');
    xlabel('X'); ylabel('Y'); zlabel('Z');
    axis equal; % Ensure equal scaling for axes
else
    title('Coordinates at Frame 1 (Data Missing)');
end

% 3D Coordinates scatter plot for the specified frame
subplot(3,2,4);
if frameNum <= size(FinalResult.Coordinates, 1) && ~isempty(FinalResult.Coordinates{frameNum,1})
    scatter3(FinalResult.Coordinates{frameNum,1}, FinalResult.Coordinates{frameNum,2}, FinalResult.Coordinates{frameNum,3}, 'filled');
    title(['Coordinates at Frame ', num2str(frameNum)]);
    xlabel('X'); ylabel('Y'); zlabel('Z');
    axis equal; % Ensure equal scaling for axes
else
    title(['Coordinates at Frame ', num2str(frameNum), ' (Data Missing)']);
end


% Displacement vector field (from Frame 1 to frameNum)
subplot(3,2,5);
if frameNum <= size(FinalResult.Coordinates, 1) && frameNum <= size(FinalResult.Displacement, 1) && ...
        ~isempty(FinalResult.Coordinates{1,1}) && ~isempty(FinalResult.Displacement{frameNum,1})
    quiver3(FinalResult.Coordinates{1,1}, FinalResult.Coordinates{1,2}, FinalResult.Coordinates{1,3}, ...
        FinalResult.Displacement{frameNum,1}, FinalResult.Displacement{frameNum,2}, FinalResult.Displacement{frameNum,3});
    title(['Tracked Displacement Vector (Frame 1 to ', num2str(frameNum), ')']);
    xlabel('X'); ylabel('Y'); zlabel('Z');
    axis equal; % Ensure equal scaling for axes
else
    title(['Tracked Displacement Vector (Frame 1 to ', num2str(frameNum), ') (Data Missing)']);
end


% Displacement components scatter plot (from Frame 1 to frameNum)
subplot(3,2,6);
if frameNum <= size(FinalResult.Coordinates, 1) && ~isempty(FinalResult.Coordinates{1,1}) && ~isempty(FinalResult.Coordinates{frameNum,1})
    scatter3(FinalResult.Coordinates{frameNum,1} - FinalResult.Coordinates{1,1}, ...
        FinalResult.Coordinates{frameNum,2} - FinalResult.Coordinates{1,2}, ...
        FinalResult.Coordinates{frameNum,3} - FinalResult.Coordinates{1,3}, 'filled');
    title(['Tracked Displacement Components (Frame 1 to ', num2str(frameNum), ')']);
    xlabel('dX'); ylabel('dY'); zlabel('dZ');
    axis equal; % Ensure equal scaling for axes
else
    title(['Tracked Displacement Components (Frame 1 to ', num2str(frameNum), ') (Data Missing)']);
end


% Display and check the deformed image results using Plotdisp_show3D_STAQ
if frameNum <= size(FinalResult.Displacement, 1) && ~isempty(FinalResult.Displacement{frameNum,1}) && ...
        isfield(RD_L, 'ResultFEMeshEachFrame') && ~isempty(RD_L.ResultFEMeshEachFrame) && ...
        ~isempty(RD_L.ResultFEMeshEachFrame{1,1}.coordinatesFEM) && ~isempty(RD_L.ResultFEMeshEachFrame{1,1}.elementsFEM)

    U3Dstereo = [FinalResult.Displacement{frameNum,1}, FinalResult.Displacement{frameNum,2}, FinalResult.Displacement{frameNum,3}]';
    % Magnitude = sqrt(U3Dstereo(1,:).^2 + U3Dstereo(2,:).^2 + U3Dstereo(3,:).^2); % Can be calculated but not directly used for plotting here

    % Check if Plotdisp_show3D_STAQ function exists
    Plotdisp_show3D_STAQ(U3Dstereo(:), RD_L.ResultFEMeshEachFrame{1,1}.coordinatesFEM, ...
        RD_L.ResultFEMeshEachFrame{1,1}.elementsFEM, 'NoEdgeColor');
    title(['3D Displacements Visualized for Frame ', num2str(frameNum)]);

else
    warning('Cannot visualize 3D displacements: Required FinalResult.Displacement or RD_L data is missing or empty.');
end

end