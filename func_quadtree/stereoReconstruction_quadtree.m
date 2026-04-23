function [FinalResult,reprojectionErrors] = stereoReconstruction_quadtree(matchedPairs, cameraParams)
% Initializition
K_left = cameraParams.cameraParamsLeft.K;
K_right = cameraParams.cameraParamsRight.K;
R_left = [1 0 0; 0 1 0; 0 0 1]; % Set left Camera Coordinate as World Coordinate
T_left = [0 0 0]'; % Set left Camera Coordinate as World Coordinate
R_right = cameraParams.rotationMatrix;
T_right = cameraParams.translationVector';

%% Triangulate to reconstruct 3D points for each frame.
% Undistort points before doing 3D reconstruction
[matchedPairs_undistort]= funUndistortPoints(matchedPairs,cameraParams);

[temp_reconstructedPoints,reprojectionErrors] = triangulatePoints(matchedPairs_undistort, K_left, K_right, R_left, T_left, R_right, T_right);

for j = 1:length(temp_reconstructedPoints)
    FinalResult.Coordinates{j,1} = temp_reconstructedPoints{j}(:,1);
    FinalResult.Coordinates{j,2} = temp_reconstructedPoints{j}(:,2);
    FinalResult.Coordinates{j,3} = temp_reconstructedPoints{j}(:,3);
end

% Calculate the 3D displacement
for i = 1:size(FinalResult.Coordinates,1)
    FinalResult.Displacement{i,1} = FinalResult.Coordinates{i,1} - FinalResult.Coordinates{1,1}; % U
    FinalResult.Displacement{i,2} = FinalResult.Coordinates{i,2} - FinalResult.Coordinates{1,2}; % V
    FinalResult.Displacement{i,3} = FinalResult.Coordinates{i,3} - FinalResult.Coordinates{1,3}; % W
end



end


function [points3D,reprojectionErrors] = triangulatePoints(matchedPairs, K_left, K_right, R_left, T_left, R_right, T_right)
% Purpose: Use the matched pixel pairs, fundamental matrix, and camera parameters to triangulate and reconstruct 3D points.
% Handles NaN points from unconverged ICGN subsets.
P_left = K_left * [R_left, T_left];
P_right = K_right * [R_right, T_right];
points3D = cell(size(matchedPairs,1),1);
reprojectionErrors = cell(size(matchedPairs,1),1);
    for i = 1:size(matchedPairs,1)
        pts = matchedPairs{i,1};
        nPts = size(pts, 1);
        finiteIdx = all(isfinite(pts), 2);

        % Initialize with NaN
        points3D{i,1} = NaN(nPts, 3);
        reprojectionErrors{i} = NaN(nPts, 1);

        if any(finiteIdx)
            [pts3D_finite, reprErr_finite] = triangulate( ...
                pts(finiteIdx, 1:2), pts(finiteIdx, 3:4), P_left, P_right);
            points3D{i,1}(finiteIdx, :) = pts3D_finite;
            reprojectionErrors{i}(finiteIdx) = reprErr_finite;
        end

        % Per-frame summary (silent)
    end

end
