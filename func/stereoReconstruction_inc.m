function [reconstructedPoints,reprojectionErrors] = stereoReconstruction_inc(cameraParams, RD_L, RD_R)
% Inputs:
% K_left, K_right: Intrinsic matrices of the left and right cameras.
% R_left, T_left, R_right, T_right: Extrinsic parameters (rotation matrices and translation vectors) of the left and right cameras.
% matchedPairs: Matched pixel pairs between the left and right images.

% Output:
% reconstructedPoints: 3D reconstructed points.

% Initializition
K_left = cameraParams.cameraParamsLeft.K;
K_right = cameraParams.cameraParamsRight.K;
R_left = [1 0 0; 0 1 0; 0 0 1]; % Set left Camera Coordinate as World Coordinate
T_left = [0 0 0]'; % Set left Camera Coordinate as World Coordinate
R_right = cameraParams.rotationMatrix;
T_right = cameraParams.translationVector';

% Compute the fundamental and essential matrices.
% [F, E] = computeFundamentalEssentialMatrix(K_left, K_right, R_left, T_left, R_right, T_right);

% Triangulate to reconstruct 3D points for each frame.
% Initializition
imageNum = size(RD_L.ResultDisp,1); % numbers except the first image
matchedPairs = cell(imageNum,1);

% task 1 2 3 4 

%%%% Since second frame %%%%%
% for i = 1:imageNum+1
% 
%     matchedPairs{i,1}(:,1:2) = RD_L.coor_Ch2{i} ;
%     matchedPairs{i,1}(:,3:4) = RD_R.coor_Ch2{i} ;
% 
% end

%% Original code
for i = 1:imageNum+1

    matchedPairs{i,1}(:,1:2) = RD_L.Coordinate_inc_ref{i}.Coor ;
    matchedPairs{i,1}(:,3:4) = RD_R.Coordinate_inc_ref{i}.Coor ;

end


% Undistort points before doing 3D reconstruction
[matchedPairs_undistort]= funUndistortPoints(matchedPairs,cameraParams);
[temp_reconstructedPoints,reprojectionErrors] = triangulatePoints(matchedPairs_undistort, K_left, K_right, R_left, T_left, R_right, T_right);
for j = 1:length(temp_reconstructedPoints)
    reconstructedPoints{j,1} = temp_reconstructedPoints{j}(:,1);
    reconstructedPoints{j,2} = temp_reconstructedPoints{j}(:,2);
    reconstructedPoints{j,3} = temp_reconstructedPoints{j}(:,3);
end
end

function [F, E] = computeFundamentalEssentialMatrix(K_left, K_right, R_left, T_left, R_right, T_right)
% Compute the essential matrix.
E = K_right' * skewSymmetric(T_right - T_left) * R_right' * R_left * K_left;
% Compute the fundamental matrix.
F = inv(K_right)' * E * inv(K_left);
end


function [points3D,reprojectionErrors] = triangulatePoints(matchedPairs, K_left, K_right, R_left, T_left, R_right, T_right)
% Purpose: Use the matched pixel pairs, fundamental matrix, and camera parameters to triangulate and reconstruct 3D points.
P_left = K_left * [R_left, T_left];
P_right = K_right * [R_right, T_right];
points3D = cell(size(matchedPairs,1),1);
reprojectionErrors = cell(size(matchedPairs,1),1);
    for i = 1:size(matchedPairs,1)
        [points3D{i,1},reprojectionErrors{i}]= triangulate(matchedPairs{i,1}(:, 1:2), matchedPairs{i,1}(:, 3:4), P_left, P_right);
    end

end


function S = skewSymmetric(v)
S = [0, -v(3), v(2); v(3), 0, -v(1); -v(2), v(1), 0];
end
