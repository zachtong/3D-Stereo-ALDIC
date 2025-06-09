function [reconstructedPoints,reprojectionErrors] = stereoReconstruction(cameraParams, RD_L, RD_R)
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
matchedPairs = cell(imageNum+1,1);

%% Challenge2
%task 1
if 0
    toBeInterpPoints =  innerPoints(1);
    matchedPairs{1,1} = [toBeInterpPoints, RD_R.coor_corr_Ch2];


    %%%%% Since second frame %%%%%
    for i = 1:imageNum

        matchedPairs{i+1,1}(:,1:2) = toBeInterpPoints + RD_L.disp_Ch2{i} ;
        matchedPairs{i+1,1}(:,3:4) = RD_R.coor_corr_Ch2 + RD_R.disp_Ch2{i};

    end
end
%%
if 0
    % task 3
    toBeInterpPoints =  innerPoints(3);
    matchedPairs{1,1} = [toBeInterpPoints, RD_R.coor_corr_Ch2];


    %%%%% Since second frame %%%%%
    for i = 1:imageNum
        % reorganize the format of matchedpairs, because the formats of RD_L
        % and RD_R are not identical.
        % Column 1,2 are Left image pixel coordinates(X_L,Y_L).
        % Column 3,4 are Right image pixel coordinates(X_R,Y_R).
        % The corresponding coordinates in L and R are named differently!
        % L = RD_L.ResultFEMesh{1, 1}.coordinatesFEM
        % R = RD_R.correspondingPixelCoors
        matchedPairs{i+1,1}(:,1:2) = toBeInterpPoints + RD_L.disp_Ch2{i} ;
        matchedPairs{i+1,1}(:,3:4) = RD_R.coor_corr_Ch2 + RD_R.disp_Ch2{i};
        % Zach modified!!! to test diff. strategies
        % matchedPairs{i+1,1}(:,3:4) = RD_R.ResultFEMesh_corr2;

    end
end
%% Original code
if 1
    %--------------------------------------------
    %%%% First frame %%%%%%
    matchedPairs{1,1} = [RD_L.ResultFEMesh{1, 1}.coordinatesFEM, RD_R.ResultFEMesh_corr{1}];
    %%%%% Since second frame %%%%%
    for i = 1:imageNum
        matchedPairs{i+1,1}(:,1:2) = RD_L.ResultFEMesh{1,1}.coordinatesFEM + reshape(RD_L.ResultDisp{i, 1}.U',[2,size(RD_L.ResultDisp{i, 1}.U,1)/2])' ;
        matchedPairs{i+1,1}(:,3:4) = RD_R.ResultFEMesh_corr{1,1} + RD_R.ResultDisp_corr{1,i}.U;

    end
    %---------------------------------------------
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
% Sub-function to compute the fundamental and essential matrices.

% Inputs:
% K_left, K_right: Intrinsic matrices of the left and right cameras.
% R_left, T_left, R_right, T_right: Extrinsic parameters of the cameras.

% Outputs:
% F: Fundamental matrix.
% E: Essential matrix.

% Purpose: Calculate the fundamental and essential matrices using the provided camera parameters to provide accurate geometric relations for subsequent triangulation.

% Compute the essential matrix.
E = K_right' * skewSymmetric(T_right - T_left) * R_right' * R_left * K_left;

% Compute the fundamental matrix.
F = inv(K_right)' * E * inv(K_left);
end

function [points3D,reprojectionErrors] = triangulatePoints(matchedPairs, K_left, K_right, R_left, T_left, R_right, T_right)
% Sub-function for triangulation to reconstruct 3D points.

% Inputs:
% matchedPairs: Matched pixel pairs between the left and right images.
% F: Fundamental matrix.
% K_left, K_right: Intrinsic matrices of the left and right cameras.
% R_left, T_left, R_right, T_right: Extrinsic parameters of the cameras.

% Output:
% points3D: 3D reconstructed points.

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
% Helper function to compute the skew-symmetric matrix of a vector.

% Input:
% v: A 3x1 vector.

% Output:
% S: Skew-symmetric matrix of the vector v.

S = [0, -v(3), v(2); v(3), 0, -v(1); -v(2), v(1), 0];
end
