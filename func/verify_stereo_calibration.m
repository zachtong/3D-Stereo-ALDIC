function [reconstructedPoints, reprojectionErrors] = verify_stereo_calibration(StereoInfo, showPlot)
%VERIFY_STEREO_CALIBRATION  Sanity check for stereo calibration + first-frame triangulation.
%
%   [pts3D, errs] = verify_stereo_calibration(StereoInfo)          % plots + prints
%   [pts3D, errs] = verify_stereo_calibration(StereoInfo, false)   % no plots
%
%   Uses the left/right matched points stored on StereoInfo (produced by
%   StereoMatch_STAQ) to triangulate frame-1 3D coordinates, and reports
%   the mean reprojection error. Useful for confirming that the
%   calibration imported from MATLAB/MatchID/MCC/DICe/OpenCorr is
%   numerically reasonable before running the full temporal pipeline.
%
%   Extracted from inline block at the bottom of main_3D_ALDIC.m
%   Section 3.2. Gated by DICpara.verifyStereoReconstruction.

if nargin < 2, showPlot = true; end

RD0_L_Pts = StereoInfo.ResultFEMeshEachFrame.coordinatesFEM;
RD0_R_Pts = StereoInfo.ResultFEMesh_corr;
matchedPairs = { [RD0_L_Pts, RD0_R_Pts] };

cameraParams = StereoInfo.cameraParams;
K_left  = cameraParams.cameraParamsLeft.K;
K_right = cameraParams.cameraParamsRight.K;
R_left  = eye(3);
T_left  = zeros(3, 1);
R_right = cameraParams.rotationMatrix;
T_right = cameraParams.translationVector';

matchedPairs_undistort = funUndistortPoints(matchedPairs, cameraParams);

P_left  = K_left  * [R_left,  T_left];
P_right = K_right * [R_right, T_right];

reconstructedPoints = cell(size(matchedPairs_undistort, 1), 1);
reprojectionErrors  = cell(size(matchedPairs_undistort, 1), 1);

for i = 1:size(matchedPairs_undistort, 1)
    [reconstructedPoints{i,1}, reprojectionErrors{i}] = triangulate( ...
        matchedPairs_undistort{i,1}(:, 1:2), ...
        matchedPairs_undistort{i,1}(:, 3:4), ...
        P_left, P_right);
end

if showPlot
    figure; plot(RD0_L_Pts(:,1), RD0_L_Pts(:,2), 'o');
    hold on; plot(RD0_R_Pts(:,1), RD0_R_Pts(:,2), 'o');
    title('Matched 2D points: left (blue) + right-corresponding (red)');

    figure; scatter3(reconstructedPoints{1,1}(:,1), ...
                     reconstructedPoints{1,1}(:,2), ...
                     reconstructedPoints{1,1}(:,3), 5, 'filled');
    title('Reconstructed 3D points (frame 1)');

    fprintf('Mean reprojection error: %.4f px\n', mean(reprojectionErrors{1}, 'omitnan'));
end
end
