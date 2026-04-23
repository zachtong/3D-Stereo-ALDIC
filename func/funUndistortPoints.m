function [undistortedPoints]...
    = funUndistortPoints(matchedPairs, cameraParams)
% Init
ImgPairsNum = size(matchedPairs,1);
PointsNum = size(matchedPairs{1,1},1);
undistortedPoints = matchedPairs;

cameraParamsLeft = cameraParameters('K', cameraParams.cameraParamsLeft.K, 'RadialDistortion', cameraParams.cameraParamsLeft.RadialDistortion, 'TangentialDistortion', cameraParams.cameraParamsLeft.TangentialDistortion);
cameraParamsRight = cameraParameters('K', cameraParams.cameraParamsRight.K, 'RadialDistortion', cameraParams.cameraParamsRight.RadialDistortion,  'TangentialDistortion', cameraParams.cameraParamsRight.TangentialDistortion);

% If your matlab version is older than 2023b, use these two line
%cameraParamsLeft = cameraParameters('IntrinsicMatrix', cameraParams.cameraParamsLeft.K, 'RadialDistortion', cameraParams.cameraParamsLeft.RadialDistortion, 'TangentialDistortion', cameraParams.cameraParamsLeft.TangentialDistortion);
%cameraParamsRight = cameraParameters('IntrinsicMatrix', cameraParams.cameraParamsRight.K, 'RadialDistortion', cameraParams.cameraParamsRight.RadialDistortion,  'TangentialDistortion', cameraParams.cameraParamsRight.TangentialDistortion);

% Sliced calculation
% chunkSize = 10000;
% numChunks = ceil(PointsNum/ chunkSize);
% for i = 1:ImgPairsNum 
%     % Init temp. sliced. matchedpairs
%     for j = 1:numChunks-1
%         matchedPairs_sliced_L{j} = matchedPairs{i,1}(chunkSize*(j-1)+1:chunkSize*j,1:2);
%         matchedPairs_sliced_R{j} = matchedPairs{i,1}(chunkSize*(j-1)+1:chunkSize*j,3:4);
%     end
%     matchedPairs_sliced_L{numChunks} = matchedPairs{i,1}(chunkSize*(numChunks-1)+1:end,1:2);
%     matchedPairs_sliced_R{numChunks} = matchedPairs{i,1}(chunkSize*(numChunks-1)+1:end,3:4);
% 
%     % Batch-wise calculation
%     for j = 1:numChunks
%         temp_undistortedPoints_L = undistortPoints(matchedPairs_sliced_L{j}, cameraParamsLeft);
%         temp_undistortedPoints_R = undistortPoints(matchedPairs_sliced_R{j}, cameraParamsRight);
%         temp_undistortedPoints{j} = [temp_undistortedPoints_L, temp_undistortedPoints_R];
%     end
% 
%     % Assembly
%     undistortedPoints{i,1} = temp_undistortedPoints{1};
%     for j = 1:numChunks-1
%         undistortedPoints{i,1} = [undistortedPoints{i,1}; temp_undistortedPoints{j+1}];
%     end
% end

% Pointwise undistorting — skip NaN points (from unconverged ICGN subsets)
for i = 1:ImgPairsNum
    temp_matchedPairs = matchedPairs{i,1};
    nPts = size(temp_matchedPairs, 1);
    temp_undistortedPoints = NaN(nPts, 4);

    % Find finite points (all 4 coords must be finite)
    finiteIdx = find(all(isfinite(temp_matchedPairs), 2));
    nFinite = length(finiteIdx);

    if nFinite == 0
        warning('funUndistortPoints: frame %d has no finite points', i);
        undistortedPoints{i,1} = temp_undistortedPoints;
        continue;
    end

    % Extract only finite points for undistortion
    finitePairs = temp_matchedPairs(finiteIdx, :);
    temp_finite_undist = zeros(nFinite, 4);
    tic;
    parfor j = 1:nFinite
        temp_undistortedPoints_L = undistortPoints(finitePairs(j,1:2), cameraParamsLeft);
        temp_undistortedPoints_R = undistortPoints(finitePairs(j,3:4), cameraParamsRight);
        temp_finite_undist(j,:) = [temp_undistortedPoints_L, temp_undistortedPoints_R];
    end
    time(i) = toc;

    % Put results back, NaN points stay NaN
    temp_undistortedPoints(finiteIdx, :) = temp_finite_undist;
    undistortedPoints{i,1} = temp_undistortedPoints;
    % Summary printed by caller if needed
end

