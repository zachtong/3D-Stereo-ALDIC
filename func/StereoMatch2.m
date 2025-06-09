function [StereoInfo, RD_L, RD_R] = StereoMatch2(Df,RD_L,RD_R,fNormalized_L,fNormalized_R,DICpara,StereoInfo)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Output: RD_R.corresponding_cooridinate:
%   
%   Input:  StereoInfo: 
%           RD_L:
%           RD_R:
%           DICpara:
%           searchStep:
%           searchRadius:
%           disparity:           
%
%   Methods:    EpipolarSearch:
%               ICGN1:
%               ICGN2:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialization
searchStep = input("searchStep:");
searchRadius = input("searchRadius:");
parallax.X = input("parallax.X:");
parallax.Y = input("parallax.Y:");
tol = 1e-2;
ICGNmethod = 'GaussNewton'; % Default

% Mesh grid for Left camera
[X_coords, Y_coords] = generateGridMatrices(DICpara);
meshGrid = MeshSetUp(X_coords,Y_coords,DICpara); % n by 2

% Calculate initial guess for every subset
subsetsNum = size(meshGrid.coordinatesFEM,1);
ZNSSDList = cell(subsetsNum,1);
pointCandidate = cell(subsetsNum,1);

for i = 1: subsetsNum
    % Obtain potiential corresponding points
    pointCandidate{i} = FindCandidate(meshGrid.coordinatesFEM(i,:),...
        parallax,searchStep,searchRadius,DICpara,StereoInfo); % pointCandidate{i,1} is 2 by n (X,Y coordinates)
    % Coarse finding the initial guess
    for j = 1:size(pointCandidate{i},2)
        ZNSSDList{i}(j) = epipolar_ICGN1(pointCandidate{i}(:,j),Df,fNormalized_L,fNormalized_R,DICpara.winsize,tol,ICGNmethod); % ZNSSDList{i,1} is 1 by n
    end
end

% Check initial guess (find mininum)
initialGuess = determineInitialGuess(pointCandidate,ZNSSDList);

% Do ICGN2 using ICGN1 results
for i = 1:subsetsNum
    [U,V] = epipolar_ICGN2(initialGuess);
    
end

% Save data (TBC)
U;
V;

end
    

function pointCandidate = FindCandidate(meshGrid,parallax,searchStep,searchRadius,DICpara,StereoInfo)  
% Determine epipolar line
epipolar = StereoInfo.cameraParams.FundamentalMatrix * [meshGrid , 1]';
line_slope = -epipolar(1) / epipolar(2);
line_intercept = -epipolar(3) / epipolar(2);

% Find searching center in epipolar line
center.X = round(((line_slope * (meshGrid(2) + parallax.Y - line_intercept)...
	+ meshGrid(1) + parallax.X) / (line_slope * line_slope + 1)));
center.Y = round(line_slope * center.X + line_intercept);

% Get Candidate points
i = 1;
flag = 0;
while (flag < searchRadius)
    flag = (i-1) * searchStep;
    % Left hand side
    X_trial = center.X + flag;
    Y_trial = round(line_slope * X_trial + line_intercept);
    % Judge whether [X_trial,Y_trial] is legal
    if (X_trial-DICpara.winsize/2 > 1) && (X_trial+DICpara.winsize/2<DICpara.ImgSize(1)) ...
           && (Y_trial-DICpara.winsize/2 > 1) && (Y_trial+DICpara.winsize/2<DICpara.ImgSize(2))
        pointCandidate(1,i) = X_trial; pointCandidate(2,i) = Y_trial;
    end
    % Right hand side
    X_trial = center.X - flag;
    Y_trial = round(line_slope * X_trial + line_intercept);
    % Judge whether [X_trial,Y_trial] is legal
    if (X_trial-DICpara.winsize/2 > 1) && (X_trial+DICpara.winsize/2<DICpara.ImgSize(1)) ...
           && (Y_trial-DICpara.winsize/2 > 1) && (Y_trial+DICpara.winsize/2<DICpara.ImgSize(2))
        pointCandidate(3,i) = X_trial; pointCandidate(4,i) = Y_trial;
    end
    i = i+1;
end

pointCandidate = [pointCandidate(1:2,:),pointCandidate(3:4,:)];

% Checking illegal numbers in pointCandidate
% TBD

end


% Determine initial guess by the ZNSSD (Find the minimum)
function initialGuess = determineInitialGuess(pointCandidate,ZNSSDList)

end

