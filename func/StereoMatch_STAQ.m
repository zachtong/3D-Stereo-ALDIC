function [StereoInfo, RD_L, RD_R] = StereoMatch_STAQ(RD_L,RD_R,Normalized_L,Normalized_R,file_name_L,NormalizedMask_L,NormalizedMask_R,DICpara,StereoInfo,stereoMatchShapeOrder)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Output:     StereoInfo.ResultFEMesh_corr: Corresponding coordinates of subsets in the first right image
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Respect caller's showImgOrNot setting (don't override to 1)
if ~isfield(DICpara,'showImgOrNot'), DICpara.showImgOrNot = 1; end
incOrNot = DICpara.DICIncOrNot; % 0: acc, 1: inc
Df = funImgGradient(Normalized_L,Normalized_L,NormalizedMask_L); % Finite difference to compute image grayscale gradients;

% ====== Integer Search ======
% DICpara.InitFFTSearchMethod = 1;
[DICpara,x0temp,y0temp,u,v,cc]= IntegerSearchQuadtree(Normalized_L,Normalized_R,file_name_L,DICpara,-1);

% ====== FEM mesh set up ======
[DICmesh] = MeshSetUp(x0temp,y0temp,DICpara); clear x0temp y0temp;
% ====== Initial Value ======
U0 = Init(u,v,cc.max,DICmesh.x0,DICmesh.y0,0);

% Set NaN at holes using mask
linearIndices1 = sub2ind(size(NormalizedMask_L), round(DICmesh.coordinatesFEM(:,1)), round(DICmesh.coordinatesFEM(:,2)));
MaskOrNot1 = NormalizedMask_L(linearIndices1);

nanIndex = find(MaskOrNot1<1);
U0(2*nanIndex) = nan;
U0(2*nanIndex-1) = nan;

nTotal = size(DICmesh.coordinatesFEM,1);
nMasked = length(nanIndex);
nValid = nTotal - nMasked;
fprintf('  Mask: %d/%d valid (%.0f%%)\n', nValid, nTotal, 100*nValid/nTotal);

if nValid < 10
    warning('StereoMatch: only %d valid points after masking! Check mask orientation.', nValid);
    % Debug: sample mask values at a few coordinate locations
    sampleIdx = round(linspace(1, nTotal, min(10,nTotal)));
    for si = 1:length(sampleIdx)
        idx = sampleIdx(si);
        cx = round(DICmesh.coordinatesFEM(idx,1));
        cy = round(DICmesh.coordinatesFEM(idx,2));
        fprintf('  coord(%d)=[%d,%d] -> mask=%d\n', idx, cx, cy, NormalizedMask_L(cx,cy));
    end
end

% Store first non-quadtree mesh info
StereoInfo.ResultFEMesh = struct( 'coordinatesFEM',DICmesh.coordinatesFEM,'elementsFEM',DICmesh.elementsFEM, ...
    'winsize',DICpara.winsize,'winstepsize',DICpara.winstepsize,'gridxyROIRange',DICpara.gridxyROIRange );
% ====== Generate a quadtree mesh considering sample's complex geometry ======
DICmesh.elementMinSize = DICpara.winsizeMin; % min element size in the refined quadtree mesh

% Notes:
% Hanging nodes and sub-elements are placed on the last
% All the void regions are generating nodes but we can ignore them
% using maskfile later.
% Quadtree mesh generation

[DICmesh,DICpara,U0] = GenerateQuadtreeMesh(U0,Df,NormalizedMask_L,DICmesh,DICpara); % Generate the quadtree mesh

% ====== Store current mesh ======
StereoInfo.ResultFEMeshEachFrame = struct( 'coordinatesFEM',...
    DICmesh.coordinatesFEM,'elementsFEM',DICmesh.elementsFEM,...
    'markCoordHoleEdge',DICmesh.markCoordHoleEdge );

%% ICGN
nQuadPts = size(DICmesh.coordinatesFEM, 1);
nQuadNaN = nnz(isnan(U0(1:2:end)));
fprintf('  Quadtree: %d pts, %.0f%% valid\n', nQuadPts, 100*(nQuadPts-nQuadNaN)/nQuadPts);

tol = 1e-3;

if stereoMatchShapeOrder == 1
    % ====== 1st-order ICGN refinement  ======
    [U,~,~,~,~,~,~] = LocalICGNQuadtree(U0,DICmesh.coordinatesFEM,Df,...
        Normalized_L,Normalized_R,DICpara,'GaussNewton',tol,1);
else
    % ====== 2nd-order ICGN refinement (not implemented — reserved) ======
    [U,~,~,~,~,~,~] = LocalICGNQuadtree(U0,DICmesh.coordinatesFEM,Df,...
        Normalized_L,Normalized_R,DICpara,'GaussNewton',tol,2);
end

% Plot (only if showImgOrNot=1)
if nnz(~isnan(U)) > 0 && DICpara.showImgOrNot == 1
    Plotdisp_show(U,DICmesh.coordinatesFEM,DICmesh.elementsFEM);
end
%Plotdisp_show(U0,DICmesh.coordinatesFEM,DICmesh.elementsFEM);


% Save data
StereoInfo.U=U;
%StereoInfo.F=F;
%StereoInfo.HtempPar=HtempPar;
%StereoInfo.LocalTime=LocalTime;
%StereoInfo.ConvItPerEle=ConvItPerEle;
%StereoInfo.LocalICGNBadPtNum=LocalICGNBadPtNum;


% Determine the corresponding ROI in the first right image
nCoords = size(StereoInfo.ResultFEMeshEachFrame.coordinatesFEM, 1);
if numel(StereoInfo.U) == 2*nCoords
    Utemp = reshape(StereoInfo.U, 2, nCoords)';
    StereoInfo.ResultFEMesh_corr = StereoInfo.ResultFEMeshEachFrame.coordinatesFEM + Utemp;
else
    warning('StereoMatch: U size (%d) != 2*coordinatesFEM (%d). Using zero disparity.', ...
        numel(StereoInfo.U), 2*nCoords);
    Utemp = zeros(nCoords, 2);
    StereoInfo.ResultFEMesh_corr = StereoInfo.ResultFEMeshEachFrame.coordinatesFEM + Utemp;
end
RD_R.Coordinates_corr = StereoInfo.ResultFEMesh_corr;

StereoInfo.CalPointsNum = size(u);


end

