function [StereoInfo, RD_L, RD_R] = StereoMatch_STAQ(RD_L,RD_R,Normalized_L,Normalized_R,file_name_L,NormalizedMask_L,NormalizedMask_R,DICpara,StereoInfo,stereoMatchShapeOrder)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Output:     StereoInfo.ResultFEMesh_corr: Corresponding coordinates of subsets in the first right image
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DICpara.showImgOrNot = 1;
incOrNot = DICpara.DICIncOrNot; % 0: acc, 1: inc
Df = funImgGradient(Normalized_L,Normalized_L,NormalizedMask_L); % Finite difference to compute image grayscale gradients;

% ====== Integer Search ======
% DICpara.InitFFTSearchMethod = 1;
[DICpara,x0temp,y0temp,u,v,cc]= IntegerSearchQuadtree(Normalized_L,Normalized_R,file_name_L,DICpara,-1);

% ====== FEM mesh set up ======
[DICmesh] = MeshSetUp(x0temp,y0temp,DICpara); clear x0temp y0temp;
% ====== Initial Value ======
U0 = Init(u,v,cc.max,DICmesh.x0,DICmesh.y0,0);

% Zach Modified
% Set zero at holes
linearIndices1 = sub2ind(size(NormalizedMask_L), DICmesh.coordinatesFEM(:,1), DICmesh.coordinatesFEM(:,2));
MaskOrNot1 = NormalizedMask_L(linearIndices1);
% Do we need to consider the mask of g? TBD
% u_inv = u'; v_inv = v';
% linearIndices2 = sub2ind(size(NormalizedMask_R), floor(DICmesh.coordinatesFEM(:,1)+u_inv(:)), floor(DICmesh.coordinatesFEM(:,2)+v_inv(:)));
% MaskOrNot2 = NormalizedMask_R(linearIndices2);
% MaskOrNot = MaskOrNot1 + MaskOrNot2;
% nanIndex = find(MaskOrNot<2);

nanIndex = find(MaskOrNot1<1);
U0(2*nanIndex) = nan;
U0(2*nanIndex-1) = nan;

% ====== Deal with incremental mode ======
% RD_L.ResultFEMesh{1+floor(1/DICpara.ImgSeqIncUnit)} = ... % To save first mesh info
%     struct( 'coordinatesFEM',DICmesh.coordinatesFEM,'elementsFEM',DICmesh.elementsFEM, ...
%     'winsize',DICpara.winsize,'winstepsize',DICpara.winstepsize,'gridxyROIRange',DICpara.gridxyROIRange );

% To save first Non-quadtree mesh info
StereoInfo.ResultFEMesh = struct( 'coordinatesFEM',DICmesh.coordinatesFEM,'elementsFEM',DICmesh.elementsFEM, ...
    'winsize',DICpara.winsize,'winstepsize',DICpara.winstepsize,'gridxyROIRange',DICpara.gridxyROIRange );
% ====== Generate a quadtree mesh considering sample's complex geometry ======
DICmesh.elementMinSize = DICpara.winsizeMin; % min element size in the refined quadtree mesh

% Notes:
% Hanging nodes and sub-elements are placed on the last
% All the void regions are generating nodes but we can ignore them
% using maskfile later.
[DICmesh,DICpara,U0] = GenerateQuadtreeMesh(U0,Df,NormalizedMask_L,DICmesh,DICpara); % Generate the quadtree mesh

% ====== Store current mesh ======
StereoInfo.ResultFEMeshEachFrame = struct( 'coordinatesFEM',...
    DICmesh.coordinatesFEM,'elementsFEM',DICmesh.elementsFEM,...
    'markCoordHoleEdge',DICmesh.markCoordHoleEdge );

%% ICGN
tol = 1e-3;

if stereoMatchShapeOrder == 1
    % ====== 1st-order ICGN refinement  ======
    [U,~,~,~,~,~,~] = LocalICGNQuadtree(U0,DICmesh.coordinatesFEM,Df,...
        Normalized_L,Normalized_R,DICpara,'GaussNewton',tol,1);
else
    % ====== 2nd-order ICGN refinement  ======
    [U,~,~,~,~,~,~] = LocalICGNQuadtree(U0,DICmesh.coordinatesFEM,Df,...
        Normalized_L,Normalized_R,DICpara,'GaussNewton',tol,2);
end

% ------ Plot ------
close all; % Plotuv(U,DICmesh.x0,DICmesh.y0);
Plotdisp_show(U,DICmesh.coordinatesFEM,DICmesh.elementsFEM);
%Plotdisp_show(U0,DICmesh.coordinatesFEM,DICmesh.elementsFEM);


% Save data
StereoInfo.U=U;
%StereoInfo.F=F;
%StereoInfo.HtempPar=HtempPar;
%StereoInfo.LocalTime=LocalTime;
%StereoInfo.ConvItPerEle=ConvItPerEle;
%StereoInfo.LocalICGNBadPtNum=LocalICGNBadPtNum;


% Determine the corresponding ROI in the first right image
Utemp = reshape(StereoInfo.U,2,size(StereoInfo.U,1)/2)';
StereoInfo.ResultFEMesh_corr = StereoInfo.ResultFEMeshEachFrame.coordinatesFEM + Utemp;  % If strategy 2, this variable could be modify as cells
RD_R.Coordinates_corr = StereoInfo.ResultFEMesh_corr;

StereoInfo.CalPointsNum = size(u);


end

