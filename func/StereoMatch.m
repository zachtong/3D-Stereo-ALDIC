function [StereoInfo, RD_L, RD_R] = StereoMatch(Df_L,RD_L,RD_R,fNormalized_L,fNormalized_R,file_name_L,DICpara,StereoInfo) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 
%   Output:     StereoInfo.ResultFEMesh_corr: Corresponding coordinates of subsets in the first right image 
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
showImgOrNot = 1; 
[DICpara,x0temp,y0temp,u,v,cc] = IntegerSearch(fNormalized_L,fNormalized_R,file_name_L,DICpara,showImgOrNot,-1); 
 
% ====== FEM mesh set up ====== 
[DICmesh] = MeshSetUp(x0temp,y0temp,DICpara); clear x0temp y0temp; 
 
% ====== Initial Value ====== 
U0 = Init(u,v,cc.max,DICmesh.x0,DICmesh.y0,0);    
 
% ====== Deal with incremental mode ====== 
% RD_L.ResultFEMesh{1+floor(1/DICpara.ImgSeqIncUnit)} = ... % To save first mesh info 
%     struct( 'coordinatesFEM',DICmesh.coordinatesFEM,'elementsFEM',DICmesh.elementsFEM, ... 
%     'winsize',DICpara.winsize,'winstepsize',DICpara.winstepsize,'gridxyROIRange',DICpara.gridxyROIRange ); 
 
RD_L.ResultFEMesh{1+ floor(0.9/DICpara.ImgSeqIncUnit)} = ... % To save first mesh info 
    struct( 'coordinatesFEM',DICmesh.coordinatesFEM,'elementsFEM',DICmesh.elementsFEM, ... 
    'winsize',DICpara.winsize,'winstepsize',DICpara.winstepsize,'gridxyROIRange',DICpara.gridxyROIRange ); 
 
tol = 1e-3; 
coordinatesFEM = RD_L.ResultFEMesh{1, 1}.coordinatesFEM; 
 
% ====== 1st-order ICGN refinement  ====== 
[U,F,~,~,~,~] = LocalICGN(U0,coordinatesFEM,Df_L,... 
    fNormalized_L,fNormalized_R,DICpara,'GaussNewton',tol,1); 
 
% ====== 2nd-order ICGN refinement  ====== 
% [U,~,~,~,~,~] = LocalICGN(U0,coordinatesFEM,Df_L,... 
%     fNormalized_L,fNormalized_R,DICpara,'GaussNewton',tol,2); 
 
% ------ Plot ------ 
close all; % Plotuv(U,DICmesh.x0,DICmesh.y0); 
Plotdisp_show(U,DICmesh.coordinatesFEM,DICmesh.elementsFEM); 
% Plotdisp_show(U0,DICmesh.coordinatesFEM,DICmesh.elementsFEM); 
% Plotstrain_show(F,DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM,[],'NoEdgeColor');

 
% Save data 
StereoInfo.U=U; 
%StereoInfo.F=F; 
%StereoInfo.HtempPar=HtempPar; 
%StereoInfo.LocalTime=LocalTime; 
%StereoInfo.ConvItPerEle=ConvItPerEle; 
%StereoInfo.LocalICGNBadPtNum=LocalICGNBadPtNum; 
 
 
% Determine the corresponding ROI in the first right image 
Utemp = reshape(StereoInfo.U,2,size(StereoInfo.U,1)/2)'; 
RD_R.ResultFEMesh_corr{1,1} = RD_L.ResultFEMesh{1, 1}.coordinatesFEM + Utemp; 
StereoInfo.ResultFEMesh_corr =  RD_R.ResultFEMesh_corr{1,1} ;
RD_R.CalPointsNum = size(u); 
 
end 
 
