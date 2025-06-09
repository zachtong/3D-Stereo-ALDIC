function [StereoInfo, RD_L, RD_R] = StereoMatchQuadtree(Df_L,RD_L,RD_R,fNormalized_L,fNormalized_R,file_name_L,...
    DICpara,StereoInfo,fNormalizedMask_L,fNormalizedMask_R)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Output:     StereoInfo.ResultFEMesh_corr: Corresponding coordinates of subsets in the first right image
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
showImgOrNot = 1;
[DICpara,x0temp_f,y0temp_f,u_f,v_f,cc]= IntegerSearch(fNormalized_L,fNormalized_R,file_name_L,DICpara,showImgOrNot,-1);


%%%%% Interpolate to f %%%%%
xnodes = max([4+0.5*DICpara.winsize,DICpara.gridxyROIRange.gridx(1)])  ...
    : DICpara.winstepsize : min([size(fNormalized_L,1)-0.5*DICpara.winsize-3,DICpara.gridxyROIRange.gridx(2)]);
ynodes = max([4+0.5*DICpara.winsize,DICpara.gridxyROIRange.gridy(1)])  ...
    : DICpara.winstepsize : min([size(fNormalized_L,2)-0.5*DICpara.winsize-3,DICpara.gridxyROIRange.gridy(2)]);

[x0temp,y0temp] = ndgrid(xnodes,ynodes);
u_f_NotNanInd = find(~isnan(u_f(:)));

op1 = rbfcreate( [x0temp_f(u_f_NotNanInd),y0temp_f(u_f_NotNanInd)]',[u_f(u_f_NotNanInd)]','RBFFunction', 'thinplate'); rbfcheck(op1);
u = rbfinterp( [x0temp(:),y0temp(:)]', op1 );
op2 = rbfcreate( [x0temp_f(u_f_NotNanInd),y0temp_f(u_f_NotNanInd)]',[v_f(u_f_NotNanInd)]','RBFFunction', 'thinplate'); rbfcheck(op2);
v = rbfinterp([x0temp(:),y0temp(:)]', op2 );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure, plot3(x0temp(:),y0temp(:),u(:),'.');
% figure, plot3(x0temp(:),y0temp(:),v(:),'.');
%%%%% Do some regularization to further decrease the noise %%%%%
u = regularizeNd([x0temp(:),y0temp(:)],u(:),{xnodes',ynodes'},1e-3);
v = regularizeNd([x0temp(:),y0temp(:)],v(:),{xnodes',ynodes'},1e-3);


% ====== FEM mesh set up ======
[DICmesh] = MeshSetUp(x0temp,y0temp,DICpara); % clear x0temp y0temp;

% ====== Initial Value ======
U0 = Init(u,v,cc.max,DICmesh.x0,DICmesh.y0,0);   

for tempi = 1:size(u,1)
    for tempj = 1:size(u,2)
        try
            if ~fNormalizedMask_L(x0temp(tempi,tempj),y0temp(tempi,tempj)) || ...
                    ~fNormalizedMask_R( floor(x0temp(tempi,tempj)+u(tempi,tempj)), floor(y0temp(tempi,tempj)+v(tempi,tempj)) )

                U0(2*(tempj+(tempi-1)*(size(u,2)))) = nan;
                U0(2*(tempj+(tempi-1)*(size(u,2)))-1) = nan;

            end
        catch
        end
    end
end
        

% ====== Deal with incremental mode ======
% RD_L.ResultFEMesh{1+floor(1/DICpara.ImgSeqIncUnit)} = ... % To save first mesh info
%     struct( 'coordinatesFEM',DICmesh.coordinatesFEM,'elementsFEM',DICmesh.elementsFEM, ...
%     'winsize',DICpara.winsize,'winstepsize',DICpara.winstepsize,'gridxyROIRange',DICpara.gridxyROIRange );

RD_L.ResultFEMesh{1+ floor(0.9/DICpara.ImgSeqIncUnit)} = ... % To save first mesh info
    struct( 'coordinatesFEM',DICmesh.coordinatesFEM,'elementsFEM',DICmesh.elementsFEM, ...
    'winsize',DICpara.winsize,'winstepsize',DICpara.winstepsize,'gridxyROIRange',DICpara.gridxyROIRange );


 


% ====== Generate a quadtree mesh considering sample's complex geometry ======
DICmesh.elementMinSize = 4; % min element size in the refined quadtree mesh

GenerateQuadtreeMeshStereo; % Generate a quadtree mesh


% ICGN stopping tolerance
tol = 1e-3;

% ====== 1st-order ICGN refinement  ======
%[U,F,HtempPar,LocalTime,ConvItPerEle,LocalICGNBadPtNum] = LocalICGN(U0,coordinatesFEM,Df_L,...
    %fNormalized_L,fNormalized_R,DICpara,'GaussNewton',tol,1);

% ====== 2nd-order ICGN refinement  ======
[U,~,~,~,~,~] = LocalICGNQuadtree(U0,DICmesh.coordinatesFEM,Df_L,...
    fNormalized_L,fNormalized_R,DICpara,'GaussNewton',tol,2);

% ------ Plot ------
close all; % Plotuv(U,DICmesh.x0,DICmesh.y0);
Plotdisp_show(U,DICmesh.coordinatesFEM,DICmesh.elementsFEM);
Plotdisp_show(U0,DICmesh.coordinatesFEM,DICmesh.elementsFEM);


% Save data
StereoInfo.U = U;
%StereoInfo.F=F;
%StereoInfo.HtempPar=HtempPar;
%StereoInfo.LocalTime=LocalTime;
%StereoInfo.ConvItPerEle=ConvItPerEle;
%StereoInfo.LocalICGNBadPtNum=LocalICGNBadPtNum;


% Determine the corresponding ROI in the first right image
Utemp = reshape(StereoInfo.U,2,size(StereoInfo.U,1)/2)';
RD_R.ResultFEMesh_corr{1,1} = DICmesh.coordinatesFEM + Utemp;
% RD_R.CalPointsNum = size(u_f);




end

