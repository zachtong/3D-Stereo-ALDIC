function [RD] = TemporalMatch_quadtree_ST1(DICpara,file_name,ImgMask,ImgNormalized,RD,StereoInfo,camera0OrNot,shapeFuncOrder)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This part is designed for temporal DIC, which means dealing with left and
% right images separately using 2D-ALDIC.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% debug options
UseGlobal = 1; % Zach 20250620
DICpara.showImgOrNot = 1;

% Init
imageNum = length(ImgNormalized);
incOrNot = DICpara.DICIncOrNot; % 0: acc, 1: inc



if length(ImgNormalized) > 7
    DICpara.dataDrivenOrNot = funParaInput('DataDrivenOrNot');
    if DICpara.dataDrivenOrNot == 0
        ImgStartDataDrivenMode = 7;
    else
        ImgStartDataDrivenMode = length(ImgNormalized) + 1;
    end
else
    ImgStartDataDrivenMode = length(ImgNormalized) + 1;
end

%% For loop
for ImgSeqNum = 2 : imageNum
    disp(['Current image frame #: ', num2str(ImgSeqNum),'/',num2str(length(ImgNormalized))]);

    if incOrNot == 1
        fNormalizedMask = double(ImgMask{ImgSeqNum-1} ) ; % Load the mask file of previous frame
        gNormalizedMask = double(ImgMask{ImgSeqNum}); % Load the mask file of current frame
        fNormalized = ImgNormalized{ImgSeqNum-1}.* fNormalizedMask; % Load previous reference image frame
        gNormalized = ImgNormalized{ImgSeqNum}.* gNormalizedMask; % Load current deformed image frame
    else
        fNormalizedMask = double(ImgMask{1} ) ; % Load the mask file of previous frame
        gNormalizedMask = ones(size(ImgNormalized{1}));
        fNormalized = ImgNormalized{1}.* fNormalizedMask; % Load previous reference image frame
        gNormalized = ImgNormalized{ImgSeqNum}.* gNormalizedMask; % Load current deformed image frame
    end

    Df = funImgGradient(fNormalized,fNormalized,fNormalizedMask);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Zach Attention: because DICpara only stores left image info, to
    % consider the right image temporal matching, we "fake" update masks
    % and ROIrange (we didn't really change the DICpara in this func)
    DICpara.ImgRefMask = fNormalizedMask;

    if strcmp(camera0OrNot, 'notCamera0') 
        [temp_rows, temp_cols] = find(fNormalizedMask);
        temp_gridx(1) = min(temp_rows);
        temp_gridy(1) = min(temp_cols); 
        temp_gridx(2) = max(temp_rows);
        temp_gridy(2) = max(temp_cols); 
        DICpara.gridxyROIRange.gridx = temp_gridx;
        DICpara.gridxyROIRange.gridy = temp_gridy;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if DICpara.showImgOrNot
        figure;
        subplot(2,2,1); imshow(fNormalized'); title('fNormalized'); colorbar;
        subplot(2,2,2); imshow(gNormalized'); title('gNormalized'); colorbar;
        subplot(2,2,3); imshow(fNormalizedMask'); title('f mask'); colorbar;
        subplot(2,2,4); imshow(gNormalizedMask'); title('g mask'); colorbar;
    end
    %% Section 3: Compute an initial guess of the unknown displacement field
    fprintf('\n'); fprintf('------------ Section 3 Start ------------ \n')

    % Default: FFT initial search for each image
    DICpara.NewFFTSearch = 1; % DICpara.InitFFTSearchMethod = [];


    if ImgSeqNum <=ImgStartDataDrivenMode || DICpara.NewFFTSearch == 1
        % ====== Integer Search ======
        %  DICpara.InitFFTSearchMethod = 1;
        [DICpara,x0temp,y0temp,u,v,cc]= IntegerSearchQuadtree(fNormalized,gNormalized,file_name,DICpara,ImgSeqNum);

        % % %%%%% Optional codes to measure more gridded measurement points %%%%%
        % [DICpara,x0temp_f,y0temp_f,u_f,v_f,cc]= IntegerSearch(fNormalized,gNormalized,file_name,DICpara,ImgSeqNum);
        % 
        % xnodes = max([1+0.5*DICpara.winsize,DICpara.gridxyROIRange.gridx(1)])  ...
        %     : DICpara.winstepsize : min([size(fNormalized,1)-0.5*DICpara.winsize-1,DICpara.gridxyROIRange.gridx(2)]);
        % ynodes = max([1+0.5*DICpara.winsize,DICpara.gridxyROIRange.gridy(1)])  ...
        %     : DICpara.winstepsize : min([size(fNormalized,2)-0.5*DICpara.winsize-1,DICpara.gridxyROIRange.gridy(2)]);
        % 
        % [x0temp,y0temp] = ndgrid(xnodes,ynodes);   u_f_NotNanInd = find(~isnan(u_f(:)));
        % 
        % op1 = rbfcreate( [x0temp_f(u_f_NotNanInd),y0temp_f(u_f_NotNanInd)]',[u_f(u_f_NotNanInd)]','RBFFunction', 'thinplate');
        % rbfcheck_maxdiff = rbfcheck(op1); % Check: rbf thin-plate interpolation
        % if rbfcheck_maxdiff > 1e-3, disp('Please check rbf interpolation! Pause here.'); pause; end
        % u = rbfinterp([x0temp(:),y0temp(:)]', op1 );
        % 
        % op2 = rbfcreate( [x0temp_f(u_f_NotNanInd),y0temp_f(u_f_NotNanInd)]',[v_f(u_f_NotNanInd)]','RBFFunction', 'thinplate');
        % rbfcheck_maxdiff = rbfcheck(op2); % Check: rbf thin-plate interpolation
        % if rbfcheck_maxdiff > 1e-3, disp('Please check rbf interpolation! Pause here.'); pause; end
        % v = rbfinterp([x0temp(:),y0temp(:)]', op2 );
        % 
        % u = regularizeNd([x0temp(:),y0temp(:)],u(:),{xnodes',ynodes'},1e-3);
        % v = regularizeNd([x0temp(:),y0temp(:)],v(:),{xnodes',ynodes'},1e-3);

        % ====== FEM mesh set up ======
        [DICmesh_nonQuadtree] = MeshSetUp(x0temp,y0temp,DICpara); clear x0temp y0temp;

        % ====== Initial Value ======
        U0 = Init(u,v,cc.max,DICmesh_nonQuadtree.x0,DICmesh_nonQuadtree.y0,0); % PlotuvInit; [x0temp,y0temp,u,v,cc]= IntegerSearchMg(fNormalized,gNormalized,file_name,DICpara);

        % Zach Modified
        % Set zero at holes
        linearIndices1 = sub2ind(size(fNormalizedMask), DICmesh_nonQuadtree.coordinatesFEM(:,1), DICmesh_nonQuadtree.coordinatesFEM(:,2));
        MaskOrNot1 = fNormalizedMask(linearIndices1);

        nanIndex = find(MaskOrNot1<1);
        U0(2*nanIndex) = nan;
        U0(2*nanIndex-1) = nan;


        % ====== Deal with incremental mode ======
        if ImgSeqNum == 2
            RD.ResultFEMesh{ImgSeqNum-1} = struct( 'coordinatesFEM',DICmesh_nonQuadtree.coordinatesFEM,'elementsFEM',DICmesh_nonQuadtree.elementsFEM, ...
                'winsize',DICpara.winsize,'winstepsize',DICpara.winstepsize,'gridxyROIRange',DICpara.gridxyROIRange );
        else
            if incOrNot == 0
            else
                RD.ResultFEMesh{ImgSeqNum-1} = struct( 'coordinatesFEM',DICmesh_nonQuadtree.coordinatesFEM,'elementsFEM',DICmesh_nonQuadtree.elementsFEM, ...
                    'winsize',DICpara.winsize,'winstepsize',DICpara.winstepsize,'gridxyROIRange',DICpara.gridxyROIRange );
            end
        end

        DICmesh_nonQuadtree.elementMinSize = DICpara.winsizeMin; % min element size in the refined quadtree mesh
        % Notes:
        % Hanging nodes and sub-elements are placed on the last
        % All the void regions are generating nodes but we can ignore them
        % using maskfile later.
        [DICmesh_quadtree,DICpara,U0] = GenerateQuadtreeMesh(U0,Df,fNormalizedMask,DICmesh_nonQuadtree,DICpara); % Generate the quadtree mesh

        % ====== Compute f(X)-g(x+u) ======
        % PlotImgDiff(x0,y0,u,v,fNormalized,gNormalized);
        if ImgSeqNum == 2
            RD.ResultFEMeshEachFrame{ImgSeqNum-1} = struct( 'coordinatesFEM',DICmesh_quadtree.coordinatesFEM,'elementsFEM',DICmesh_quadtree.elementsFEM);
        else
            if incOrNot == 0
            else
                RD.ResultFEMeshEachFrame{ImgSeqNum-1} = struct( 'coordinatesFEM',DICmesh_quadtree.coordinatesFEM,'elementsFEM',DICmesh_quadtree.elementsFEM);
            end
        end

    else % imgSeqNum > ImgStartDataDrivenMode and NewFFT == 0 (For acc. mode.)
        % switch camera0OrNot
        % case 'camera0'
        %     inv_U0 = RD.ResultDisp{ImgSeqNum-2}.U';
        %     U0 = inv_U0(:);
        %     U0(2*nanIndex) = 0;
        %     U0(2*nanIndex-1) = 0;
        % case 'notCamera0'
        %     inv_U0 = RD.ResultDisp_corr{ImgSeqNum-2}.U';
        %     U0 = inv_U0(:);
        %     U0(2*nanIndex) = 0;
        %     U0(2*nanIndex-1) = 0;
    end

fprintf('------------ Section 3 Done ------------ \n \n')
%% Section 4: ALDIC Subproblem 1 -or- Local ICGN Subset DIC
fprintf('------------ Section 4 Start ------------ \n')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section is to solve the first local step in ALDIC: Subproblem 1
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ====== ALStep 1 Subproblem1: Local Subset DIC ======
mu=0; beta=0; tol=1e-2; ALSolveStep=1; ALSub1Time=zeros(6,1); ALSub2Time=zeros(6,1);
ConvItPerEle=zeros(size(DICmesh_quadtree.coordinatesFEM,1),6); ALSub1BadPtNum=zeros(6,1);
disp(['***** Start step',num2str(ALSolveStep),' Subproblem1 *****'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------ Start Local DIC IC-GN iteration ------
tic;
[USubpb1,FSubpb1,HtempPar,ALSub1Timetemp,ConvItPerEletemp,LocalICGNBadPtNumtemp,markCoordHoleStrain] = ...
    LocalICGNQuadtree(U0,DICmesh_quadtree.coordinatesFEM,Df,fNormalized,gNormalized,DICpara,'GaussNewton',tol,shapeFuncOrder);
ALSub1Time(ALSolveStep) = ALSub1Timetemp; ConvItPerEle(:,ALSolveStep) = ConvItPerEletemp; ALSub1BadPtNum(ALSolveStep) = LocalICGNBadPtNumtemp;
fprintf('Local: %.4f s\n', toc);

% % ====== Start: Zach-Optional: Thin-plate interpolate bad points =====
% coordinatesFEM = DICmesh_quadtree.coordinatesFEM;
% U = USubpb1; F = FSubpb1;
% nanindexU = find(isnan(U(1:2:end))==1); notnanindex = setdiff([1:1:size(coordinatesFEM,1)],nanindexU);
% %nanindexF = find(isnan(F(1:4:end))==1); notnanindexF = setdiff([1:1:size(coordinatesFEM,1)],nanindexF);
% 
% %%%%% Ux %%%%%%
% % Zach's debug
% % fi1 = rbfsplit([coordinatesFEM(notnanindex,1:2)],[U(2*notnanindex-1)],[coordinatesFEM(:,1:2)],2000,20);
% 
% op1 = rbfcreate( [coordinatesFEM(notnanindex,1:2)]',[U(2*notnanindex-1)]','RBFFunction', 'thinplate');
% % rbfcheck_maxdiff = rbfcheck(op1); % check rbf interpolation
% % if rbfcheck_maxdiff > 1e-3, disp('Please check rbf interpolation! Pause here.'); pause; end
% fi1 = rbfinterp([coordinatesFEM(:,1:2)]', op1 );
% 
% %%%%% Uy %%%%%%
% % Zach's debug
% % fi2 = rbfsplit([coordinatesFEM(notnanindex,1:2)],[U(2*notnanindex)],[coordinatesFEM(:,1:2)],2000,20);
% 
% op2 = rbfcreate( [coordinatesFEM(notnanindex,1:2)]',[U(2*notnanindex)]','RBFFunction', 'thinplate');
% % rbfcheck_maxdiff = rbfcheck(op2); % check rbf interpolation
% % if rbfcheck_maxdiff > 1e-3, disp('Please check rbf interpolation! Pause here.'); pause; end
% fi2 = rbfinterp([coordinatesFEM(:,1:2)]', op2 );
% 
% %%%%% Assemble [Ux, Uy] %%%%%
% U_rbf_thinplate = [fi1(:),fi2(:)]';  U_rbf_thinplate = U_rbf_thinplate(:);
% 
% %%%%% F11 %%%%%
% op =rbfcreate([coordinatesFEM(notnanindex,1:2)]', ...
%     F(4*notnanindex-3)','RBFFunction', 'thinplate');
% rbfcheck_maxdiff = rbfcheck(op);
% if rbfcheck_maxdiff > 1e-3, disp('Please check rbf interpolation! Pause here.'); pause; end
% fi11 = rbfinterp([coordinatesFEM(:,1:2)]', op );
% 
% %%%%% F21 %%%%%
% op =rbfcreate([coordinatesFEM(notnanindex,1:2)]', ...
%     F(4*notnanindex-2)','RBFFunction', 'thinplate');
% rbfcheck_maxdiff = rbfcheck(op);
% if rbfcheck_maxdiff > 1e-3, disp('Please check rbf interpolation! Pause here.'); pause; end
% fi21 = rbfinterp([coordinatesFEM(:,1:2)]', op );
% 
% %%%%% F12 %%%%%
% op =rbfcreate([coordinatesFEM(notnanindex,1:2)]', ...
%     F(4*notnanindex-1)','RBFFunction', 'thinplate');
% rbfcheck_maxdiff = rbfcheck(op);
% if rbfcheck_maxdiff > 1e-3, disp('Please check rbf interpolation! Pause here.'); pause; end
% fi12 = rbfinterp([coordinatesFEM(:,1:2)]', op );
% 
% %%%%% F22 %%%%%
% op =rbfcreate([coordinatesFEM(notnanindex,1:2)]', ...
%     F(4*notnanindex-0)','RBFFunction', 'thinplate');
% rbfcheck_maxdiff = rbfcheck(op);
% if rbfcheck_maxdiff > 1e-3, disp('Please check rbf interpolation! Pause here.'); pause; end
% fi22 = rbfinterp([coordinatesFEM(:,1:2)]', op );
% 
% %%%%% Assemble [F11,F21,F12,F22] %%%%%
% F_rbf_thinplate = [fi11(:),fi21(:),fi12(:),fi22(:)]';  F_rbf_thinplate = F_rbf_thinplate(:);
% USubpb1 = U_rbf_thinplate;   % FSubpb1 = F_rbf_thinplate;
% % ====== End: Zach-Optional: Thin-plate interpolate bad points =====



% ------ Plot ------
USubpb1World = USubpb1; USubpb1World(2:2:end) = -USubpb1(2:2:end); % FSubpb1World = FSubpb1;
% Plotdisp_show(USubpb1World,DICmesh_quadtree.coordinatesFEMWorld,DICmesh_quadtree.elementsFEM(:,1:4),DICpara,'EdgeColor');
% Plotstrain_show(FSubpb1World,DICmesh_quadtree.coordinatesFEMWorld,DICmesh_quadtree.elementsFEM(:,1:4),DICpara,'EdgeColor');
save(['Subpb1_step',num2str(ALSolveStep)],'USubpb1','FSubpb1');
fprintf('------------ Section 4 Done ------------ \n \n')


if UseGlobal
    fprintf('------------ Section 5 Start ------------ \n'); tic;
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This section is to solve the global step in ALDIC Subproblem 2
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % ======= ALStep 1 Subproblem 2: Global constraint =======
    % ------ Smooth displacements for a better F ------
    DICpara.DispFilterSize=0; DICpara.DispFilterStd=0; DICpara.StrainFilterSize=0; DICpara.StrainFilterStd=0; LevelNo=1;
    DICpara.DispSmoothness = 0; DICpara.StrainSmoothness = 0;
    if DICpara.DispSmoothness>1e-6, USubpb1 = funSmoothDispQuadtree(USubpb1,DICmesh_quadtree,DICpara); end
    % if DICpara.StrainSmoothness>1e-6, FSubpb1 = funSmoothStrainQuadtree(FSubpb1,DICmesh_quadtree,DICpara); end

    % ====== Define penalty parameter ======
    mu = 1e-3; udual = 0*FSubpb1; vdual = 0*USubpb1;
    betaList = [1e-3,1e-2,1e-1]*mean(DICpara.winstepsize).^2.*mu; % Tune beta in the betaList
    Err1 = zeros(length(betaList),1); Err2 = Err1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp(['***** Start step',num2str(ALSolveStep),' Subproblem2 *****']);
    DICpara.GaussPtOrder = 2; alpha = 0;  % No regularization added
    % ====== Solver using finite element method ======
    if ImgSeqNum == 2
        for tempk = 1:length(betaList)
            beta = betaList(tempk); display(['Try #',num2str(tempk),' beta = ',num2str(beta)]);
            GaussPtOrder=3; alpha=0; 
            [USubpb2] = Subpb2Quadtree(DICmesh_quadtree,DICpara.GaussPtOrder,beta,mu,USubpb1,FSubpb1,udual,vdual,alpha,mean(DICpara.winstepsize),0);
            % [FSubpb2,~,~] = funGlobalNodalStrainQuadtree(DICmesh_quadtree,USubpb2,DICpara.GaussPtOrder,0);

            Err1(tempk) = norm(USubpb1-USubpb2,2);
            % Err2(tempk) = norm(FSubpb1-FSubpb2,2);
        end

        Err1Norm = (Err1-mean(Err1))/std(Err1); % figure, plot(Err1Norm);
        %Err2Norm = (Err2-mean(Err2))/std(Err2); % figure, plot(Err2Norm);
        ErrSum = Err1Norm; % +Err2Norm; % figure, plot(ErrSum); title('Tune the best \beta in the subproblem 2');
        [~,indexOfbeta] = min(ErrSum);

        try % Tune the best beta by a quadratic polynomial 0fitting
            [fitobj] = fit(log10(betaList(indexOfbeta-1:1:indexOfbeta+1))',ErrSum(indexOfbeta-1:1:indexOfbeta+1),'poly2');
            p = coeffvalues(fitobj); beta = 10^(-p(2)/2/p(1));
        catch, beta = betaList(indexOfbeta);
        end
        display(['Best beta = ',num2str(beta)]);
    else
        try beta = DICpara.beta;
        catch, beta = 1e-3*mean(DICpara.winstepsize).^2.*mu;
        end
    end

    % Using the optimal beta to solve the ALDIC Subproblem 2 again
    if abs(beta-betaList(end))>abs(eps)
        [USubpb2] = Subpb2Quadtree(DICmesh_quadtree,DICpara.GaussPtOrder,beta,mu,USubpb1,FSubpb1,udual,vdual,alpha,mean(DICpara.winstepsize),0);
        %[FSubpb2,~,~] = funGlobalNodalStrainQuadtree(DICmesh_quadtree,USubpb2,DICpara.GaussPtOrder,0);
        ALSub2Time(ALSolveStep) = toc;
    end

    % ------- Smooth strain field --------
    if DICpara.DispSmoothness>1e-6, USubpb2 = funSmoothDispQuadtree(USubpb2,DICmesh_quadtree,DICpara); end
    % ------- Don't smooth strain fields near the boundary --------
    %for tempk=0:3, FSubpb2(4*DICmesh_quadtree.markCoordHoleEdge-tempk) = FSubpb1(4*DICmesh_quadtree.markCoordHoleEdge-tempk); end
    %if DICpara.StrainSmoothness>1e-6, FSubpb2 = funSmoothStrainQuadtree(0.1*FSubpb2+0.9*FSubpb1,DICmesh_quadtree,DICpara); end
    for tempk=0:1, USubpb2(2*markCoordHoleStrain-tempk) = USubpb1(2*markCoordHoleStrain-tempk); end
    %for tempk=0:3, FSubpb2(4*markCoordHoleStrain-tempk) = FSubpb1(4*markCoordHoleStrain-tempk); end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % ------- Save data ------
    FSubpb2 = FSubpb1;
    save(['Subpb2_step',num2str(ALSolveStep)],'USubpb2','FSubpb2');

    % ------ Plot ------
    USubpb2World = USubpb2; USubpb2World(2:2:end) = -USubpb2(2:2:end); % FSubpb2World = FSubpb2;
    close all; Plotdisp_show(USubpb2World,DICmesh_quadtree.coordinatesFEMWorld,DICmesh_quadtree.elementsFEM(:,1:4),DICpara,'EdgeColor');
    % Plotstrain_show(FSubpb2World,DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM(:,1:4),DICpara,'EdgeColor');

    % ======= Update dual variables =======
    udual = FSubpb2 - FSubpb1; vdual = USubpb2 - USubpb1;
    save(['uvdual_step',num2str(ALSolveStep)],'udual','vdual');
    fprintf('------------ Section 5 Done ------------ \n \n')


    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Section 6: ADMM iterations
    fprintf('------------ Section 6 Start ------------ \n')
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This section is the ADMM iteration, where both Subproblems 1 & 2 are solved iteratively.
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % ==================== ADMM AL Loop ==========================
    ALSolveStep = 1; tol2 = 1e-2; UpdateY = 1e4;
    HPar = cell(21,1); for tempj = 1:21, HPar{tempj} = HtempPar(:,tempj); end

    while (ALSolveStep < 3)
        ALSolveStep = ALSolveStep + 1;  % Update using the last step

        %%%%%%%% These lines can be used to further update each DIC subset window size %%%%%%%
        % Ftemp1 = FSubpb2(1:2:end); Ftemp2 = FSubpb2(2:2:end);
        % [DFtemp1,~,~] = funGlobalNodalStrainQuadtree(DICmesh,Ftemp1,DICpara.GaussPtOrder,0);
        % [DFtemp2,~,~] = funGlobalNodalStrainQuadtree(DICmesh,Ftemp2,DICpara.GaussPtOrder,0);
        %
        % winsize_x_ub1 = abs(2*FSubpb2(1:4:end)./DFtemp1(1:4:end)); winsize_x_ub2 = abs(2*FSubpb2(3:4:end)./DFtemp1(3:4:end));
        % winsize_y_ub1 = abs(2*FSubpb2(1:4:end)./DFtemp1(2:4:end)); winsize_y_ub2 = abs(2*FSubpb2(3:4:end)./DFtemp1(4:4:end));
        % winsize_x_ub3 = abs(2*FSubpb2(2:4:end)./DFtemp2(1:4:end)); winsize_x_ub4 = abs(2*FSubpb2(4:4:end)./DFtemp2(3:4:end));
        % winsize_y_ub3 = abs(2*FSubpb2(2:4:end)./DFtemp2(2:4:end)); winsize_y_ub4 = abs(2*FSubpb2(4:4:end)./DFtemp2(4:4:end));
        %
        % winsize_x_ub = round(min([winsize_x_ub1,winsize_x_ub2,winsize_x_ub3,winsize_x_ub4,DICpara.winsize*ones(length(winsize_x_ub1),1)],[],2));
        % winsize_x_List = max([winsize_x_ub, 10*ones(length(winsize_x_ub1),1)],[],2);
        % winsize_y_ub = round(min([winsize_y_ub1,winsize_y_ub2,winsize_y_ub3,winsize_y_ub4,DICpara.winsize*ones(length(winsize_y_ub1),1)],[],2));
        % winsize_y_List = max([winsize_y_ub, 10*ones(length(winsize_y_ub1),1)],[],2);
        % winsize_List = 2*ceil([winsize_x_List,winsize_y_List]/2);
        winsize_List = DICpara.winsize*ones(size(DICmesh_quadtree.coordinatesFEM,1),2);
        DICpara.winsize_List = winsize_List;


        %%%%%%%%%%%%%%%%%%%%%%% Subproblem 1 %%%%%%%%%%%%%%%%%%%%%%%%%
        disp(['***** Start step',num2str(ALSolveStep),' Subproblem1 *****']);
        tic; [USubpb1,~,ALSub1Timetemp,ConvItPerEletemp,LocalICGNBadPtNumtemp] = Subpb1Quadtree(...
            USubpb2,FSubpb2,udual,vdual,DICmesh_quadtree.coordinatesFEM,...
            Df,fNormalized,gNormalized,mu,beta,HPar,ALSolveStep,DICpara,'GaussNewton',tol);
        % FSubpb1 = FSubpb2; toc
        % for tempk=0:1, USubpb1(2*markCoordHoleStrain-tempk) = USubpb2(2*markCoordHoleStrain-tempk); end
        ALSub1Time(ALSolveStep) = ALSub1Timetemp; ConvItPerEle(:,ALSolveStep) = ConvItPerEletemp; ALSub1BadPtNum(ALSolveStep) = LocalICGNBadPtNumtemp;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ------  Manually find some bad points from Local Subset ICGN step ------
        % disp('--- Start to manually remove bad points --- \n')
        % disp('    Comment codes here if you do not have bad local points \n')
        % %%%%% Comment START %%%%%
        %  [USubpb1,FSubpb1] = funRemoveOutliersQuadtree(DICmesh,DICpara,USubpb1,FSubpb1);
        %  disp('--- Remove bad points done ---')
        % %%%%% Comment END %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        save(['Subpb1_step',num2str(ALSolveStep)],'USubpb1','FSubpb1');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ============== Subproblem 2 ==============
        disp(['***** Start step',num2str(ALSolveStep),' Subproblem2 *****'])
        tic; [USubpb2] = Subpb2Quadtree(DICmesh_quadtree,DICpara.GaussPtOrder,beta,mu,USubpb1,FSubpb1,udual,vdual,alpha,mean(DICpara.winstepsize),0);
        % [FSubpb2,~,~] = funGlobalNodalStrainQuadtree(DICmesh_quadtree,USubpb2,DICpara.GaussPtOrder,0);
        ALSub2Time(ALSolveStep) = toc; toc

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ------- Smooth strain field --------
        if DICpara.DispSmoothness>1e-6, USubpb2 = funSmoothDispQuadtree(USubpb2,DICmesh_quadtree,DICpara); end
        % ------- Don't change strain fields near the boundary --------
        %for tempk=0:3, FSubpb2(4*DICmesh_quadtree.markCoordHoleEdge-tempk) = FSubpb1(4*DICmesh_quadtree.markCoordHoleEdge-tempk); end
        %if DICpara.StrainSmoothness>1e-6, FSubpb2 = funSmoothStrainQuadtree(0.1*FSubpb2+0.9*FSubpb1,DICmesh_quadtree,DICpara); end
        for tempk=0:1, USubpb2(2*markCoordHoleStrain-tempk) = USubpb1(2*markCoordHoleStrain-tempk); end
        %for tempk=0:3, FSubpb2(4*markCoordHoleStrain-tempk) = FSubpb1(4*markCoordHoleStrain-tempk); end

        save(['Subpb2_step',num2str(ALSolveStep)],'USubpb2','FSubpb2');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Compute norm of UpdateY
        USubpb2_Old = load(['Subpb2_step',num2str(ALSolveStep-1)],'USubpb2');
        USubpb2_New = load(['Subpb2_step',num2str(ALSolveStep)],'USubpb2');
        USubpb1_Old = load(['Subpb1_step',num2str(ALSolveStep-1)],'USubpb1');
        USubpb1_New = load(['Subpb1_step',num2str(ALSolveStep)],'USubpb1');
        if (mod(ImgSeqNum-2,DICpara.ImgSeqIncUnit) ~= 0 && (ImgSeqNum>2)) || (ImgSeqNum < DICpara.ImgSeqIncUnit)
            UpdateY = norm((USubpb2_Old.USubpb2 - USubpb2_New.USubpb2), 2)/sqrt(size(USubpb2_Old.USubpb2,1));
            try
                UpdateY2 = norm((USubpb1_Old.USubpb1 - USubpb1_New.USubpb1), 2)/sqrt(size(USubpb1_Old.USubpb1,1));
            catch
            end
        end
        try
            disp(['Update local step  = ',num2str(UpdateY2)]);
            disp(['Update global step = ',num2str(UpdateY)]);
        catch
        end
        fprintf('*********************************** \n \n');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Update dual variables------------------------------
        udual = FSubpb2 - FSubpb1; vdual = USubpb2 - USubpb1;

        save(['uvdual_step',num2str(ALSolveStep)],'udual','vdual');
        try
            if UpdateY < tol2 || UpdateY2 < tol2
                break
            end
        catch
        end

    end
    fprintf('------------ Section 6 Done ------------ \n \n')
end

switch incOrNot
    case 0 % acc. mode
        % Step1: Save data
        if UseGlobal % AL or not
            RD.ResultDisp_acc{ImgSeqNum-1,1}.U = full(USubpb2);
            RD.ResultDisp{ImgSeqNum-1,1}.ALSub1BadPtNum = ALSub1BadPtNum;
            % RD.ResultDefGrad{ImgSeqNum-1,1}.F = full(FSubpb2); % tempFoamAL;
            % RD.DICmesh{ImgSeqNum-1} = DICmesh;
        else % regardless subproblem2
            RD.ResultDisp_acc{ImgSeqNum-1,1}.U = full(USubpb1);
            RD.ResultDisp{ImgSeqNum-1,1}.ALSub1BadPtNum = ALSub1BadPtNum;
            % RD.ResultDefGrad{ImgSeqNum-1}.F = full(FSubpb1); % tempFoamAL;
            % RD.DICmesh{ImgSeqNum-1} = DICmesh;
        end

        % Step2: interpolation due to Reference updating
        switch camera0OrNot
            case 'camera0'
                RD.ResultDisp{ImgSeqNum-1,1}.U = reshape(RD.ResultDisp_acc{ImgSeqNum-1}.U,2,size(RD.ResultDisp_acc{ImgSeqNum-1}.U,1)*0.5)';
            case 'notCamera0'
                tempX = RD.ResultDisp_acc{ImgSeqNum-1, 1}.U(1:2:end);
                tempY = RD.ResultDisp_acc{ImgSeqNum-1, 1}.U(2:2:end);

                % tempResultDisp(:,1) = rbfsplit(RD.ResultFEMeshEachFrame{1, 1}.coordinatesFEM,tempX,RD.Coordinates_corr,2000,20);
                % tempResultDisp(:,2) = rbfsplit(RD.ResultFEMeshEachFrame{1, 1}.coordinatesFEM,tempY,RD.Coordinates_corr,2000,20);

                rbfInterpX = rbfcreate(RD.ResultFEMeshEachFrame{1, 1}.coordinatesFEM',tempX');
                rbfInterpY = rbfcreate(RD.ResultFEMeshEachFrame{1, 1}.coordinatesFEM',tempY');

                tempResultDisp(:,1) = rbfinterp(RD.Coordinates_corr',rbfInterpX);
                tempResultDisp(:,2) = rbfinterp(RD.Coordinates_corr',rbfInterpY);

                % Based on Left first frame
                RD.ResultDisp{ImgSeqNum-1,1}.U = tempResultDisp;

                % Based on Right first frame
                RD.ResultDisp_corr{ImgSeqNum-1,1}.U = [RD.ResultDisp_acc{ImgSeqNum-1,1}.U(1:2:end),RD.ResultDisp_acc{ImgSeqNum-1,1}.U(2:2:end)];

        end

    case 1 % inc. mode
        % Step1: Save data
        if UseGlobal
            RD.ResultDisp_inc{ImgSeqNum-1,1}.U = full(USubpb2);
            RD.ResultDisp{ImgSeqNum-1,1}.ALSub1BadPtNum = ALSub1BadPtNum;
            % RD.ResultDefGrad{ImgSeqNum-1,1}.F = full(FSubpb2); % tempFoamAL;
            % RD.DICmesh{ImgSeqNum-1} = DICmesh;
        else % regardless subproblem2
            RD.ResultDisp_inc{ImgSeqNum-1}.U = full(USubpb1);
            RD.ResultDisp{ImgSeqNum-1}.ALSub1BadPtNum = ALSub1BadPtNum;
            % RD.ResultDefGrad{ImgSeqNum-1}.F = full(FSubpb1); % tempFoamAL;
            % RD.DICmesh{ImgSeqNum-1} = DICmesh;
        end

        % Step2: interpolation due to Reference updating
        switch camera0OrNot
            case 'camera0'
                if ImgSeqNum == 2
                    RD.ResultDisp{1,1}.U = reshape(RD.ResultDisp_inc{ImgSeqNum-1}.U,2,size(RD.ResultDisp_inc{ImgSeqNum-1}.U,1)*0.5)';

                else
                    tempX = RD.ResultDisp_inc{ImgSeqNum-1, 1}.U(1:2:end);
                    tempY = RD.ResultDisp_inc{ImgSeqNum-1, 1}.U(2:2:end);
                    % tempResultDisp(:,1) = rbfsplit(RD.ResultFEMeshEachFrame{ImgSeqNum-1, 1}.coordinatesFEM,tempX,RD.ResultFEMeshEachFrame{1,1}.coordinatesFEM,2000,20);
                    % tempResultDisp(:,2) = rbfsplit(RD.ResultFEMeshEachFrame{ImgSeqNum-1, 1}.coordinatesFEM,tempY,RD.ResultFEMeshEachFrame{1,1}.coordinatesFEM,2000,20);
                    %
                    rbfInterpX = rbfcreate(RD.ResultFEMeshEachFrame{ImgSeqNum-1, 1}.coordinatesFEM',tempX');
                    rbfInterpY = rbfcreate(RD.ResultFEMeshEachFrame{ImgSeqNum-1, 1}.coordinatesFEM',tempY');

                    tempResultDisp(1,:) = rbfinterp(RD.ResultFEMeshEachFrame{1,1}.coordinatesFEM',rbfInterpX);
                    tempResultDisp(2,:) = rbfinterp(RD.ResultFEMeshEachFrame{1,1}.coordinatesFEM',rbfInterpY);
                    RD.ResultDisp{ImgSeqNum-1,1}.U = RD.ResultDisp{ImgSeqNum-2,1}.U  + tempResultDisp'; % disp_inc 是相较于最初第一张ROI的位移
                end
            case 'notCamera0'
                tempX = RD.ResultDisp_inc{ImgSeqNum-1, 1}.U(1:2:end);
                tempY = RD.ResultDisp_inc{ImgSeqNum-1, 1}.U(2:2:end);
                % tic;
                % tempResultDisp(:,1) = rbfsplit(RD.ResultFEMeshEachFrame{ImgSeqNum-1, 1}.coordinatesFEM,tempX,RD.Coordinates_corr,min(2000,size(tempY,1)),20);
                % tempResultDisp(:,2) = rbfsplit(RD.ResultFEMeshEachFrame{ImgSeqNum-1, 1}.coordinatesFEM,tempY,RD.Coordinates_corr,min(2000,size(tempY,1)),20);
                % % toc;
                rbfInterpX = rbfcreate(RD.ResultFEMeshEachFrame{ImgSeqNum-1, 1}.coordinatesFEM',tempX');
                rbfInterpY = rbfcreate(RD.ResultFEMeshEachFrame{ImgSeqNum-1, 1}.coordinatesFEM',tempY');

                tempResultDisp(:,1) = rbfinterp(RD.Coordinates_corr',rbfInterpX);
                tempResultDisp(:,2) = rbfinterp(RD.Coordinates_corr',rbfInterpY);

                if ImgSeqNum == 2
                    RD.ResultDisp{ImgSeqNum-1,1}.U = tempResultDisp; % disp_inc 是相较于最初第一张ROI的位移
                else
                    RD.ResultDisp{ImgSeqNum-1,1}.U = RD.ResultDisp{ImgSeqNum-2,1}.U  + tempResultDisp; % disp_inc 是相较于最初第一张ROI的位移
                end

        end
end

if DICpara.showImgOrNot
    % ------ Plot ------
    % if DICpara.DICIncOrNot == 1
    % U_inc = reshape(RD.ResultDisp_inc{ImgSeqNum,1}.',[],1);
    % close all; Plotuv(U_inc,DICmesh_init.x0,DICmesh_init.y0World);
    % Plotdisp_show(U_inc,DICmesh_init.coordinatesFEMWorld,DICmesh_init.elementsFEM);
    % else

    if UseGlobal
        USubpb2World = USubpb2; USubpb2World(2:2:end) = -USubpb2(2:2:end);
        FSubpb2World = FSubpb2; FSubpb2World(2:4:end) = -FSubpb2World(2:4:end); FSubpb2World(3:4:end) = -FSubpb2World(3:4:end);
        close all; % plotuv(USubpb2World,DICmesh_quadtree.x0,DICmesh_quadtree.y0World);
        Plotdisp_show(USubpb2World,DICmesh_quadtree.coordinatesFEMWorld,DICmesh_quadtree.elementsFEM);
        Plotdisp_show(USubpb1World,DICmesh_quadtree.coordinatesFEMWorld,DICmesh_quadtree.elementsFEM(:,1:4),DICpara,'EdgeColor');
    else
        USubpb2World = USubpb1; USubpb2World(2:2:end) = -USubpb1(2:2:end);
        FSubpb2World = FSubpb1; FSubpb2World(2:4:end) = -FSubpb2World(2:4:end); FSubpb2World(3:4:end) = -FSubpb2World(3:4:end);
        close all; % plotuv(USubpb2World,DICmesh_quadtree.x0,DICmesh_quadtree.y0World);
        Plotdisp_show(USubpb2World,DICmesh_quadtree.coordinatesFEMWorld,DICmesh_quadtree.elementsFEM);
        Plotdisp_show(USubpb1World,DICmesh_quadtree.coordinatesFEMWorld,DICmesh_quadtree.elementsFEM(:,1:4),DICpara,'EdgeColor');
    end
    
    % end
else
end


% ------ Save results ------
% Find img name and save all the results

% [~,imgname,imgext] = fileparts(file_name{1,end});
% results_name = ['results_',imgname,'_ws',num2str(DICpara.winsize),'_st',num2str(DICpara.winstepsize),'.mat'];
% save(results_name, 'file_name','DICpara','DICmesh','RD','ALSub1Time','ALSub2Time','ALSolveStep');
end
end









