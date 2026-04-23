function [RD] = TemporalMatch_quadtree_ST1(DICpara,file_name,ImgMask,ImgNormalized,RD,StereoInfo,camera0OrNot,shapeFuncOrder)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This part is designed for temporal DIC, which means dealing with left and
% right images separately using 2D-ALDIC.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% debug options
UseGlobal = 1; % 0=Local DIC only (fast), 1=ALDIC with global constraints
DICpara.showImgOrNot = 0;  % Set to 1 for debugging, 0 for batch processing

% Auto-enable parallel computing if Parallel Computing Toolbox available
if DICpara.ClusterNo <= 1
    try
        nCores = feature('numcores');
        if nCores > 1 && license('test', 'Distrib_Computing_Toolbox')
            DICpara.ClusterNo = nCores;
            fprintf('Auto-enabled parallel ICGN with %d workers\n', nCores);
        end
    catch
        % Keep sequential if detection fails
    end
end

% Init
imageNum = length(ImgNormalized);
incOrNot = DICpara.DICIncOrNot; % 0: acc, 1: inc



% Data-driven mode (use previous-frame result as initial guess after frame
% 7) removed — the alternate branch was fully commented out. The "smart
% FFT skip" logic in the frame loop below already covers the inc-mode
% reuse case without needing a mode flag.
ImgStartDataDrivenMode = length(ImgNormalized) + 1;

%% For loop
for ImgSeqNum = 2 : imageNum
    fprintf('--- Frame %d/%d ---\n', ImgSeqNum, length(ImgNormalized));

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
    % Section 3: Initial guess

    % Decide FFT search strategy for initial guess
    useFFTSearch = true;
    if incOrNot == 1 && ImgSeqNum >= 3
        try
            prevDisp = RD.ResultDisp_inc{ImgSeqNum-2, 1}.U;  % [2N_prev x 1]
            prevCoords = RD.ResultFEMeshEachFrame{ImgSeqNum-2, 1}.coordinatesFEM;  % [N_prev x 2]
            prevU_x = prevDisp(1:2:end); prevU_y = prevDisp(2:2:end);
            validFrac = 1 - nnz(isnan(prevU_x)) / numel(prevU_x);
            if validFrac > 0.5, useFFTSearch = false; end
        catch
            % Fall back to FFT search if previous result unavailable
        end
    end

    DICpara.NewFFTSearch = 1;

    if ImgSeqNum <=ImgStartDataDrivenMode || DICpara.NewFFTSearch == 1
        if useFFTSearch
            % ====== Integer Search (FFT) ======
            [DICpara,x0temp,y0temp,u,v,cc]= IntegerSearchQuadtree(fNormalized,gNormalized,file_name,DICpara,ImgSeqNum);
        else
            % ====== Skip FFT: use previous frame displacement as initial guess ======
            fprintf('  Skipping FFT: reusing previous frame result (%.0f%% valid)\n', 100*validFrac);
            winsize_v = DICpara.winsize; if isscalar(winsize_v), winsize_v = [winsize_v winsize_v]; end
            winstep_v = DICpara.winstepsize; if isscalar(winstep_v), winstep_v = [winstep_v winstep_v]; end
            gridxROI = DICpara.gridxyROIRange.gridx;
            gridyROI = DICpara.gridxyROIRange.gridy;
            XList = gridxROI(1) : winstep_v(1) : gridxROI(2) - winsize_v(1) - winstep_v(1);
            YList = gridyROI(1) : winstep_v(2) : gridyROI(2) - winsize_v(2) - winstep_v(2);
            [x0temp, y0temp] = ndgrid(XList, YList);

            % Interpolate previous frame's displacement to current grid
            prevNonNan = ~isnan(prevU_x);
            Fi_u = scatteredInterpolant(prevCoords(prevNonNan,1), prevCoords(prevNonNan,2), ...
                prevU_x(prevNonNan), 'natural', 'none');
            Fi_v = scatteredInterpolant(prevCoords(prevNonNan,1), prevCoords(prevNonNan,2), ...
                prevU_y(prevNonNan), 'natural', 'none');
            u = reshape(Fi_u(x0temp(:), y0temp(:)), size(x0temp));
            v = reshape(Fi_v(x0temp(:), y0temp(:)), size(x0temp));
            u(isnan(u)) = 0; v(isnan(v)) = 0;  % Safe fallback for points outside convex hull
            cc.max = ones(size(u));  % Dummy CC (not from FFT search)
        end


        % ====== FEM mesh set up ======
        [DICmesh_nonQuadtree] = MeshSetUp(x0temp,y0temp,DICpara); clear x0temp y0temp;

        % ====== Initial Value ======
        U0 = Init(u,v,cc.max,DICmesh_nonQuadtree.x0,DICmesh_nonQuadtree.y0,0); % PlotuvInit; [x0temp,y0temp,u,v,cc]= IntegerSearchMg(fNormalized,gNormalized,file_name,DICpara);

        % Zach Modified
        % Set zero at holes
        linearIndices1 = sub2ind(size(fNormalizedMask), round(DICmesh_nonQuadtree.coordinatesFEM(:,1)), round(DICmesh_nonQuadtree.coordinatesFEM(:,2)));
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

    end

% Section 4: Local ICGN
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section is to solve the first local step in ALDIC: Subproblem 1
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ====== ALStep 1 Subproblem1: Local Subset DIC ======
mu=0; beta=0; tol=1e-2; ALSolveStep=1; ALSub1Time=zeros(6,1); ALSub2Time=zeros(6,1);
ConvItPerEle=zeros(size(DICmesh_quadtree.coordinatesFEM,1),6); ALSub1BadPtNum=zeros(6,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------ Start Local DIC IC-GN iteration ------
tic;
[USubpb1,FSubpb1,HtempPar,ALSub1Timetemp,ConvItPerEletemp,LocalICGNBadPtNumtemp,markCoordHoleStrain] = ...
    LocalICGNQuadtree(U0,DICmesh_quadtree.coordinatesFEM,Df,fNormalized,gNormalized,DICpara,'GaussNewton',tol,shapeFuncOrder);
ALSub1Time(ALSolveStep) = ALSub1Timetemp; ConvItPerEle(:,ALSolveStep) = ConvItPerEletemp; ALSub1BadPtNum(ALSolveStep) = LocalICGNBadPtNumtemp;

%% Remove bad points
coordinatesFEM = DICmesh_quadtree.coordinatesFEM;
U = USubpb1; 
F = FSubpb1;

nanindexU = find(isnan(U(1:2:end))==1);
notnanindex  = setdiff(1:size(coordinatesFEM,1), nanindexU);

nanindexF = find(isnan(F(1:4:end))==1);
notnanindexF = setdiff(1:size(coordinatesFEM,1), nanindexF);

% ===== U interpolate (only if needed) =====
if ~isempty(nanindexU) && numel(notnanindex) > 10
    src_xy = coordinatesFEM(notnanindex, 1:2);
    query_xy = coordinatesFEM(:, 1:2);
    Fi_u = scatteredInterpolant(src_xy(:,1), src_xy(:,2), U(2*notnanindex-1), 'natural', 'none');
    Fi_v = scatteredInterpolant(src_xy(:,1), src_xy(:,2), U(2*notnanindex),   'natural', 'none');
    fi1 = Fi_u(query_xy(:,1), query_xy(:,2));
    fi2 = Fi_v(query_xy(:,1), query_xy(:,2));
    U_interp = [fi1(:), fi2(:)]';
    USubpb1 = U_interp(:);
end

% ===== F interpolate (only if needed) =====
if ~isempty(nanindexF) && numel(notnanindexF) > 10
    src_xy = coordinatesFEM(notnanindexF, 1:2);
    query_xy = coordinatesFEM(:, 1:2);
    Fi_11 = scatteredInterpolant(src_xy(:,1), src_xy(:,2), F(4*notnanindexF-3), 'natural', 'none');
    Fi_21 = scatteredInterpolant(src_xy(:,1), src_xy(:,2), F(4*notnanindexF-2), 'natural', 'none');
    Fi_12 = scatteredInterpolant(src_xy(:,1), src_xy(:,2), F(4*notnanindexF-1), 'natural', 'none');
    Fi_22 = scatteredInterpolant(src_xy(:,1), src_xy(:,2), F(4*notnanindexF-0), 'natural', 'none');
    fi11 = Fi_11(query_xy(:,1), query_xy(:,2));
    fi21 = Fi_21(query_xy(:,1), query_xy(:,2));
    fi12 = Fi_12(query_xy(:,1), query_xy(:,2));
    fi22 = Fi_22(query_xy(:,1), query_xy(:,2));
    F_interp = [fi11(:), fi21(:), fi12(:), fi22(:)]';
    FSubpb1 = F_interp(:);
end

% sanity (silent unless NaN remain)
nNanU = nnz(isnan(USubpb1)); nNanF = nnz(isnan(FSubpb1));
if nNanU > 0 || nNanF > 0, fprintf('  Post-interp NaN: U=%d, F=%d\n', nNanU, nNanF); end



% Section 4 done


if UseGlobal
    tic; % Section 5: Global constraint
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This section is to solve the global step in ALDIC Subproblem 2
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % ======= ALStep 1 Subproblem 2: Global constraint =======
    % ------ Smooth displacements for a better F ------
    % 2D smoothing inside ADMM, controlled by DICpara.Smooth2DTimes
    % (0 = disabled, default). Filter size/std also come from DICpara.
    if ~isfield(DICpara,'Smooth2DTimes') || isempty(DICpara.Smooth2DTimes)
        DICpara.Smooth2DTimes = 0;
    end
    LevelNo = 1;
    for smoothN = 1:DICpara.Smooth2DTimes
        USubpb1 = funSmoothDispQuadtree(USubpb1,DICmesh_quadtree,DICpara);
    end

    % ====== Define penalty parameter ======
    mu = 1e-3; udual = zeros(size(FSubpb1)); vdual = zeros(size(USubpb1));
    betaList = [1e-3,1e-2,1e-1]*mean(DICpara.winstepsize).^2.*mu; % Tune beta in the betaList
    Err1 = zeros(length(betaList),1); Err2 = Err1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Subproblem2 step
    DICpara.GaussPtOrder = 2; alpha = 0;  % No regularization added
    % ====== Solver using finite element method ======
    if ImgSeqNum == 2
        for tempk = 1:length(betaList)
            beta = betaList(tempk);
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
        fprintf('  beta=%.2e\n', beta);
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

    % ======= Update dual variables =======
    udual = FSubpb2 - FSubpb1; vdual = USubpb2 - USubpb1;


    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Section 6: ADMM iterations
    % Section 6: ADMM iterations
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This section is the ADMM iteration, where both Subproblems 1 & 2 are solved iteratively.
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % ==================== ADMM AL Loop ==========================
    ALSolveStep = 1; tol2 = 1e-2; UpdateY = 1e4;
    HPar = cell(21,1); for tempj = 1:21, HPar{tempj} = HtempPar(:,tempj); end
    USubpb1_prev = USubpb1; USubpb2_prev = USubpb2;  % In-memory convergence tracking

    while (ALSolveStep < 3)
        ALSolveStep = ALSolveStep + 1;  % Update using the last step

        % Uniform winsize across mesh (adaptive per-node winsize not used).
        winsize_List = DICpara.winsize*ones(size(DICmesh_quadtree.coordinatesFEM,1),2);
        DICpara.winsize_List = winsize_List;

        % ====== Subproblem 1 (local ICGN) ======
        tic; [USubpb1,~,ALSub1Timetemp,ConvItPerEletemp,LocalICGNBadPtNumtemp] = Subpb1Quadtree(...
            USubpb2,FSubpb2,udual,vdual,DICmesh_quadtree.coordinatesFEM,...
            Df,fNormalized,gNormalized,mu,beta,HPar,ALSolveStep,DICpara,'GaussNewton',tol);
        ALSub1Time(ALSolveStep) = ALSub1Timetemp; ConvItPerEle(:,ALSolveStep) = ConvItPerEletemp; ALSub1BadPtNum(ALSolveStep) = LocalICGNBadPtNumtemp;

        % ====== Subproblem 2 (global FEM) ======
        tic; [USubpb2] = Subpb2Quadtree(DICmesh_quadtree,DICpara.GaussPtOrder,beta,mu,USubpb1,FSubpb1,udual,vdual,alpha,mean(DICpara.winstepsize),0);
        ALSub2Time(ALSolveStep) = toc;

        % ------- Optional smoothing --------
        for smoothN = 1:DICpara.Smooth2DTimes
            USubpb2 = funSmoothDispQuadtree(USubpb2,DICmesh_quadtree,DICpara);
        end
        % Preserve hole-edge displacements from Subpb1 (no global blur)
        for tempk=0:1, USubpb2(2*markCoordHoleStrain-tempk) = USubpb1(2*markCoordHoleStrain-tempk); end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Compute norm of UpdateY (in-memory, no file I/O)
        if (mod(ImgSeqNum-2,DICpara.ImgSeqIncUnit) ~= 0 && (ImgSeqNum>2)) || (ImgSeqNum < DICpara.ImgSeqIncUnit)
            UpdateY = norm((USubpb2_prev - USubpb2), 2)/sqrt(size(USubpb2_prev,1));
            try
                UpdateY2 = norm((USubpb1_prev - USubpb1), 2)/sqrt(size(USubpb1_prev,1));
            catch ME
                warning('%s: %s', ME.identifier, ME.message);
            end
        end
        USubpb1_prev = USubpb1; USubpb2_prev = USubpb2;
        try
            fprintf('  ADMM step%d: local=%.2e, global=%.2e\n', ALSolveStep, UpdateY2, UpdateY);
        catch
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Update dual variables------------------------------
        udual = FSubpb2 - FSubpb1; vdual = USubpb2 - USubpb1;

        try
            if UpdateY < tol2 || UpdateY2 < tol2
                break
            end
        catch ME
            warning('%s: %s', ME.identifier, ME.message);
        end

    end
    % Section 6 done
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

                src_coords = RD.ResultFEMeshEachFrame{1, 1}.coordinatesFEM;
                query_coords = RD.Coordinates_corr;
                Fi_X = scatteredInterpolant(src_coords(:,1), src_coords(:,2), tempX, 'natural', 'none');
                Fi_Y = scatteredInterpolant(src_coords(:,1), src_coords(:,2), tempY, 'natural', 'none');
                tempResultDisp(:,1) = Fi_X(query_coords(:,1), query_coords(:,2));
                tempResultDisp(:,2) = Fi_Y(query_coords(:,1), query_coords(:,2));

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
                    src_coords = RD.ResultFEMeshEachFrame{ImgSeqNum-1, 1}.coordinatesFEM;
                    Fi_X = scatteredInterpolant(src_coords(:,1), src_coords(:,2), tempX, 'natural', 'none');
                    Fi_Y = scatteredInterpolant(src_coords(:,1), src_coords(:,2), tempY, 'natural', 'none');
                    % FIX: Evaluate at DEFORMED positions (Lagrangian tracking)
                    deformed_pos_L = RD.ResultFEMeshEachFrame{1,1}.coordinatesFEM + RD.ResultDisp{ImgSeqNum-2,1}.U;
                    tempResultDisp(1,:) = Fi_X(deformed_pos_L(:,1), deformed_pos_L(:,2))';
                    tempResultDisp(2,:) = Fi_Y(deformed_pos_L(:,1), deformed_pos_L(:,2))';
                    RD.ResultDisp{ImgSeqNum-1,1}.U = RD.ResultDisp{ImgSeqNum-2,1}.U  + tempResultDisp'; % disp_inc accumulated relative to the first-frame ROI
                end
            case 'notCamera0'
                tempX = RD.ResultDisp_inc{ImgSeqNum-1, 1}.U(1:2:end);
                tempY = RD.ResultDisp_inc{ImgSeqNum-1, 1}.U(2:2:end);
                % tic;
                % tempResultDisp(:,1) = rbfsplit(RD.ResultFEMeshEachFrame{ImgSeqNum-1, 1}.coordinatesFEM,tempX,RD.Coordinates_corr,min(2000,size(tempY,1)),20);
                % tempResultDisp(:,2) = rbfsplit(RD.ResultFEMeshEachFrame{ImgSeqNum-1, 1}.coordinatesFEM,tempY,RD.Coordinates_corr,min(2000,size(tempY,1)),20);
                % % toc;
                src_coords = RD.ResultFEMeshEachFrame{ImgSeqNum-1, 1}.coordinatesFEM;
                Fi_X = scatteredInterpolant(src_coords(:,1), src_coords(:,2), tempX, 'natural', 'none');
                Fi_Y = scatteredInterpolant(src_coords(:,1), src_coords(:,2), tempY, 'natural', 'none');

                if ImgSeqNum == 2
                    tempResultDisp(:,1) = Fi_X(RD.Coordinates_corr(:,1), RD.Coordinates_corr(:,2));
                    tempResultDisp(:,2) = Fi_Y(RD.Coordinates_corr(:,1), RD.Coordinates_corr(:,2));
                    RD.ResultDisp{ImgSeqNum-1,1}.U = tempResultDisp; % disp_inc relative to first frame
                else
                    % FIX: Evaluate at DEFORMED positions (Lagrangian tracking)
                    deformed_pos_R = RD.Coordinates_corr + RD.ResultDisp{ImgSeqNum-2,1}.U;
                    tempResultDisp(:,1) = Fi_X(deformed_pos_R(:,1), deformed_pos_R(:,2));
                    tempResultDisp(:,2) = Fi_Y(deformed_pos_R(:,1), deformed_pos_R(:,2));
                    RD.ResultDisp{ImgSeqNum-1,1}.U = RD.ResultDisp{ImgSeqNum-2,1}.U + tempResultDisp; % accumulated disp
                end

        end
end

if DICpara.showImgOrNot
    if UseGlobal
        UDispWorld = USubpb2; UDispWorld(2:2:end) = -USubpb2(2:2:end);
    else
        UDispWorld = USubpb1; UDispWorld(2:2:end) = -USubpb1(2:2:end);
    end
    close all;
    Plotdisp_show(UDispWorld,DICmesh_quadtree.coordinatesFEMWorld,DICmesh_quadtree.elementsFEM(:,1:4),DICpara,'EdgeColor');
end

end  % for ImgSeqNum
end  % function









