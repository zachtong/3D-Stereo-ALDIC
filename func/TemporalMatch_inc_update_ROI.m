function [RD] = TemporalMatch_inc_update_ROI(DICpara,file_name,ImgNormalized,RD,StereoInfo,camera0OrNot)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This part is designed for temporal DIC, which means dealing with left and
% right images separately using 2D-ALDIC.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% debug options
UseGlobal = 1;
showImgOrNot = 1;

% Init
imageNum = length(ImgNormalized);

switch camera0OrNot
    case 'camera0'

    case 'notCamera0'
        RD.Coordinates_corr = StereoInfo.ResultFEMesh_corr; % Backup: The first picture from the right camera corresponds
        % to the picture from the left camera and needs to be interpolated back.
        % [DICpara] = extendROI(DICpara,RD);
end


%% For loop
for ImgSeqNum = 2 : imageNum
    disp(['Current image frame #: ', num2str(ImgSeqNum),'/',num2str(length(ImgNormalized))]);

    %% Task 1: manually renew the ROI
    switch camera0OrNot
        case 'camera0'
            if    ImgSeqNum > 2
                %% Task 1234: manually renew the ROI
                ImgTemp = imread(file_name{1,ImgSeqNum-1});
                imshow(ImgTemp);
                gridx = zeros(1,2); gridy = zeros(1,2);
                [gridx(1), gridy(1)] = ginput(1);
                [gridx(2), gridy(2)] = ginput(1);
                DICpara.gridxyROIRange.gridx = round(gridx);
                DICpara.gridxyROIRange.gridy = round(gridy);
            end

        case 'notCamera0'
            %% Task 1234: manually renew the ROI
            ImgTemp = imread(file_name{1,ImgSeqNum-1});
            imshow(ImgTemp);
            gridx = zeros(1,2); gridy = zeros(1,2);
            [gridx(1), gridy(1)] = ginput(1);
            [gridx(2), gridy(2)] = ginput(1);
            DICpara.gridxyROIRange.gridx = round(gridx);
            DICpara.gridxyROIRange.gridy = round(gridy);
    end


    fNormalized = ImgNormalized{ImgSeqNum-1}; % Load previous reference image frame
    gNormalized = ImgNormalized{ImgSeqNum}; % Load current deformed image frame
    Df = funImgGradient(fNormalized,fNormalized);

    %% Section 3: Compute an initial guess of the unknown displacement field
    fprintf('\n'); fprintf('------------ Section 3 Start ------------ \n')
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This section is to find or update an initial guess of the unknown displacements.
    % The key idea is to either to use a new FFT-based cross correlation peak fitting,
    % or use the results from previous frames to compute a new initial guess for the next frame;
    % Particularly Calin the incremental mode DIC, the reference image can also be updated, e.g.,
    % " fNormalized = ImgNormalized{ImgSeqNum-mod(ImgSeqNum-1,ImgSeqIncUnit)}; "
    %
    % DICpara.NewFFTSearch = 0; % If you want to apply the FFT-based cross correlation to
    % compute the initial guess for each frame, please make sure that "DICpara.NewFFTSearch = 0".
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%% One practical strategy is to let first 7 frames do the FFT-based
    %%%%% cross correlation and then using the data driven method to estimate
    %%%%% initial guesses for other frames
    if ImgSeqNum == 2
        DICpara.NewFFTSearch = 1; DICpara.InitFFTSearchMethod = [];
    elseif ImgSeqNum < 7
        DICpara.NewFFTSearch = 1; % Use FFT-based cross correlation to compute the initial guess
    else
        DICpara.NewFFTSearch = 0; % Apply data driven method to estimate initial guesses for later frames
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% task --------------------ZACH!!!
    DICpara.NewFFTSearch = 1; % Do not use the data driven method

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ImgSeqNum == 2 || DICpara.NewFFTSearch == 1 % Apply FFT-based cross correlation to compute the initial guess

        % ====== Integer Search ======
        [DICpara,x0temp,y0temp,u,v,cc]= IntegerSearch(fNormalized,gNormalized,file_name,DICpara,showImgOrNot,ImgSeqNum);

        % ====== FEM mesh set up ======
        [DICmesh] = MeshSetUp(x0temp,y0temp,DICpara); clear x0temp y0temp;
        % ====== Initial Value ======
        U0 = Init(u,v,cc.max,DICmesh.x0,DICmesh.y0,0); % PlotuvInit; [x0temp,y0temp,u,v,cc]= IntegerSearchMg(fNormalized,gNormalized,file_name,DICpara);
        % ====== Deal with incremental mode ======
        fNormalizedNewIndex = ImgSeqNum-mod(ImgSeqNum-2,DICpara.ImgSeqIncUnit)-1;
        if DICpara.ImgSeqIncUnit == 1, fNormalizedNewIndex = fNormalizedNewIndex-1; end
        RD.ResultFEMesh{1+floor(fNormalizedNewIndex/DICpara.ImgSeqIncUnit)} = ... % To save first mesh info
            struct( 'coordinatesFEM',DICmesh.coordinatesFEM,'elementsFEM',DICmesh.elementsFEM, ...
            'winsize',DICpara.winsize,'winstepsize',DICpara.winstepsize,'gridxyROIRange',DICpara.gridxyROIRange );

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif mod(ImgSeqNum-2,DICpara.ImgSeqIncUnit) == 0 % To update ref image in incremental mode
        fNormalizedNewIndex = ImgSeqNum-mod(ImgSeqNum-2,DICpara.ImgSeqIncUnit)-1;
        if DICpara.ImgSeqIncUnit == 1,  fNormalizedNewIndex = fNormalizedNewIndex-1; end
        fNormalized = ImgNormalized{fNormalizedNewIndex}; % Update reference
        [DICpara,DICmesh] = ReadImageRefUpdate(file_name,ImgSeqNum,RD.ResultDisp{ImgSeqNum-2}.U,DICpara,DICmesh); % Update reference image if needed;
        U0 = zeros(2*size(DICmesh.coordinatesFEM,1),1); % PlotuvInit;
        RD.ResultFEMesh{1+floor(fNormalizedNewIndex/DICpara.ImgSeqIncUnit)} = ... % To save first mesh info
            struct( 'coordinatesFEM',DICmesh.coordinatesFEM,'elementsFEM',DICmesh.elementsFEM, ...
            'winsize',DICpara.winsize,'winstepsize',DICpara.winstepsize,'gridxyROIRange',DICpara.gridxyROIRange );

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else % Use the solved results from the last frame as the new initial guess
        if ImgSeqNum < 7 % Import previous U for ImgSeqNum [2,6]
            U0 = RD.ResultDisp{ImgSeqNum-2}.U;

        else % When ImgSeqNum > 6: POD predicts next disp U0 from previous results of (ImgSeqNum+[-5:1:-1])
            nTime = 5; np = length(RD.ResultDisp{ImgSeqNum-2}.U)/2; % "nTime" value 5 is an empirical value, can be changed.
            T_data_u = zeros(nTime,np); T_data_v = zeros(nTime,np);
            for tempi = 1:nTime
                T_data_u(tempi,:) = RD.ResultDisp{ImgSeqNum-(2+nTime)+tempi, 1}.U(1:2:np*2)';
                T_data_v(tempi,:) = RD.ResultDisp{ImgSeqNum-(2+nTime)+tempi, 1}.U(2:2:np*2)';

            end
            nB = 3; t_train = (ImgSeqNum-1-nTime:ImgSeqNum-2)'; t_pre = (ImgSeqNum-1)';
            [u_pred,~,~,~] = funPOR_GPR(T_data_u,t_train,t_pre,nB);
            [v_pred,~,~,~] = funPOR_GPR(T_data_v,t_train,t_pre,nB);
            tempu = u_pred(1,:); tempv = v_pred(1,:);
            U0 = [tempu(:),tempv(:)]'; U0 = U0(:);

            % %%%%% After running the new ImgSeqNum, you can uncomment these
            % %%%%% lines to compare how the initial guess has been improved.
            % Plotdisp_show(U0-RD.ResultDisp{ImgSeqNum-1}.U,DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM(:,1:4));
            % Plotdisp_show(RD.ResultDisp{ImgSeqNum-2}.U-RD.ResultDisp{ImgSeqNum-1}.U,DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM(:,1:4));
        end
    end

    % ====== Compute f(X)-g(x+u) ======
    % PlotImgDiff(x0,y0,u,v,fNormalized,gNormalized);
    RD.ResultFEMeshEachFrame{ImgSeqNum-1} = struct( 'coordinatesFEM',DICmesh.coordinatesFEM,'elementsFEM',DICmesh.elementsFEM);
    fprintf('------------ Section 3 Done ------------ \n \n')

    %% Section 4: Subproblem 1 -or- Local ICGN Subset DIC
    fprintf('------------ Section 4 Start ------------ \n')
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This section is to solve the first local step in ALDIC: Subproblem 1
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % ====== ALStep 1 Subproblem1: Local Subset DIC ======
    mu=0; beta=0; tol=1e-3; ALSolveStep=1; ALSub1Time=zeros(6,1); ALSub2Time=zeros(6,1);
    ConvItPerEle=zeros(size(DICmesh.coordinatesFEM,1),6); ALSub1BadPtNum=zeros(6,1);
    disp(['***** Start step',num2str(ALSolveStep),' Subproblem1 *****'])
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ------ Start Local DIC IC-GN iteration ------
    [USubpb1,FSubpb1,HtempPar,ALSub1Timetemp,ConvItPerEletemp,LocalICGNBadPtNumtemp] = ...
        LocalICGN(U0,DICmesh.coordinatesFEM,Df,fNormalized,gNormalized,DICpara,'GaussNewton',tol,1);

    % Close next line for 2nd-order shape function (Line 323 as well)
    ALSub1Time(ALSolveStep) = ALSub1Timetemp;
    ConvItPerEle(:,ALSolveStep) = ConvItPerEletemp;
    ALSub1BadPtNum(ALSolveStep) = LocalICGNBadPtNumtemp;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ------  Manually find some bad points from Local Subset ICGN step ------
    % Comment these lines below if you don't need to manually remove local bad points
    % %%%%% Comment START %%%%%
    % close all; USubpb1World = USubpb1; USubpb1World(2:2:end) = -USubpb1(2:2:end);
    % %Plotuv(USubpb1World,DICmesh.x0,DICmesh.y0World);
    % disp('--- Start to manually remove bad points ---')
    % u = reshape(USubpb1(1:2:end),size(DICmesh.x0,1),size(DICmesh.x0,2));
    % v = reshape(USubpb1(2:2:end),size(DICmesh.x0,1),size(DICmesh.x0,2));
    % [u,v,~,Local_BadptRow,Local_BadptCol,RemoveOutliersList] = funRemoveOutliers(u',v',[],0.5,100); u=u';v=v';
    % USubpb1(1:2:end) = reshape(u,size(DICmesh.coordinatesFEM,1),1); USubpb1(2:2:end) = reshape(v,size(DICmesh.coordinatesFEM,1),1);
    % f11 = reshape(FSubpb1(1:4:end),size(DICmesh.x0,1),size(DICmesh.x0,2));
    % f21 = reshape(FSubpb1(2:4:end),size(DICmesh.x0,1),size(DICmesh.x0,2));
    % f12 = reshape(FSubpb1(3:4:end),size(DICmesh.x0,1),size(DICmesh.x0,2));
    % f22 = reshape(FSubpb1(4:4:end),size(DICmesh.x0,1),size(DICmesh.x0,2));
    % f11=f11'; f11(RemoveOutliersList) = NaN; f11 = inpaint_nans(f11,4); f11=f11';
    % f21=f21'; f21(RemoveOutliersList) = NaN; f21 = inpaint_nans(f21,4); f21=f21';
    % f12=f12'; f12(RemoveOutliersList) = NaN; f12 = inpaint_nans(f12,4); f12=f12';
    % f22=f22'; f22(RemoveOutliersList) = NaN; f22 = inpaint_nans(f22,4); f22=f22';
    % FSubpb1(1:4:end) = reshape(f11,size(DICmesh.coordinatesFEM,1),1);
    % FSubpb1(2:4:end) = reshape(f21,size(DICmesh.coordinatesFEM,1),1);
    % FSubpb1(3:4:end) = reshape(f12,size(DICmesh.coordinatesFEM,1),1);
    % FSubpb1(4:4:end) = reshape(f22,size(DICmesh.coordinatesFEM,1),1);
    % disp('--- Remove bad points done ---')
    % %%%%% Comment END %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % ------ Plot ------
    USubpb1World = USubpb1; USubpb1World(2:2:end) = -USubpb1(2:2:end);
    FSubpb1World = FSubpb1; FSubpb1World(2:4:end) = -FSubpb1World(2:4:end); FSubpb1World(3:4:end) = -FSubpb1World(3:4:end);
    close all; Plotuv(USubpb1World,DICmesh.x0,DICmesh.y0World);
    Plotdisp_show(USubpb1World,DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM,[],'NoEdgeColor');
    Plotstrain_show(FSubpb1World,DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM,[],'NoEdgeColor');


    save(['Subpb1_step',num2str(ALSolveStep)],'USubpb1','FSubpb1');
    % fprintf('------------ Section 4 Done ------------ \n \n')

    if UseGlobal
        %% Section 5: Subproblem 2 -- solve the global compatible displacement field
        fprintf('------------ Section 5 Start ------------ \n');
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % This section is to solve the global step in ALDIC Subproblem 2
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % ======= ALStep 1 Subproblem 2: Global constraint =======
        % ------ Smooth displacements for a better F ------
        DICpara.DispFilterSize=0; DICpara.DispFilterStd=0; DICpara.StrainFilterSize=0; DICpara.StrainFilterStd=0; LevelNo=1;
        FSubpb1 = funSmoothStrain(FSubpb1,DICmesh,DICpara);

        % ====== Define penalty parameter ======
        mu = 1e-3; udual = 0*FSubpb1; vdual = 0*USubpb1;
        betaList = [1e-3,sqrt(1e-5),1e-2,sqrt(1e-3),1e-1,sqrt(1e-1)]*mean(DICpara.winstepsize).^2.*mu; % Tune beta in the betaList
        Err1 = zeros(length(betaList),1); Err2 = Err1;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp(['***** Start step',num2str(ALSolveStep),' Subproblem2 *****'])
        % ====== Check to use FD or FE methods to solve Subpb2 step ======
        if DICpara.Subpb2FDOrFEM == 1 % Using FD method
            % ====== Build sparse finite difference operator ======
            disp('Assemble finite difference operator D');
            M = size(DICmesh.x0,1); N = size(DICmesh.x0,2);
            tic; Rad = 1; D = funDerivativeOp((M-2*Rad),(N-2*Rad),DICpara.winstepsize); % D = sparse(4*(M-2*Rad)*(N-2*Rad), 2*(M-2*Rad)*(N-2*Rad));
            D2 = funDerivativeOp(M,N,DICpara.winstepsize); toc
            disp('Finish assembling finite difference operator D');
            % ===== Solver using finite difference approximation ======
            tic; a = FSubpb1-udual; b = USubpb1-vdual;
            Rad = 1; [temp3,temp4] = funFDNeumannBCInd(size(DICmesh.coordinatesFEM,1),M,N,Rad); % Find coordinatesFEM that belong to (x(Rad+1:M-Rad,Rad+1:N-Rad),y(Rad+1:M-Rad,Rad+1:N-Rad))
            atemp = a(temp3); btemp = b(temp4); hbar = waitbar(0,'Please wait for Subproblem 2 global step!');
            for tempk = 1:length(betaList)
                beta = betaList(tempk);
                tempAMatrixSub2 = (beta*(D')*D) + mu*speye(2*(M-2*Rad)*(N-2*Rad));
                USubpb2temp = (tempAMatrixSub2) \ (beta*D'*atemp + mu*btemp ) ;
                USubpb2 = USubpb1; USubpb2(temp4) = USubpb2temp;
                FSubpb2 = D2*USubpb2;

                Err1(tempk) = norm(USubpb1-USubpb2,2);
                Err2(tempk) = norm(FSubpb1-FSubpb2,2);
                waitbar(tempk/(length(betaList)+1));
            end
            ErrSum = Err1+Err2*mean(DICpara.winstepsize)^2; [~,indexOfbeta] = min(ErrSum);

            try % Tune the best beta by a quadratic polynomial fitting
                [fitobj] = fit(log10(betaList(indexOfbeta-1:1:indexOfbeta+1))',ErrSum(indexOfbeta-1:1:indexOfbeta+1),'poly2');
                p = coeffvalues(fitobj); beta = 10^(-p(2)/2/p(1));
            catch, beta = betaList(indexOfbeta);
            end

            % Using the optimal beta to solve the ALDIC Subproblem 2 again
            tempAMatrixSub2 = (beta*(D')*D) + mu*speye(2*(M-2*Rad)*(N-2*Rad));
            USubpb2temp = (tempAMatrixSub2) \ (beta*D'*atemp + mu*btemp) ;
            USubpb2 = USubpb1; USubpb2(temp4) = USubpb2temp;
            waitbar(1); close(hbar);
            %%%%%%%%%%%%%% End of using finite difference approximation %%%%%%%%%%%%%%

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        else %Subpb2FDOrFEM: Using FE method
            M = size(DICmesh.x0,1); N = size(DICmesh.x0,2); GaussPtOrder = 2; alpha = 0;
            close all; hbar = waitbar(0,'Please wait for Subproblem 2 global step!');
            % ====== Solver using finite element method ======
            for tempk = 1:length(betaList)
                beta = betaList(tempk);
                GaussPtOrder = 2; alpha = 0; [USubpb2] = Subpb2(DICmesh,beta,mu,USubpb1,FSubpb1,udual,vdual,alpha,GaussPtOrder);
                [FSubpb2] = funGlobal_NodalStrainAvg(DICmesh.coordinatesFEM,DICmesh.elementsFEM,USubpb2,GaussPtOrder);
                Err1(tempk) = norm(USubpb1-USubpb2,2);
                Err2(tempk) = norm(FSubpb1-FSubpb2,2);
                waitbar(tempk/(length(betaList)+1));
            end
            Err1Norm = (Err1-mean(Err1))/std(Err1); figure, plot(Err1Norm);
            Err2Norm = (Err2-mean(Err2))/std(Err2); figure, plot(Err2Norm);
            ErrSum = Err1Norm+Err2Norm; figure,plot(ErrSum); [~,indexOfbeta] = min(ErrSum);

            try % Tune the best beta by a quadratic polynomial fitting
                [fitobj] = fit(log10(betaList(indexOfbeta-1:1:indexOfbeta+1))',ErrSum(indexOfbeta-1:1:indexOfbeta+1),'poly2');
                p = coeffvalues(fitobj); beta = 10^(-p(2)/2/p(1));
            catch, beta = betaList(indexOfbeta);
            end

            % Using the optimal beta to solve the ALDIC Subproblem 2 again
            [USubpb2] = Subpb2(DICmesh,beta,mu,USubpb1,FSubpb1,udual,vdual,alpha,GaussPtOrder);
            USubpb2 = full(USubpb2);
            waitbar(1); close(hbar);
        end
        ALSub2Time(ALSolveStep) = toc; toc
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % ------- Before computing strain, we smooth the displacement field -------
        % USubpb2 = funSmoothDisp(USubpb2,DICmesh,DICpara);
        % ------- Compute strain field --------
        if DICpara.Subpb2FDOrFEM == 1 %FD
            FSubpb2 = D2*USubpb2; % D2 = funDerivativeOp(M,N,winstepsize);
        else %FEM
            [FSubpb2] = funGlobal_NodalStrainAvg(DICmesh.coordinatesFEM,DICmesh.elementsFEM,USubpb2,GaussPtOrder);
        end

        % ------- Smooth strain field --------
        FSubpb2 = funSmoothStrain(FSubpb2,DICmesh,DICpara);

        % ------- Save data ------
        save(['Subpb2_step',num2str(ALSolveStep)],'USubpb2','FSubpb2');

        % ------ Plot ------
        % USubpb2World = USubpb2; USubpb2World(2:2:end) = -USubpb2(2:2:end);
        % FSubpb2World = FSubpb2; FSubpb2World(2:4:end) = -FSubpb2World(2:4:end); FSubpb2World(3:4:end) = -FSubpb2World(3:4:end);
        % close all; Plotuv(USubpb2World,DICmesh.x0,DICmesh.y0World);
        % Plotdisp_show(USubpb2World,DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM,[],'NoEdgeColor');
        % Plotstrain_show(FSubpb2World,DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM,[],'NoEdgeColor');
        % Plotstrain_show(FSubpb2,DICmesh.coordinatesFEM ,DICmesh.elementsFEM,[],'NoEdgeColor');

        % convtimes = [ConvItPerEle(:,2),ConvItPerEle(:,2)]'; convtimes=convtimes(:);
        % Plotdisp_show(convtimes,DICmesh.coordinatesFEM,DICmesh.elementsFEM,[],'NoEdgeColor');


        % ======= Update dual variables =======
        if DICpara.Subpb2FDOrFEM == 1 %FD
            udualtemp1 = (FSubpb2 - FSubpb1); udualtemp2 = udualtemp1(temp3);
            vdualtemp1 = (USubpb2 - USubpb1); vdualtemp2 = vdualtemp1(temp4);
            udual = zeros(4*M*N,1); vdual = zeros(2*M*N,1);
            udual(temp3) = udualtemp2; vdual(temp4) = vdualtemp2;
        else  % FEM or other methods
            udual = FSubpb2 - FSubpb1; vdual = USubpb2 - USubpb1;
        end
        save(['uvdual_step',num2str(ALSolveStep)],'udual','vdual');
        fprintf('------------ Section 5 Done ------------ \n \n')


        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Section 6: ADMM iterations
        fprintf('------------ Section 6 Start ------------ \n')
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % This section is to run ADMM iteration: Subproblem 1 & 2
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % ==================== ADMM AL Loop ==========================
        ALSolveStep = 1; tol2 = 1e-4; UpdateY = 1e4; CrackOrNot = 0; CrackPath1 = [0,0]; CrackPath2 = [0,0]; CrackTip = [0,0];
        HPar = cell(21,1); for tempj = 1:21, HPar{tempj} = HtempPar(:,tempj); end

        while (ALSolveStep < 3)
            ALSolveStep = ALSolveStep + 1;  % Update using the last step
            %%%%%%%%%%%%%%%%%%%%%%% Subproblem 1 %%%%%%%%%%%%%%%%%%%%%%%%%
            disp(['***** Start step',num2str(ALSolveStep),' Subproblem1 *****']);
            tic;[USubpb1,~,ALSub1Timetemp,ConvItPerEletemp,LocalICGNBadPtNumtemp] = Subpb1( ...
                USubpb2,FSubpb2,udual,vdual,DICmesh.coordinatesFEM,...
                Df,fNormalized,gNormalized,mu,beta,HPar,ALSolveStep,DICpara,'GaussNewton',tol);

            FSubpb1 = FSubpb2; toc
            ALSub1Time(ALSolveStep) = ALSub1Timetemp; ConvItPerEle(:,ALSolveStep) = ConvItPerEletemp;  ALSub1BadPtNum(ALSolveStep) = LocalICGNBadPtNumtemp;

            % [row1,~] = find(ConvItPerEletemp(:)<0);
            % [row2,~] = find(ConvItPerEletemp(:)>99);
            % LocalICGNBadPt = unique(union(row1,row2));
            % USubpb1(2*LocalICGNBadPt-1) = nan;
            % USubpb1(2*LocalICGNBadPt) = nan;
            % Plotdisp_show(USubpb1,DICmesh.coordinatesFEM,DICmesh.elementsFEM,[],'NoEdgeColor');

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % ------  Manually find some bad points from Local Subset ICGN step ------
            % disp('--- Start to manually remove bad points --- \n')
            % disp('    Comment codes here if you do not have bad local points \n')
            % %%%%% Comment START %%%%%
            %   close all; Plotuv(USubpb1,DICmesh.x0,DICmesh.y0);
            %   u = reshape(USubpb1(1:2:end),M,N); v = reshape(USubpb1(2:2:end),M,N);
            %   [u,v,~,Subpb1_BadptRow,Subpb1_BadptCol] = funRemoveOutliers(u',v',[],0.5,100,Local_BadptRow,Local_BadptCol); u=u';v=v';
            %   disp('--- Remove bad points done ---')
            %   USubpb1(1:2:end) = reshape(u,size(DICmesh.coordinatesFEM,1),1);
            %   USubpb1(2:2:end) = reshape(v,size(DICmesh.coordinatesFEM,1),1);
            %   close all; Plotuv(USubpb1,DICmesh.x0,DICmesh.y0); Plotdisp_show(USubpb1,DICmesh.coordinatesFEM,DICmesh.elementsFEM);
            % %%%%% Comment END %%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            save(['Subpb1_step',num2str(ALSolveStep)],'USubpb1','FSubpb1');
            USubpb1 = funSmoothDisp(USubpb1,DICmesh,DICpara);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % ============== Subproblem 2 ==============
            disp(['***** Start step',num2str(ALSolveStep),' Subproblem2 *****'])
            if DICpara.Subpb2FDOrFEM == 1 %FD
                % ------- using finite difference approximation --------
                tic; a = FSubpb1-udual; b = USubpb1-vdual; atemp = a(temp3); btemp = b(temp4);
                USubpb2temp = (tempAMatrixSub2) \ (beta*D'*atemp + mu*btemp ) ;
                USubpb2 = USubpb1; USubpb2(temp4) = USubpb2temp; %toc
                % ------- End of using finite difference approximation --------
            else % FEM
                tic; [USubpb2] = Subpb2(DICmesh,beta,mu,USubpb1,FSubpb1,udual,vdual,alpha,GaussPtOrder);
                USubpb2 = full(USubpb2);
            end
            ALSub2Time(ALSolveStep) = toc; toc
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % ------- Before computing strain, we smooth the displacement field -------
            USubpb2 = funSmoothDisp(USubpb2,DICmesh,DICpara);
            % ------- Compute strain field --------
            if DICpara.Subpb2FDOrFEM == 1 %FD
                FSubpb2 = D2*USubpb2; % D2 = funDerivativeOp(M,N,winstepsize);
            else %FEM
                GaussPtOrder = 2; [FSubpb2] = funGlobal_NodalStrainAvg(DICmesh.coordinatesFEM,DICmesh.elementsFEM,USubpb2,GaussPtOrder);
            end
            % ------- Smooth strain field --------
            FSubpb2 = funSmoothStrain(FSubpb2,DICmesh,DICpara);

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
            if DICpara.Subpb2FDOrFEM == 1 %FD
                udualtemp1 =  (FSubpb2 - FSubpb1); udualtemp2 = udualtemp1(temp3);
                vdualtemp1 =  (USubpb2 - USubpb1); vdualtemp2 = vdualtemp1(temp4);
                udual(temp3) = udual(temp3) + udualtemp2;
                vdual(temp4) = vdual(temp4) + vdualtemp2;
            else %FEM
                udual = FSubpb2 - FSubpb1; vdual = USubpb2 - USubpb1;
            end

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

    %% Step1: Save data
    % disp_inc1 is the disp of nodes based on the renewed ROI
    % disp_inc2 is the disp of nodes based on the interpolated ROI
    % coor_inc_ref is the coor of nodes based on the ref(first) ROI

    if UseGlobal
        RD.ResultDisp_inc1{ImgSeqNum-1,1}.U = full(USubpb2);
        RD.ResultDisp_inc1{ImgSeqNum-1,1}.ALSub1BadPtNum = ALSub1BadPtNum;
        RD.ResultDisp_inc1{ImgSeqNum-1,1}.F = full(FSubpb2); % tempFoamAL;
        RD.DICmesh{ImgSeqNum-1} = DICmesh;
    else % regardless subproblem2
        RD.ResultDisp_inc1{ImgSeqNum-1}.U = full(USubpb1);
        RD.ResultDisp_inc1{ImgSeqNum-1}.ALSub1BadPtNum = ALSub1BadPtNum;
        RD.ResultDisp_inc1{ImgSeqNum-1}.F = full(FSubpb1); % tempFoamAL;
        RD.DICmesh{ImgSeqNum-1} = DICmesh;
    end

    %% Step2: interpolation due to Reference updating
    switch camera0OrNot
        case 'camera0'
            % Build the coordinates in every frame
            RD.Coordinate_inc_ref{1,1}.Coor = RD.ResultFEMeshEachFrame{1,1}.coordinatesFEM;

            % Find the disp_inc in the ref config.
            if ImgSeqNum == 2
                RD.ResultDisp_inc2{1,1}.U = [RD.ResultDisp_inc1{1}.U(1:2:end),RD.ResultDisp_inc1{1}.U(2:2:end)];

            else
                tempX = RD.ResultDisp_inc1{ImgSeqNum-1, 1}.U(1:2:end);
                tempY = RD.ResultDisp_inc1{ImgSeqNum-1, 1}.U(2:2:end);

                rbfInterpX = rbfcreate(RD.ResultFEMeshEachFrame{ImgSeqNum-1, 1}.coordinatesFEM',tempX');
                rbfInterpY = rbfcreate(RD.ResultFEMeshEachFrame{ImgSeqNum-1, 1}.coordinatesFEM',tempY');

                tempResultDisp(1,:) = rbfinterp(RD.Coordinate_inc_ref{ImgSeqNum-1,1}.Coor',rbfInterpX);
                tempResultDisp(2,:) = rbfinterp(RD.Coordinate_inc_ref{ImgSeqNum-1,1}.Coor',rbfInterpY);
                RD.ResultDisp_inc2{ImgSeqNum-1,1}.U =  tempResultDisp'; 
            end

            % Build the coordinates in every frame
            RD.Coordinate_inc_ref{ImgSeqNum,1}.Coor = RD.Coordinate_inc_ref{ImgSeqNum-1,1}.Coor + RD.ResultDisp_inc2{ImgSeqNum-1,1}.U;

           

        case 'notCamera0'

            % Build the coordinates in every frame 
            RD.Coordinate_inc_ref{1,1}.Coor = RD.Coordinates_corr; 

            % Find the disp_inc in the ref config.
            tempX = RD.ResultDisp_inc1{ImgSeqNum-1, 1}.U(1:2:end);
            tempY = RD.ResultDisp_inc1{ImgSeqNum-1, 1}.U(2:2:end);

            rbfInterpX = rbfcreate(RD.ResultFEMeshEachFrame{ImgSeqNum-1, 1}.coordinatesFEM',tempX');
            rbfInterpY = rbfcreate(RD.ResultFEMeshEachFrame{ImgSeqNum-1, 1}.coordinatesFEM',tempY');

            tempResultDisp(1,:) = rbfinterp(RD.Coordinate_inc_ref{ImgSeqNum-1,1}.Coor',rbfInterpX);
            tempResultDisp(2,:) = rbfinterp(RD.Coordinate_inc_ref{ImgSeqNum-1,1}.Coor',rbfInterpY);
            RD.ResultDisp_inc2{ImgSeqNum-1,1}.U = tempResultDisp'; 
           
            % Build the coordinates in every frame
            RD.Coordinate_inc_ref{ImgSeqNum,1}.Coor = RD.Coordinate_inc_ref{ImgSeqNum-1,1}.Coor + RD.ResultDisp_inc2{ImgSeqNum-1,1}.U;

    end




    % if 0
    %     % ------ Plot ------
    %     if DICpara.DICIncOrNot == 1
    %         U_inc = reshape(RD.ResultDisp_inc{ImgSeqNum,1}.',[],1);
    %         close all; Plotuv(U_inc,DICmesh_init.x0,DICmesh_init.y0World); %%%%%%%%更新不更新roi，coorEach变，最终位移不变。
    %         Plotdisp_show(U_inc,DICmesh_init.coordinatesFEMWorld,DICmesh_init.elementsFEM);
    %     else
    %         USubpb2World = USubpb2; USubpb2World(2:2:end) = -USubpb2(2:2:end);
    %         FSubpb2World = FSubpb2; FSubpb2World(2:4:end) = -FSubpb2World(2:4:end); FSubpb2World(3:4:end) = -FSubpb2World(3:4:end);
    %         close all; Plotuv(USubpb2World,DICmesh_init.x0,DICmesh_init.y0World);
    %         Plotdisp_show(USubpb2World,DICmesh_init.coordinatesFEMWorld,DICmesh_init.elementsFEM);
    %     end
    % else
    % end
    % ------ Save results ------
    % Find img name and save all the results

    % [~,imgname,imgext] = fileparts(file_name{1,end});
    % results_name = ['results_',imgname,'_ws',num2str(DICpara.winsize),'_st',num2str(DICpara.winstepsize),'.mat'];
    % save(results_name, 'file_name','DICpara','DICmesh','RD','ALSub1Time','ALSub2Time','ALSolveStep');

end



end






