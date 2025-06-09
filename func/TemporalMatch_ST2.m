function RD = TemporalMatch_ST2(DICpara,file_name,ImgNormalized,RD,camera0OrNot,varargin)

% debug
UseGlobal = 0;
showImgOrNot = 1;

% Init
if strcmp(camera0OrNot,'notCamera0')
    ImgNormalized_reference = varargin;
    ImgNormalized = [ImgNormalized_reference,ImgNormalized]; % insert the ref. image into the first row
end

fNormalized = ImgNormalized{1};
Df = funImgGradient(fNormalized,fNormalized);

for ImgSeqNum = 2 : length(ImgNormalized)

    disp(['Current image frame #: ', num2str(ImgSeqNum),'/',num2str(length(ImgNormalized))]);
    gNormalized = ImgNormalized{ImgSeqNum}; % Load current deformed image frame

    %% Section 3: Compute an initial guess of the unknown displacement field
    fprintf('\n'); fprintf('------------ Section 3 Start ------------ \n')
    %%%%% One practical strategy is to let first 7 frames do the FFT-based
    %%%%% cross correlation and then using the data driven method to estimate
    %%%%% initial guesses for other frames
    if ImgSeqNum == 2
        DICpara.NewFFTSearch = 1; DICpara.InitFFTSearchMethod = [];
    elseif ImgSeqNum < 7          %  Zach revised!!!!!!!!!!!!!!!!!!
        DICpara.NewFFTSearch = 1; % Use FFT-based cross correlation to compute the initial guess
    else
        DICpara.NewFFTSearch = 0; % Apply data driven method to estimate initial guesses for later frames
    end

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
    mu=0; beta=0; tol=1e-2; ALSolveStep=1; ALSub1Time=zeros(6,1); ALSub2Time=zeros(6,1);
    ConvItPerEle=zeros(size(DICmesh.coordinatesFEM,1),6); ALSub1BadPtNum=zeros(6,1);
    disp(['***** Start step',num2str(ALSolveStep),' Subproblem1 *****'])
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ------ Start Local DIC IC-GN iteration ------
    [USubpb1,FSubpb1,HtempPar,ALSub1Timetemp,ConvItPerEletemp,LocalICGNBadPtNumtemp] = ...
        LocalICGN(U0,DICmesh.coordinatesFEM,Df,fNormalized,gNormalized,DICpara,'GaussNewton',tol,1);

    % Close next line for 2nd-order shape function (Line 323 as well)
    ALSub1Time(ALSolveStep) = ALSub1Timetemp; ConvItPerEle(:,ALSolveStep) = ConvItPerEletemp; ALSub1BadPtNum(ALSolveStep) = LocalICGNBadPtNumtemp; toc

    save(['Subpb1_step',num2str(ALSolveStep)],'USubpb1','FSubpb1');
    % fprintf('------------ Section 4 Done ------------ \n \n')

    if UseGlobal
        %% Section 5: Subproblem 2 -- solve the global compatible displacement field
        fprintf('------------ Section 5 Start ------------ \n'); tic;
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
        %USubpb2 = funSmoothDisp(USubpb2,coordinatesFEM,elementsFEM,x0,y0,winstepsize,DispFilterSize,DispFilterStd);
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
        USubpb2World = USubpb2; USubpb2World(2:2:end) = -USubpb2(2:2:end);
        FSubpb2World = FSubpb2; FSubpb2World(2:4:end) = -FSubpb2World(2:4:end); FSubpb2World(3:4:end) = -FSubpb2World(3:4:end);
        % close all; Plotuv(USubpb2World,DICmesh.x0,DICmesh.y0World);
        % Plotdisp_show(USubpb2World,DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM,[],'NoEdgeColor');
        % Plotstrain_show(FSubpb2World,DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM,[],'NoEdgeColor');
        % Plotstrain_show(FSubpb2,DICmesh.coordinatesFEM ,DICmesh.elementsFEM,[],'NoEdgeColor');

        % ======= Update dual variables =======
        if DICpara.Subpb2FDOrFEM == 1 %FD
            udualtemp1 = (FSubpb2 - FSubpb1); udualtemp2 = udualtemp1(temp3);
            vdualtemp1 = (USubpb2 - USubpb1); vdualtemp2 = vdualtemp1(temp4);
            udual = zeros(4*M*N,1); vdual = zeros(2*M*N,1);
            udual(temp3) = udualtemp2; vdual(temp4) = vdualtemp2;

            FSubpb2(temp3) = FSubpb1(temp3);
            USubpb2(temp4) = USubpb1(temp4);

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
            % USubpb1 = funSmoothDisp(USubpb1,DICmesh,DICpara);

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
            % USubpb2 = funSmoothDisp(USubpb2,DICmesh,DICpara);
            % ------- Compute strain field --------
            if DICpara.Subpb2FDOrFEM == 1 %FD
                FSubpb2 = D2*USubpb2; % D2 = funDerivativeOp(M,N,winstepsize);
            else %FEM
                GaussPtOrder = 2; [FSubpb2] = funGlobal_NodalStrainAvg(DICmesh.coordinatesFEM,DICmesh.elementsFEM,USubpb2,GaussPtOrder);
            end
            % ------- Smooth strain field --------
            % FSubpb2 = funSmoothStrain(FSubpb2,DICmesh,DICpara);

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

                FSubpb2(temp3) = FSubpb1(temp3);
                USubpb2(temp4) = USubpb1(temp4);

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

    %% Save data
    if strcmp(camera0OrNot,'notCamera0')

        if UseGlobal
            RD.ResultDisp{ImgSeqNum-1}.U = full(USubpb2);
            RD.ResultDefGrad{ImgSeqNum-1}.F = full(FSubpb2); % tempFoamAL;
        else % regardless subproblem2
            RD.ResultDisp{ImgSeqNum-1}.U = full(USubpb1);
            RD.ResultDefGrad{ImgSeqNum-1}.F = full(FSubpb1); % tempFoamAL;
        end
        RD.ResultDisp{ImgSeqNum-1}.ALSub1BadPtNum = ALSub1BadPtNum;
        RD.DICmesh{ImgSeqNum-1} = DICmesh;

    else

        if UseGlobal
            RD.ResultDisp{ImgSeqNum-1}.U = full(USubpb2);
            RD.ResultDefGrad{ImgSeqNum-1}.F = full(FSubpb2); % tempFoamAL;
        else % regardless subproblem2
            RD.ResultDisp{ImgSeqNum-1}.U = full(USubpb1);
            RD.ResultDefGrad{ImgSeqNum-1}.F = full(FSubpb1); % tempFoamAL;
        end
        RD.ResultDisp{ImgSeqNum-1}.ALSub1BadPtNum = ALSub1BadPtNum;
        RD.DICmesh{ImgSeqNum-1} = DICmesh;

    end
    
    if showImgOrNot == 1
        if UseGlobal
            % ------ Plot ------
            USubpb2World = USubpb2; USubpb2World(2:2:end) = -USubpb2(2:2:end);
            FSubpb2World = FSubpb2; FSubpb2World(2:4:end) = -FSubpb2World(2:4:end); FSubpb2World(3:4:end) = -FSubpb2World(3:4:end);
            %close all;
            Plotuv(USubpb2World,DICmesh.x0,DICmesh.y0World);
            Plotdisp_show(USubpb2World,DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM);
            %Plotstrain_show(FSubpb2World,DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM);
        else
            % ------ Plot ------
            USubpb2World = USubpb1; USubpb2World(2:2:end) = -USubpb1(2:2:end);
            FSubpb2World = FSubpb1; FSubpb2World(2:4:end) = -FSubpb2World(2:4:end); FSubpb2World(3:4:end) = -FSubpb2World(3:4:end);
            close all; Plotuv(USubpb2World,DICmesh.x0,DICmesh.y0World);
            Plotdisp_show(USubpb2World,DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM);
            %Plotstrain_show(FSubpb2World,DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM);
        end
    else
    end

    % ------ Save results ------
    % Find img name and save all the results

    % [~,imgname,imgext] = fileparts(file_name{1,end});
    % results_name = ['results_',imgname,'_ws',num2str(DICpara.winsize),'_st',num2str(DICpara.winstepsize),'.mat'];
    % save(results_name, 'file_name','DICpara','DICmesh','RD','ALSub1Time','ALSub2Time','ALSolveStep');
end
end






