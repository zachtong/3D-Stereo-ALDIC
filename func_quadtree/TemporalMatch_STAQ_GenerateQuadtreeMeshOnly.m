function TemporalMatch_STAQ_GenerateQuadtreeMeshOnly(ImgSeqNum,DICpara,file_name,ImgMask,ImgNormalized,RD,StereoInfo,camera0OrNot,shapeFuncOrder)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This part is designed for temporal DIC, which means dealing with left and
% right images separately using 2D-ALDIC.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% debug options
UseGlobal = 1;
DICpara.showImgOrNot = 0;

% Init
imageNum = length(ImgNormalized);
incOrNot = DICpara.DICIncOrNot; % 0: acc, 1: inc

switch camera0OrNot
    case 'camera0'

    case 'notCamera0'
        RD.Coordinates_corr = StereoInfo.ResultFEMesh_corr; % Backup: The first picture from the right camera corresponds
        % to the picture from the left camera and needs to be interpolated back.
end


%% For loop
    disp(['Current image frame #: ', num2str(ImgSeqNum),'/',num2str(length(ImgNormalized))]);

    if incOrNot == 1
        fNormalizedMask = double(ImgMask{ImgSeqNum-1} ) ; % Load the mask file of previous frame
        gNormalizedMask = double(ImgMask{ImgSeqNum}); % Load the mask file of current frame
    else
        fNormalizedMask = double(ImgMask{1} ) ; % Load the mask file of previous frame
        gNormalizedMask = ones(size(ImgNormalized{1}));
    end


    fNormalized = ImgNormalized{ImgSeqNum-1}.* fNormalizedMask; % Load previous reference image frame
    gNormalized = ImgNormalized{ImgSeqNum}.* gNormalizedMask; % Load current deformed image frame

    Df = funImgGradient(fNormalized,fNormalized,fNormalizedMask);

    DICpara.ImgRefMask = fNormalizedMask;

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
    % DICpara.NewFFTSearch = 1; % DICpara.InitFFTSearchMethod = [];

        % ====== Integer Search ======
        %DICpara.InitFFTSearchMethod = 1;
        %[DICpara,x0temp,y0temp,u,v,cc]= IntegerSearchQuadtree(fNormalized,gNormalized,file_name,DICpara,showImgOrNot,ImgSeqNum);
        DICpara.fixSearchDistanceOrNot = 1;
        % %%%%% Optional codes to measure more gridded measurement points %%%%%
        [DICpara,x0temp_f,y0temp_f,u_f,v_f,cc]= IntegerSearch(fNormalized,gNormalized,file_name,DICpara,ImgSeqNum);

        xnodes = 1+0.5*DICpara.winsize : DICpara.winstepsize : DICpara.ImgSize(1);
        ynodes = 1+0.5*DICpara.winsize : DICpara.winstepsize : DICpara.ImgSize(2);

        [x0temp,y0temp] = ndgrid(xnodes,ynodes);   u_f_NotNanInd = find(~isnan(u_f(:)));

        op1 = rbfcreate( [x0temp_f(u_f_NotNanInd),y0temp_f(u_f_NotNanInd)]',[u_f(u_f_NotNanInd)]','RBFFunction', 'thinplate');
        rbfcheck_maxdiff = rbfcheck(op1); % Check: rbf thin-plate interpolation
        if rbfcheck_maxdiff > 1e-3, disp('Please check rbf interpolation! Pause here.'); pause; end
        u = rbfinterp([x0temp(:),y0temp(:)]', op1 );

        op2 = rbfcreate( [x0temp_f(u_f_NotNanInd),y0temp_f(u_f_NotNanInd)]',[v_f(u_f_NotNanInd)]','RBFFunction', 'thinplate');
        rbfcheck_maxdiff = rbfcheck(op2); % Check: rbf thin-plate interpolation
        if rbfcheck_maxdiff > 1e-3, disp('Please check rbf interpolation! Pause here.'); pause; end
        v = rbfinterp([x0temp(:),y0temp(:)]', op2 );

        u = regularizeNd([x0temp(:),y0temp(:)],u(:),{xnodes',ynodes'},1e-3);
        v = regularizeNd([x0temp(:),y0temp(:)],v(:),{xnodes',ynodes'},1e-3);

        % ====== FEM mesh set up ======
        [DICmesh_nonQuadtree] = MeshSetUp(x0temp,y0temp,DICpara); %clear x0temp y0temp;

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

        DICmesh_nonQuadtree.elementMinSize = DICpara.winsizeMin; % min element size in the refined quadtree mesh

        GenerateQuadtreeMesh(U0,Df,fNormalizedMask,DICmesh_nonQuadtree,DICpara); % Generate the quadtree mesh
                                                        
    fprintf('------------ Section 3 Done ------------ \n \n')
   
end









