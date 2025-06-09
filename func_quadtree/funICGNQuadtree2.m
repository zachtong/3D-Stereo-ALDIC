function [U,F,stepwithinwhile,HGlobal] = funICGNQuadtree2(U0,x0,y0,Df,ImgRef,ImgDef,winsize,tol,ICGNmethod)
%FUNCTION [U,F,stepwithinwhile,HGlobal] = funICGN(U0,x0,y0,Df,ImgRef,ImgDef,winsize,tol,method)
% The Local ICGN subset solver (part II): ICGN iteration 
% (see part I: ./func/LocalICGN.m)
% ----------------------------------------------
%   INPUT: U0                   Initial guess of the displacement fields
%          x0,y0                FE mesh nodal coordinates
%          Df                   Image grayscale value gradients
%          ImgRef               Reference image
%          ImgDef               Deformed image
%          winsize              DIC parameter subset size  
%          ICGNmethod           ICGN iteration scheme: 'GaussNewton' -or- 'LevenbergMarquardt'
%          tol                  ICGN iteration stopping threshold
%
%   OUTPUT: U                   Disp vector: [Ux_node1, Uy_node1, ... , Ux_nodeN, Uy_nodeN]';
%           F                   Deformation gradient tensor
%                               F = [F11_node1, F21_node1, F12_node1, F22_node1, ... , F11_nodeN, F21_nodeN, F12_nodeN, F22_nodeN]';
%           stepwithinwhile     ICGN iteration step for convergence
%           HGlobal             Hessian matrix for each local subset
%
% ----------------------------------------------
% Author: Jin Yang.  
% Contact and support: jyang526@wisc.edu -or- aldicdvc@gmail.com
% Last time updated: 02/2020.
% ==============================================

%% Initialization
warning('off');
DfCropWidth = Df.DfCropWidth;
imgSize = Df.imgSize;
winsize0 = winsize;


%% ---------------------------
% Find local subset region
x = [x0-winsize/2 ; x0+winsize/2 ; x0+winsize/2 ; x0-winsize/2];  % [coordinates(elements(j,:),1)];
y = [y0-winsize/2 ; y0+winsize/2 ; y0+winsize/2 ; y0-winsize/2];  % [coordinates(elements(j,:),2)];

% ---------------------------
% Initialization: Get P0
P0 = [U0(1) zeros(1,5), U0(2), zeros(1,5)]';
P = P0;
 

% ---------------------------
% Find region for f
[XX,YY] = ndgrid([x(1):1:x(3)], [y(1):1:y(3)]);
%tempf = imgfNormalizedbc.eval(XX,YY); 
%DfDx = imgfNormalizedbc.eval_Dx(XX,YY);
%DfDy = imgfNormalizedbc.eval_Dy(XX,YY);

%%%%%%%%%%%% !!!Mask START %%%%%%%%%%%%%%
tempfImgMask = Df.ImgRefMask([x(1):1:x(3)], [y(1):1:y(3)]);
tempf = ImgRef([x(1):1:x(3)],[y(1):1:y(3)]) .* tempfImgMask;
%%%%%%%%%%%% !!!Mask END %%%%%%%%%%%%%%

DfDx = Df.DfDx((x(1)-DfCropWidth):1:(x(3)-DfCropWidth), (y(1)-DfCropWidth):1:(y(3)-DfCropWidth));
DfDy = Df.DfDy((x(1)-DfCropWidth):1:(x(3)-DfCropWidth), (y(1)-DfCropWidth):1:(y(3)-DfCropWidth));

%% %%%%%%%% If there are >60% of the subset are painted with patterns %%%%%%%%%%%%
maxIterNum = 100;

%%%%%%%%%%%% !!!Mask START %%%%%%%%%%%%%%
DfDxImgMaskIndCount = sum(double(1-logical(tempf(:))));
% [DfDxImgMaskIndRow,~] = find(abs(tempf)<1e-10);
% DfDxImgMaskIndCount = length(DfDxImgMaskIndRow);
%%%%%%%%%%%% !!!Mask END %%%%%%%%%%%%%%
  
if DfDxImgMaskIndCount < 0.4*(winsize+1)^2
     
    if DfDxImgMaskIndCount > 0.0*(winsize+1)^2 % For those subsets where are 0's in the image mask file
         
        winsize = 2*max(ceil(sqrt((winsize+1)^2/((winsize+1)^2-DfDxImgMaskIndCount))*winsize/2)); % Increase the subset size a bit to guarantuee there there enough pixels
        x = [x0-winsize/2 ; x0+winsize/2 ; x0+winsize/2 ; x0-winsize/2]; % Update x
        y = [y0-winsize/2 ; y0+winsize/2 ; y0+winsize/2 ; y0-winsize/2]; % Update y
        [XX,YY] = ndgrid([x(1):1:x(3)],[y(1):1:y(3)]); 
        tempfImgMask = Df.ImgRefMask([x(1):1:x(3)],[y(1):1:y(3)]);
        tempf = ImgRef([x(1):1:x(3)],[y(1):1:y(3)]) .* tempfImgMask;
        DfDx = Df.DfDx((x(1)-DfCropWidth):1:(x(3)-DfCropWidth), (y(1)-DfCropWidth):1:(y(3)-DfCropWidth));
        DfDy = Df.DfDy((x(1)-DfCropWidth):1:(x(3)-DfCropWidth), (y(1)-DfCropWidth):1:(y(3)-DfCropWidth));
         
    end
    
    %%%%%%%%%% !!!Mask: START %%%%%%%%%%%%
    %%%%% Find connected region to deal with possible continuities %%%%%
    tempf_BW2 = bwselect(logical(tempfImgMask), floor((winsize+1)/2), floor((winsize+1)/2), 4 );
    DfDx_BW2 = bwselect(logical(tempfImgMask), floor((winsize+1)/2), floor((winsize+1)/2), 4 );
    DfDy_BW2 = bwselect(logical(tempfImgMask), floor((winsize+1)/2), floor((winsize+1)/2), 4 );
    tempf = tempf .* double(tempf_BW2);
    DfDx = DfDx .* double(tempf_BW2);
    DfDy = DfDy .* double(tempf_BW2);
    %%%%%%%%%% !!!Mask: END %%%%%%%%%%%%
    
    dX = XX(:)-x0; dY = YY(:)-y0;

    DfDx_column = DfDx(:);
    DfDy_column = DfDy(:);

    dFdWdP = [DfDx_column, DfDx_column.*dX, DfDx_column.*dY, (DfDx_column.*dX.^2)/2, DfDx_column.*dX.*dY,...
    (DfDx_column.*dY.^2)/2, DfDy_column, DfDy_column.*dX, DfDy_column.*dY, (DfDy_column.*dX.^2)/2, ...
    DfDy_column.*dX.*dY, (DfDy_column.*dY.^2)/2 ];

    H = dFdWdP'*dFdWdP;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%%% Old codes: to compute H matrix by a for loop %%%%
    % tempCoordx = XX(:); tempCoordy = YY(:);
    % for tempij = 1:size(tempCoordx,1)
    %
    %         H = H + ([DfDx(tempCoordx(tempij)-DfCropWidth,tempCoordy(tempij)-DfCropWidth) DfDy(tempCoordx(tempij)-DfCropWidth,tempCoordy(tempij)-DfCropWidth)]*...
    %             [tempCoordx(tempij)-x0 0 tempCoordy(tempij)-y0 0 1 0; 0 tempCoordx(tempij)-x0 0 tempCoordy(tempij)-y0 0 1])'* ...
    %             ([DfDx(tempCoordx(tempij)-DfCropWidth,tempCoordy(tempij)-DfCropWidth) DfDy(tempCoordx(tempij)-DfCropWidth,tempCoordy(tempij)-DfCropWidth)]*...
    %             [tempCoordx(tempij)-x0 0 tempCoordy(tempij)-y0 0 1 0; 0 tempCoordx(tempij)-x0 0 tempCoordy(tempij)-y0 0 1]);
    %
    % end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%% !!!Mask: START %%%%%%%%%%%%
    meanf = mean(tempf(abs(tempf)>1e-10));
    bottomf = sqrt((length(tempf(abs(tempf)>1e-10))-1)*var(tempf(abs(tempf)>1e-10)));
    %%%%%%%%%% !!!Mask: END %%%%%%%%%%%%
    
    % --------------------------
    % Initialize while loop
    normOfWOld=2; normOfWNew=1; normOfWNewAbs=1; stepwithinwhile=0;
    switch ICGNmethod   % For Gauss-Newton method
        case 'LevenbergMarquardt'
            delta = 0.001; % For Levenberg-Marquardt method
            KappaOld=1e10; KappaNew=1e10; KappaStore=zeros(10,1); PStore=zeros(10,6);
        otherwise % 'GaussNewton'
            delta = 0;
    end
    
    while( (stepwithinwhile<maxIterNum) && (normOfWNew>tol) && (normOfWNewAbs>tol) )
        
        stepwithinwhile = stepwithinwhile+1;
        
        if stepwithinwhile>1 && DfDxImgMaskIndCount>0.0*(winsize0+1)^2
            winsize = 2*max(ceil(sqrt((winsize0+1)^2/(sum(double(tempg_BW2(:)))))*winsize/2)); % Increase the subset size a bit to guarantuee there there enough pixels
            x = [x0-winsize/2 ; x0+winsize/2 ; x0+winsize/2 ; x0-winsize/2]; % Update x
            y = [y0-winsize/2 ; y0+winsize/2 ; y0+winsize/2 ; y0-winsize/2]; % Update y
            [XX,YY] = ndgrid([x(1):1:x(3)],[y(1):1:y(3)]);
            tempfImgMask = Df.ImgRefMask([x(1):1:x(3)],[y(1):1:y(3)]);
        end
        
        %%%%%%%% Find region for g 
        % %[tempCoordy, tempCoordx] = meshgrid(y(1):y(3),x(1):x(3));
        tempCoordxMat = XX - x0*ones(winsize+1,winsize+1);
        tempCoordyMat = YY - y0*ones(winsize+1,winsize+1);

        % u22 and v22 are current x and y coordiantes
        u22 = 0.5*P(4)*tempCoordxMat.^2 + P(5)*tempCoordxMat.*tempCoordyMat + 0.5*P(6)*tempCoordyMat.^2 + ...
            (1+P(2))*tempCoordxMat + P(3)*tempCoordyMat + (x0+P(1))*ones(winsize+1,winsize+1);

        v22 = 0.5*P(10)*tempCoordxMat.^2 + P(11)*tempCoordxMat.*tempCoordyMat + 0.5*P(12)*tempCoordyMat.^2 + ...
            (P(8))*tempCoordxMat + (1+P(9))*tempCoordyMat + (y0+P(7))*ones(winsize+1,winsize+1);
        
        row1 = find(u22<3); row2 = find(u22>imgSize(1)-2); row3 = find(v22<3); row4 = find(v22>imgSize(2)-2);
        if ~isempty([row1; row2; row3; row4])
            normOfWNew = 1e6; % warning('Out of image boundary!')
            break;
        else
            
            %tempg = imggNormalizedbc.eval(u22,v22)
            %tempg = ba_interp2(ImgDef, v22, u22, 'cubic');
            
            % BicubicBspline interpolation
            tempg = ba_interp2_spline(ImgDef, v22, u22, 'cubicspline');

            DgDxImgMaskIndCount = sum(double(1-logical(tempg(:))));
            
            
            %%%%%%%%%% !!!Mask: START %%%%%%%%%%%%
            if DfDxImgMaskIndCount>0.0*(winsize0+1)^2 || DgDxImgMaskIndCount>0.0*(winsize0+1)^2
                  
                %%%%% Find connected region to deal with possible continuities %%%%%
                % tempg_BW2 = logical(tempg);
                tempg_BW2 = bwselect(logical(tempg), floor((winsize+1)/2), floor((winsize+1)/2), 8 );
                
                [rowtemp,~] = find(tempg_BW2==0);
                if isempty(rowtemp)
                    tempg_BW2 = tempfImgMask;
                    tempg_BW2 = bwselect(tempg_BW2, floor((winsize+1)/2), floor((winsize+1)/2), 8 );
                end
                tempg = tempg .* double(tempg_BW2);

                % tempf = ImgRef([x(1):1:x(3)],[y(1):1:y(3)]) .* tempfImgMask;
                % tempf = tempf .* double(tempg_BW2);
                % 
                % DfDx = Df.DfDx((x(1)-DfCropWidth):1:(x(3)-DfCropWidth), (y(1)-DfCropWidth):1:(y(3)-DfCropWidth));
                % DfDy = Df.DfDy((x(1)-DfCropWidth):1:(x(3)-DfCropWidth), (y(1)-DfCropWidth):1:(y(3)-DfCropWidth));
                % 
                % DfDx = DfDx .* double(tempg_BW2);
                % DfDy = DfDy .* double(tempg_BW2);
                % 
                % %%%%%%%%%%%%%%%%%%%% not sure why I updated H here JY!!!!
                % %%%%%%%%%%%%%%%%%%%% need to test here again TBD!!!!!!
                % dX = XX(:)-x0; dY = YY(:)-y0;
                % 
                % DfDx_column = DfDx(:);
                % DfDy_column = DfDy(:);
                % 
                % dFdWdP = [DfDx_column, DfDx_column.*dX, DfDx_column.*dY, (DfDx_column.*dX.^2)/2, DfDx_column.*dX.*dY,...
                %     (DfDx_column.*dY.^2)/2, DfDy_column, DfDy_column.*dX, DfDy_column.*dY, (DfDy_column.*dX.^2)/2, ...
                %     DfDy_column.*dX.*dY, (DfDy_column.*dY.^2)/2 ];
                % 
                % H = dFdWdP'*dFdWdP;
                % 
                % %%%%%%%%%% !!!Mask: START %%%%%%%%%%%%
                % meanf = mean(tempf(abs(tempf)>1e-10));
                % bottomf = sqrt((length(tempf(abs(tempf)>1e-10))-1)*var(tempf(abs(tempf)>1e-10)));
                % %%%%%%%%%% !!!Mask: END %%%%%%%%%%%%

                % figure(1); clf;
                % subplot(3,2,1), surf(DfDx,'edgecolor','none'); title('DfDx');view(2); axis equal; axis tight;
                % subplot(3,2,2), surf(DfDy,'edgecolor','none'); title('DfDy');view(2); axis equal; axis tight;
                % subplot(3,2,3), surf(tempf0,'edgecolor','none'); title('tempf'); view(2); axis equal; axis tight;caxis([-1,1]);
                % subplot(3,2,4), surf(tempg0,'edgecolor','none'); title('tempg'); view(2); axis equal; axis tight; caxis([-1,1]);
                % subplot(3,2,5), imshow(flipud( tempfImgMask) ); title('im f mask');
                % subplot(3,2,6), imshow(flipud( tempg_BW2) ); title('im g mask');
                % 
                % pause;

            end
            %%%%%%%%%% !!!Mask: END %%%%%%%%%%%

    
            % figure(1); clf;
            % subplot(3,2,1), surf(DfDx,'edgecolor','none'); title('DfDx');view(2); axis equal; axis tight;
            % subplot(3,2,2), surf(DfDy,'edgecolor','none'); title('DfDy');view(2); axis equal; axis tight;
            % subplot(3,2,3), surf(tempf,'edgecolor','none'); title('tempf'); view(2); axis equal; axis tight;caxis([-1,1]);
            % subplot(3,2,4), surf(tempg,'edgecolor','none'); title('tempg'); view(2); axis equal; axis tight; caxis([-1,1]);
            % subplot(3,2,5), imshow(flipud( tempfImgMask) ); title('im f mask');
            % 
            % 
            % pause;
                
                
            % ====== Old version codes ======
            % tempg = zeros(size(tempf,1)*size(tempf,2),1);
            % [tempCoordy, tempCoordx] = meshgrid(1:winsize+1,1:winsize+1);
            % tempCoordx = tempCoordx(:); tempCoordy = tempCoordy(:);
            %
            % for tempij = 1:size(tempCoordx,1)
            %     tempg(tempij)= ...
            %         fungInterpolation_g(u22(tempCoordx(tempij),tempCoordy(tempij)), v22(tempCoordx(tempij),tempCoordy(tempij)), ...
            %         g(floor(u22(tempCoordx(tempij),tempCoordy(tempij)))-1:floor(u22(tempCoordx(tempij),tempCoordy(tempij)))+2, ...
            %         floor(v22(tempCoordx(tempij),tempCoordy(tempij)))-1:floor(v22(tempCoordx(tempij),tempCoordy(tempij)))+2));
            % end
            %
            % tempg = reshape(tempg, winsize+1, winsize+1);
            % ===============================
            
            % A = [1+P(1) P(2) 0; P(3) 1+P(4) 0; P(5) P(6) 1];
            % tform = affine2d((A));
            %
            % tempg2 = g((x(1)-winsize/2):(x(3)+winsize/2), (y(1)-winsize/2):(y(3)+winsize/2));
            % tempg3 = imwarp(tempg2,tform,'cubic');
            %
            % figure; imshow(tempf,[]);
            % figure; imshow(tempg2,[]);
            % figure; imshow(tempg3,[]);
            %
            % [M,N] = size(tempg3)
            % tempg = tempg3(ceil((M+1)/2)-winsize/2:ceil((M+1)/2)+winsize/2, ceil((N+1)/2)-winsize/2:ceil((N+1)/2)+winsize/2);
            % figure; imshow(tempg,[]);
            
            %%%%%%%%%% !!!Mask: START %%%%%%%%%%%%
            %%%%% Find connected region to deal with possible continuities %%%%%
            meang = mean(tempg(abs(tempg)>1e-10));
            bottomg = sqrt((length(tempg(abs(tempg)>1e-10))-1)*var(tempg(abs(tempg)>1e-10)));
            %%%%%%%%%% !!!Mask: END %%%%%%%%%%%%
            
            % ============ For Levenberg-Marquardt method ============
            switch ICGNmethod
                case 'LevenbergMarquardt'

%%% TBD 
                    % % Compute functinoal error
                    % KappaOld = KappaNew;
                    % Kappatemp = (tempf-meanf)/bottomf - (tempg-meang)/bottomg;
                    % Kappatemp = Kappatemp.*Kappatemp;
                    % KappaNew = sum(Kappatemp(:));
                    % 
                    % if KappaNew < 1.02*KappaOld
                    %     delta = delta/10;
                    % else
                    %     delta = delta*10;
                    %     % Perform P inverse
                    %     DeltaP = -DeltaP;
                    %     tempP1 =  (-DeltaP(1)-DeltaP(1)*DeltaP(4)+DeltaP(2)*DeltaP(3))/temp;
                    %     tempP2 =  -DeltaP(2)/temp;
                    %     tempP3 =  -DeltaP(3)/temp;
                    %     tempP4 =  (-DeltaP(4)-DeltaP(1)*DeltaP(4)+DeltaP(2)*DeltaP(3))/temp;
                    %     tempP5 =  (-DeltaP(5)-DeltaP(4)*DeltaP(5)+DeltaP(3)*DeltaP(6))/temp;
                    %     tempP6 =  (-DeltaP(6)-DeltaP(1)*DeltaP(6)+DeltaP(2)*DeltaP(5))/temp;
                    % 
                    %     tempMatrix = [1+P(1) P(3) P(5); P(2) 1+P(4) P(6); 0 0 1]*...
                    %         [1+tempP1 tempP3 tempP5; tempP2 1+tempP4 tempP6; 0 0 1];
                    % 
                    %     P1 = tempMatrix(1,1)-1;
                    %     P2 = tempMatrix(2,1);
                    %     P3 = tempMatrix(1,2);
                    %     P4 = tempMatrix(2,2)-1;
                    %     P5 = tempMatrix(1,3);
                    %     P6 = tempMatrix(2,3);
                    %     P = [P1 P2 P3 P4 P5 P6]';
                    % end
                    % 
                    % % Find region for g
                    % % [tempCoordy, tempCoordx] = meshgrid(y(1):y(3),x(1):x(3));
                    % %Repeated! tempCoordx = XX - x0*ones(winsize+1,winsize+1);
                    % %Repeated! tempCoordy = YY - y0*ones(winsize+1,winsize+1);
                    % u22 = (1+P(1))*tempCoordxMat + P(3)*tempCoordyMat + (x0+P(5))*ones(winsize+1,winsize+1);
                    % v22 = P(2)*tempCoordxMat + (1+P(4))*tempCoordyMat + (y0+P(6))*ones(winsize+1,winsize+1);
                    % 
                    % tempg = ImgDef.eval(u22,v22);
                    % 
                    % %%%%%%%%%% !!!Mask: START %%%%%%%%%%%%
                    % %%%%% Find connected region to deal with possible continuities %%%%%
                    % tempg_BW2 = bwselect(logical(tempg), floor((winsize+1)/2), floor((winsize+1)/2), 4 );
                    % tempg = tempg .* double(tempg_BW2);
                    % %%%%%%%%%% !!!Mask: END %%%%%%%%%%%%
                    % 
                    % % ====== Old version codes ======
                    % % tempg = zeros(size(tempf,1)*size(tempf,2),1);
                    % %
                    % % [tempCoordy, tempCoordx] = meshgrid(1:winsize+1,1:winsize+1);
                    % % tempCoordx = tempCoordx(:); tempCoordy = tempCoordy(:);
                    % %
                    % % parfor tempij = 1:size(tempCoordx,1)
                    % %     tempg(tempij)= ...
                    % %         fungInterpolation_g(u22(tempCoordx(tempij),tempCoordy(tempij)), v22(tempCoordx(tempij),tempCoordy(tempij)), ...
                    % %         g(floor(u22(tempCoordx(tempij),tempCoordy(tempij)))-1:floor(u22(tempCoordx(tempij),tempCoordy(tempij)))+2, ...
                    % %         floor(v22(tempCoordx(tempij),tempCoordy(tempij)))-1:floor(v22(tempCoordx(tempij),tempCoordy(tempij)))+2));
                    % % end
                    % %
                    % % tempg = reshape(tempg, winsize+1, winsize+1);
                    % % ==================================
                    % 
                    % %%%%%%%%%% !!!Mask: START %%%%%%%%%%%%
                    % %%%%% Find connected region to deal with possible continuities %%%%%
                    % meang = mean(tempg(abs(tempg)>1e-10));
                    % bottomg = sqrt((length(tempg(abs(tempg)>1e-10))-1)*var(tempg(abs(tempg)>1e-10)));
                    %%%%%%%%%% !!!Mask: END %%%%%%%%%%%%
                    
                otherwise
            end
            
            % % ============ End of Levenberg-Marquardt method ============
            
            % ====== Assemble b vector old version ======
            % b = zeros(6,1);
            %
            % %[tempCoordy, tempCoordx] = meshgrid(y(1):y(3),x(1):x(3));
            % %tempCoordx = tempCoordx(:); tempCoordy = tempCoordy(:);
            %
            % for tempij = 1:size(tempCoordx,1)
            %     b = b + bottomf*([DfDx(tempCoordx(tempij)-DfCropWidth,tempCoordy(tempij)-DfCropWidth) DfDy(tempCoordx(tempij)-DfCropWidth,tempCoordy(tempij)-DfCropWidth)]*...
            %             [tempCoordx(tempij)-x0 0 tempCoordy(tempij)-y0 0 1 0; 0 tempCoordx(tempij)-x0 0 tempCoordy(tempij)-y0 0 1])'* ...
            %             ((tempf(tempCoordx(tempij)+1-x(1), tempCoordy(tempij)+1-y(1))-meanf)/bottomf - ...
            %             (tempg(tempCoordx(tempij)+1-x(1), tempCoordy(tempij)+1-y(1))-meang)/bottomg);
            % end
            % ====== Assemble b vector fast version ======

            b = dFdWdP'*(tempf(:)-meanf-bottomf/bottomg*(tempg(:)-meang));  
 


            normOfWOld = normOfWNew;
            normOfWNew = norm(b(:)); normOfWNewAbs = normOfWNew;
            
            if stepwithinwhile ==1
                normOfWNewInit = normOfWNew;
            end
            if normOfWNewInit > tol
                normOfWNew = normOfWNew/normOfWNewInit;
            else
                normOfWNew = 0;
            end
            
            StopCrit=@(dP,Zeta) sqrt(sum((dP'.*[1,Zeta,Zeta,0.5* Zeta.^2,Zeta.^2,0.5*Zeta.^2,1,Zeta,Zeta,0.5*Zeta.^2,Zeta.^2 ,0.5*Zeta.^2]).^2));
            StopVal = 1;


            if (normOfWNew<tol) || (normOfWNewAbs<tol)
                break
            else
                 
                DeltaP = -(H + delta*diag(diag(H))) \ b;
                DeltaPMatrix = funSFP2(DeltaP);
 
                if (det(DeltaPMatrix) ~= 0)
 
                    PMatrix = funSFP2(P);

                    tempMatrix = PMatrix/(DeltaPMatrix);

                    P1 = tempMatrix(4,6);
                    P2 = tempMatrix(4,4)-1;
                    P3 = tempMatrix(4,5);
                    P4 = tempMatrix(4,1)*2;
                    P5 = tempMatrix(4,2);
                    P6 = tempMatrix(4,3)*2;
                    P7 = tempMatrix(5,6);
                    P8 = tempMatrix(5,4);
                    P9 = tempMatrix(5,5)-1;
                    P10 = tempMatrix(5,1)*2;
                    P11 = tempMatrix(5,2);
                    P12 = tempMatrix(5,3)*2;

                    P = [P1 P2 P3 P4 P5 P6 P7 P8 P9 P10 P11 P12]';
 

                    StopVal=StopCrit(DeltaP,(winsize)/2);
 
                else
                    disp(['Det(DeltaP)==0!'])
                    break
                end
                
            end
        end
    end % end of while
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    if (normOfWNew<tol) || (normOfWNewAbs<tol)
        % elementsLocalMethodConvergeOrNot = 1;
    else
        stepwithinwhile = maxIterNum+1; %101
    end
    
    if (isnan(normOfWNew)==1)
        stepwithinwhile = maxIterNum+1; %101
    end
    if sum(abs(tempf(:))) < 1e-6
        stepwithinwhile = maxIterNum+3; %103
    end
        
    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
else % if norm(diag(H)) > abs(eps)
    
    H = zeros(12,12);
    stepwithinwhile = maxIterNum+2; %102
    
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% P = [u, u_x, u_y, u_xx, u_xy, u_yy, v, v_x, v_y, v_xx, v_xy, v_yy]';

% Save U
U(1) = P(1); U(2) = P(7);

% Save F
F(1) = P(2); F(2) = P(8); F(3) = P(3); F(4) = P(9); % U,x V,y U,x V,y

% Save HGlobal matrix
HGlobal = [H(1:12) H(14:24) H(27:36) H(40:48) H(53:60) H(66:72) ...
        H(79:84) H(92:96) H(105:108) H(118:120) H(131:132) H(144)];
 

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PMatrix = funSFP2(P)

    A1 = 2*P(2) + P(2)^2 + P(1)*P(4);
    A2 = 2*P(1)*P(5) + 2*(1+P(2))*P(3);
    A3 = P(3)^2 + P(1)*P(6);
    A4 = 2*P(1)*(1+P(2));
    A5 = 2*P(1)*P(3);
    A6 = P(1)^2;
    A7 = 0.5*( P(7)*P(4) + 2*(1+P(2))*P(8) + P(1)*P(10) );
    A8 = P(3)*P(8) + P(2)*P(9) + P(7)*P(5) + P(1)*P(11) + P(9) + P(2);
    A9 = 0.5*( P(7)*P(6) + 2*(1+P(9))*P(3) + P(1)*P(12) );
    A10 = P(7) + P(7)*P(2) + P(1)*P(8);
    A11 = P(1) + P(7)*P(3) + P(1)*P(9);
    A12 = P(1)*P(7);
    A13 = P(8)^2 + P(7)*P(10);
    A14 = 2*P(7)*P(11) + 2*P(8)*(1+P(9));
    A15 = 2*P(9) + P(9)^2 + P(7)*P(12);
    A16 = 2*P(7)*P(8);
    A17 = 2*P(7)*(1+P(9));
    A18 = P(7)^2;
    
    
    PMatrix = [ 1+A1, A2, A3, A4, A5, A6;
        A7, 1+A8, A9, A10, A11, A12;
        A13, A14, 1+A15, A16, A17, A18;
        0.5*P(4), P(5), 0.5*P(6), 1+P(2), P(3), P(1);
        0.5*P(10), P(11), 0.5*P(12), P(8), 1+P(9), P(7);
        0, 0, 0, 0, 0, 1];

end
