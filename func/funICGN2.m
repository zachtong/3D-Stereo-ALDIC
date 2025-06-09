function [U,F,stepwithinwhile,HGlobal,varargout] = funICGN2(...
    U0,x0,y0,Df,ImgRef,ImgDef,winsize,tol,ICGNmethod,varargin)
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


if nargin >= 10
    method = 'BicubicBspline';
    BicubicBspline = varargin{1};
else
    method = 'Bicubic';
end

%% ---------------------------
% Find local subset region
x = [x0-winsize/2 ; x0+winsize/2 ; x0+winsize/2 ; x0-winsize/2];  % [coordinates(elements(j,:),1)];
y = [y0-winsize/2 ; y0+winsize/2 ; y0+winsize/2 ; y0-winsize/2];  % [coordinates(elements(j,:),2)];

% ---------------------------
% Initialization: P = [u, u_x, u_y, u_xx, u_xy, u_yy, v, v_x, v_y, v_xx, v_xy, v_yy]';
P0 = [U0(1) zeros(1,5), U0(2), zeros(1,5)]';
P = P0;

% ---------------------------
% Find region for f
[XX,YY] = ndgrid([x(1):1:x(3)],[y(1):1:y(3)]);

tempf = ImgRef([x(1):1:x(3)],[y(1):1:y(3)]);
DfDx = Df.DfDx((x(1)-DfCropWidth):1:(x(3)-DfCropWidth), (y(1)-DfCropWidth):1:(y(3)-DfCropWidth));
DfDy = Df.DfDy((x(1)-DfCropWidth):1:(x(3)-DfCropWidth), (y(1)-DfCropWidth):1:(y(3)-DfCropWidth));


%% %%%%%%%% If there are >50% of the subset are painted with patterns %%%%%%%%%%%%
[DfDxImgMaskIndRow,~] = find(DfDx==0);
if length(DfDxImgMaskIndRow)<0.50*(winsize+1)^2

    if length(DfDxImgMaskIndRow)>0.1*(winsize+1)^2
        winsize = max(round(sqrt(2)*winsize));
        x = [x0-winsize/2 ; x0+winsize/2 ; x0+winsize/2 ; x0-winsize/2]; % Update x
        y = [y0-winsize/2 ; y0+winsize/2 ; y0+winsize/2 ; y0-winsize/2]; % Update y
        [XX,YY] = ndgrid([x(1):1:x(3)],[y(1):1:y(3)]);
        tempf = ImgRef([x(1):1:x(3)],[y(1):1:y(3)]);
        DfDx = Df.DfDx((x(1)-DfCropWidth):1:(x(3)-DfCropWidth), (y(1)-DfCropWidth):1:(y(3)-DfCropWidth));
        DfDy = Df.DfDy((x(1)-DfCropWidth):1:(x(3)-DfCropWidth), (y(1)-DfCropWidth):1:(y(3)-DfCropWidth));
    end

    dX = XX(:)-x0; dY = YY(:)-y0;

    DfDx_column = DfDx(:);
    DfDy_column = DfDy(:);

    dFdWdP = [DfDx_column, DfDx_column.*dX, DfDx_column.*dY, (DfDx_column.*dX.^2)/2, DfDx_column.*dX.*dY,...
        (DfDx_column.*dY.^2)/2, DfDy_column, DfDy_column.*dX, DfDy_column.*dY, (DfDy_column.*dX.^2)/2, ...
        DfDy_column.*dX.*dY, (DfDy_column.*dY.^2)/2 ];

    H = dFdWdP'*dFdWdP;

    meanf = mean(tempf(:));
    bottomf = sqrt((length(tempf(:))-1)*var(tempf(:)));

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

    while( (stepwithinwhile<=100) && (normOfWNew>tol) && (normOfWNewAbs>tol) )

        stepwithinwhile = stepwithinwhile+1;

        % Find region for g
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

            %tempg = ba_interp2(ImgDef, v22, u22, 'cubic');

            % BicubicBspline interpolation
            switch method
                case 'Bicubic'
                    tempg = ba_interp2(ImgDef, v22, u22, 'cubic');
                case 'BicubicBspline'
                    % BicubicBspline interpolation
                    [tempg,BicubicBspline] = BicubicBsplineInterp(BicubicBspline,ImgDef,v22,u22);
            end



            % subplot(1,2,1); surf(tempg,'EdgeColor','none'); view(2);
            % subplot(1,2,2); surf(tempf,'EdgeColor','none'); view(2);


            meang = mean(tempg(:));
            bottomg = sqrt((length(tempg(:))-1)*var(tempg(:)));

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

            if (normOfWNew<tol) || (normOfWNewAbs<tol) || (StopVal<tol)
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
        stepwithinwhile = 101;
    end

    if (isnan(normOfWNew)==1)
        stepwithinwhile = 101;
    end


    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else % if norm(diag(H)) > abs(eps)

    H = zeros(12,12);
    stepwithinwhile = 102;

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

if nargout >= 5
    varargout{1} = BicubicBspline;
end


end


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
