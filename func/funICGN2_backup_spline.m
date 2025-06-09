function [U,F,stepwithinwhile,HGlobal] = funICGN2_backup_spline(U0,x0,y0,Df,ImgRef,ImgDef,winsize,tol,ICGNmethod,InterpCoef)
%FUNCTION [U,F,stepwithinwhile,HGlobal] = funICGN(U0,x0,y0,Df,ImgRef,ImgDef,winsize,tol,method)
% Refer to the SubCorr.m in ADIC3D.
% There are some image coordinate system transformation because ADIC3D and
% 3D-ALDIC have different image coordinate index orders. Codes for
% transformation are with "transformation" comments.
% ==============================================

%% Convert coordinate
[x0, y0] = swap(x0, y0); % transformation
ImgRef = ImgRef'; % transformation
ImgDef = ImgDef'; % transformation


%% Initialization
warning('off');
DfDxStartx = Df.DfAxis(1); DfDxStarty = Df.DfAxis(3);
imgSize = Df.imgSize;

%% Declare anonymous function
SFPVec2Mat=@(P) reshape([P(2)*2+P(1)*P(4)+P(2)^2+1,P(1)*P(10)*1/2+P(4)*P(7)*(1/2)+P(8)*(P(2)*2+2)*1/2,...
        P(7)*P(10)+P(8)^2,P(4)*1/2,P(10)*1/2,0,P(1)*P(5)*2+P(3)*(P(2)*2+2),...
        P(2)+P(9)+P(2)*P(9)+P(3)*P(8)+P(1)*P(11)+P(5)*P(7)+1,P(7)*P(11)*2.0+P(8)*(P(9)+1)*2,P(5),P(11)...
        0,P(1)*P(6)+P(3)^2,P(1)*P(12)*1/2+P(6)*P(7)*1/2+P(3)*(P(9)+1),P(9)*2+P(7)*P(12)+P(9)^2+1,...
        P(6)*1/2,P(12)*1/2,0,P(1)*(P(2)+1)*2,P(7)+P(1)*P(8)+P(2)*P(7),P(7)*P(8)*2,P(2)+1,...
        P(8),0,P(1)*P(3)*2,P(1)+P(1)*P(9)+P(3)*P(7),P(7)*(P(9)+1)*2,P(3),P(9)+1,0,P(1)^2,...
        P(1)*P(7),P(7)^2,P(1),P(7),1],[6,6]);

Mat2SFPVec=@(W) [W(34),W(22)-1,W(28),W(4).*2,W(10),W(16) .*2,W(35), W(23),W(29)-1,W(5).*2,W(11),W(17).*2];
StopCrit=@(dP,Zeta) sqrt(sum((dP'.*[1,Zeta,Zeta,0.5* Zeta.^2,Zeta.^2,0.5*Zeta.^2,1,Zeta,Zeta,0.5*Zeta.^2,Zeta.^2 ,0.5*Zeta.^2]).^2));
W=@(dX,dY,P) [P(1)+P(3).*dY+P(4).*dX.^2.*(1/2)+P(6).* dY.^2.*(1/2)+dX.*(P(2)+1)+P(5).*dX.*dY,...
        P(7)+P(8).*dX+ P(10).*dX.^2.*(1/2)+P(12).*dY.^2.*(1/2)+dY.*(P(9)+1)+ P(11).*dX.*dY];   

%% ---------------------------
% Find local subset region

x = [x0-winsize/2 ; x0+winsize/2 ; x0+winsize/2 ; x0-winsize/2];  % [coordinates(elements(j,:),1)];
y = [y0-winsize/2 ; y0+winsize/2 ; y0+winsize/2 ; y0-winsize/2];  % [coordinates(elements(j,:),2)];

% ---------------------------
% Initialization: Get P0
P = [U0(1) 0 0 0 0 0 U0(2) 0 0 0 0 0 ]';

% ---------------------------
% Find region for f
tempf = ImgRef([x(1):1:x(3)],[y(1):1:y(3)]);
meanf = mean(tempf(:));
bottomf = sqrt((length(tempf(:))-1)*var(tempf(:)));

DfDx = Df.DfDx((y(1)-DfDxStarty):1:(y(3)-DfDxStarty), (x(1)-DfDxStartx):1:(x(3)-DfDxStartx));
DfDy = Df.DfDy((y(1)-DfDxStarty):1:(y(3)-DfDxStarty), (x(1)-DfDxStartx):1:(x(3)-DfDxStartx));
DfDx = DfDx';% transformation
DfDy = DfDy';% transformation


%% Hessian matrix
[dX,dY]=meshgrid(-(winsize)/2:(winsize)/2,-(winsize)/2: (winsize)/2);
dX=dX(:); dY=dY(:); DfDx_column = DfDx(:); DfDy_column = DfDy(:);
dFdWdP = [DfDx_column,DfDx_column.*dX,DfDx_column.*dY, (DfDx_column.*dX.^2)/2,DfDx_column.*dX.*dY,...
    (DfDx_column.*dY.^2)/2,DfDy_column,DfDy_column.*dX ,DfDy_column.*dY,(DfDy_column.*dX.^2)/2,...
    DfDy_column.*dX.*dY,(DfDy_column.*dY.^2)/2];

H = dFdWdP'*dFdWdP;
H2=inv(H);

%% Iterate
% Initialize while loop

stepwithinwhile = 0; StopVal = 1;
while( (stepwithinwhile<=100) && StopVal>tol ) % .

    stepwithinwhile = stepwithinwhile+1;

    dXY = W(dX,dY,P) ; % the displacements of each pixel within window

    tempg = InterpCoef(x0.*ones(size(dXY,1),1)+dXY(:,2), y0.*ones(size(dXY,1),1)+dXY(:,1));
    % tempg = reshape(tempg,[winsize+1,winsize+1]);

    meang = mean(tempg(:));
    bottomg = sqrt((length(tempg(:))-1)*var(tempg(:)));

    J=dFdWdP'*(tempf(:)-meanf-bottomf/bottomg*(tempg(:)-meang)); % Jacobian in ICGN
    dP=-H2*J; 
    
    % Update P
    P=Mat2SFPVec(SFPVec2Mat(P)/SFPVec2Mat(dP))'; % update

    StopVal=StopCrit(dP,(winsize)/2);


end % end of while
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save U
U(1) = P(1); U(2) = P(7);

% Save F and global_H
F(1) = P(2); F(2) = P(3); F(3) = P(8); F(4) = P(9); % U,x U,y V,x V,y
HGlobal = [H(1:12) H(14:24) H(27:36) H(40:48) H(53:60) H(66:72) H(79:84) H(92:96) H(105:108) H(118:120) H(131:132) H(144)]; %  Be changed to 12*12


end


function [a,b] = swap(a,b)
temp = a;
a =b;
b = temp;
end
