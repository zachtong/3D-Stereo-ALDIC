function Plotdisp(U_2D,U_3D,x,y,sizeOfImg,CurrentImg,DICpara)
%FUNCTION Plotdisp(U,x,y,sizeOfImg,CurrentImg,DICpara)
% To plot DIC solved displacement components
% ----------------------------------------------
%
%   INPUT: U                 Displacement vector: U = [Ux_node1, Uy_node1, Ux_node2, Uy_node2, ... , Ux_nodeN, Uy_nodeN]';
%          x,y               FE mesh x- and y-coordinates
%          sizeOfImg         DIC deformed image size
%          CurrentImg        Current deformed image
%          DICpara           DIC paramters
%
%   OUTPUT: Plots of x-displacement field and y-displacement field.
%
%   TODO: users could change caxis range based on their own choices.
%
% ----------------------------------------------
% Reference
% [1] RegularizeNd. Matlab File Exchange open source.
% https://www.mathworks.com/matlabcentral/fileexchange/61436-regularizend
% [2] Gridfit. Matlab File Exchange open source.
% https://www.mathworks.com/matlabcentral/fileexchange/8998-surface-fitting-using-gridfit
% ----------------------------------------------
% Author: Jin Yang.
% Contact and support: jyang526@wisc.edu -or- aldicdvc@gmail.com
% Last time updated: 11/2020.
% ==============================================


%% Initialization
warning off; load('./plotFiles/colormap_RdYlBu.mat','cMap');

OrigDICImgTransparency = DICpara.OrigDICImgTransparency; % Original raw DIC image transparency
Image2PlotResults = DICpara.Image2PlotResults; % Choose image to plot over (first only, second and next images)


%% For task1, I commented it
% if DICpara.DICIncOrNot == 0
%     temp(:,1) = U_2D.U(1:2:end,:);
%     temp(:,2) = U_2D.U(2:2:end,:);
%     U_2D = temp;
% end
%
% M = size(x,1); N = size(x,2);
% u = U_3D{1,1}; v = U_3D{1,2}; w = U_3D{1,3};
% u = reshape(u,M,N); v = reshape(v,M,N); w = reshape(w,M,N);
%     if M < 4, x2 = x(:,1)'; else x2 = linspace(x(1,1),x(end,1),2*(length(x(:,1))-1)+1); x2=x2(:)'; end % x2 = interp(x(:,1)',4);
%     if N < 4, y2 = y(1,:); else y2 = linspace(y(1,1),y(1,end),2*(length(y(1,:))-1)+1); y2=y2(:)'; end % y2 = interp(y(1,:),4);
%
%
%     %% Compute displacement components to manipulate the reference image
%     disp_u = gridfit(reshape(x,M*N,1),reshape(y,M*N,1),reshape(u,M*N,1),x2,y2);
%     disp_v = gridfit(reshape(x,M*N,1),reshape(y,M*N,1),reshape(v,M*N,1),x2,y2);
%     disp_w = gridfit(reshape(x,M*N,1),reshape(y,M*N,1),reshape(w,M*N,1),x2,y2);
%     disp_x = gridfit(reshape(x,M*N,1),reshape(y,M*N,1),reshape(U_2D(:,1),M*N,1),x2,y2);
%     disp_y = gridfit(reshape(x,M*N,1),reshape(y,M*N,1),reshape(U_2D(:,2),M*N,1),x2,y2);
%
%
% % Please don't delete this line, to deal with the image and physical world coordinates
% [x2,y2]=ndgrid(x2,y2); x2=x2'; y2=y2';
% 


%% task 1 ------------------------
M = size(x,1); N = size(x,2);
u = U_3D{1,1}; v = U_3D{1,2}; w = U_3D{1,3};
x2 = x'; y2 = y';
disp_x = reshape(U_2D(:,1),N,M);
disp_y = reshape(U_2D(:,2),N,M);

disp_u = reshape(u,M,N); disp_u = disp_u';
disp_v = reshape(v,M,N); disp_v = disp_v';
disp_w = reshape(w,M,N); disp_w = disp_w';




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== 1) dispx u ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig1=figure; ax1=axes;
try h1 = imshow(  flipud(imread(CurrentImg)),'InitialMagnification','fit'); % why?????
catch h1 = surf(  flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
end

axis equal; axis tight; box on; set(gca,'fontSize',18); view(2); set(gca,'ydir','normal');
hold on; ax2=axes; h2 = surf( (x2+Image2PlotResults*disp_x), ...
    (sizeOfImg(2)+1)-((y2+Image2PlotResults*disp_y)), disp_u, 'EdgeColor','none','LineStyle','none');

% Debug
% xtemp = (x2+Image2PlotResults*disp_x);
% ytemp = (sizeOfImg(2)+1)-((y2+Image2PlotResults*disp_y));
% hold on; ax2=axes; h2 = surf( xtemp(3:end-3,7:end-7), ytemp(3:end-3,7:end-7) , disp_u(3:end-3,7:end-7), 'EdgeColor','none','LineStyle','none');

set(gca,'fontSize',18); view(2); box on; % set(gca,'ydir','normal');
alpha(h2,OrigDICImgTransparency); axis equal;  axis tight; colormap(cMap);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% TODO: manually modify colormap and caxis %%%%%%
% colormap(jet); caxis([-1.7,.2]); % caxis([-0.025,0.025]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

linkaxes([ax1,ax2]);  % Link axes together
ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = []; % Hide the top axes
colormap(ax1,'gray'); % Give each one its own colormap
set([ax1,ax2],'Position',[.17 .11 .685 .815]);
ax1.Visible = 'on'; ax1.TickLabelInterpreter = 'latex';
%%%%% convert pixel unit to the physical world unit %%%%%
xticklabels(ax1, num2cell(round(ax1.XTick*100)/100, length(ax1.XTick) )' );
yticklabels(ax1, num2cell(round(ax1.YTick*100)/100, length(ax1.YTick) )' );
cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'latex';

% xlabel( '$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
% title('$x-$displacement $u$','FontWeight','Normal','Interpreter','latex'); set(gcf,'color','w');



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== 2) dispx v ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig1=figure; ax1=axes;
try h1=imshow( flipud(imread(CurrentImg)),'InitialMagnification','fit');
catch h1=surf(  flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
end

axis equal; axis tight; box on; set(gca,'fontSize',18); view(2); set(gca,'ydir','normal');
hold on; ax2=axes; h2 = surf( (x2+Image2PlotResults*disp_x)...
    , (sizeOfImg(2)+1)-((y2+Image2PlotResults*disp_y)), disp_v, 'EdgeColor','none','LineStyle','none');
set(gca,'fontSize',18); view(2); box on; %set(gca,'ydir','normal');
alpha(h2,OrigDICImgTransparency); axis equal;  axis tight; colormap(cMap);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% TODO: manually modify colormap and caxis %%%%%%
% colormap(jet); caxis([5,12]); % caxis([-0.025,0.025]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

linkaxes([ax1,ax2]);  % Link axes together
ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = []; % Hide the top axes
colormap(ax1,'gray'); % Give each one its own colormap
set([ax1,ax2],'Position',[.17 .11 .685 .815]);
ax1.Visible = 'on'; ax1.TickLabelInterpreter = 'latex';
%%%%% convert pixel unit to the physical world unit %%%%%
xticklabels(ax1, num2cell(round(ax1.XTick*100)/100, length(ax1.XTick) )' );
yticklabels(ax1, num2cell(round(ax1.YTick*100)/100, length(ax1.YTick) )' );
cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'latex';

% xlabel( '$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
% title('$y-$displacement $v$','FontWeight','Normal','Interpreter','latex'); set(gcf,'color','w');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== 3) dispx w ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig1=figure; ax1=axes;
try h1=imshow( flipud(imread(CurrentImg)),'InitialMagnification','fit');
catch h1=surf(  flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
end

axis equal; axis tight; box on; set(gca,'fontSize',18); view(2); set(gca,'ydir','normal');
hold on; ax2=axes; h2 = surf( (x2+Image2PlotResults*disp_x)...
    , (sizeOfImg(2)+1)-((y2+Image2PlotResults*disp_y)), disp_w, 'EdgeColor','none','LineStyle','none');
set(gca,'fontSize',18); view(2); box on; %set(gca,'ydir','normal');
alpha(h2,OrigDICImgTransparency); axis equal;  axis tight; colormap(cMap);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% TODO: manually modify colormap and caxis %%%%%%
% colormap(jet); caxis([5,12]); % caxis([-0.025,0.025]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

linkaxes([ax1,ax2]);  % Link axes together
ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = []; % Hide the top axes
colormap(ax1,'gray'); % Give each one its own colormap
set([ax1,ax2],'Position',[.17 .11 .685 .815]);
ax1.Visible = 'on'; ax1.TickLabelInterpreter = 'latex';
%%%%% convert pixel unit to the physical world unit %%%%%
xticklabels(ax1, num2cell(round(ax1.XTick*100)/100, length(ax1.XTick) )' );
yticklabels(ax1, num2cell(round(ax1.YTick*100)/100, length(ax1.YTick) )' );
cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'latex';

% xlabel( '$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
% title('$y-$displacement $v$','FontWeight','Normal','Interpreter','latex'); set(gcf,'color','w');




