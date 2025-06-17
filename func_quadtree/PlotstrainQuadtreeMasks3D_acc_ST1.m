function [strain_exx,strain_exy,strain_eyy,strain_principal_max,strain_principal_min, ...
    strain_maxshear,strain_vonMises,dwdx,dwdy] = PlotstrainQuadtreeMasks3D_acc_ST1(U_2D_L,coefficients,voidIndex,FirstFEM,FirstImg,CurrentImg,DICpara)
%PLOTSTRAINQUADTREE: to plot DIC solved strain fields on a quadtree mesh 
% and overlaid with the original DIC images 
%   [strain_exx,strain_exy,strain_eyy,strain_principal_max,strain_principal_min, ...
%    strain_maxshear,strain_vonMises] = PlotstrainQuadtree(U,F,coordinatesFEMWorld,elementsFEM,CurrentImg,DICpara)
% ----------------------------------------------
%
%   INPUT: F                    DIC solved deformation gradient tensor
%          coordinatesFE        FE mesh coordinates
%          elementsFEM          FE mesh elements
%
%   OUTPUT: strain_exx              strain xx-compoent
%           strain_exy              strain xy-compoent
%           strain_eyy              strain yy-compoent
%           strain_principal_max    max principal strain on the xy-plane
%           strain_principal_min    min principal strain on the xy-plane
%           strain_maxshear         max shear strain on the xy-plane
%           strain_vonMises         equivalent von Mises strain
%
%   Plots:       
%       1) strain sxx
%       2) strain sxy
%       3) strain syy
%       4) max principal strain on the xy-plane 
%       5) min principal strain on the xy-plane
%       6) max shear strain on the xy-plane
%       7) equivalent von Mises strain
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
% Last time updated: 2020.12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Initialization
% warning off; load('./plotFiles/colormap_RdYlBu.mat','cMap');
% run('./plotFiles/Black_rainbow.m');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% convert pixel unit to the physical world unit %%%%%
try um2px = DICpara.um2px; 
catch um2px = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OrigDICImgTransparency = DICpara.OrigDICImgTransparency; % Original raw DIC image transparency
Image2PlotResults = DICpara.Image2PlotResults; % Choose image to plot over (first only, second and next images)

%disp_u = U{1}; disp_v = U{2}; disp_w = U{3};
%coordinatesFEMWorldDef = [coordinatesFEMWorld(:,1)+Image2PlotResults*disp_u, coordinatesFEMWorld(:,2)+Image2PlotResults*disp_v];

%%%%%%%%%%% JY!!!Mask START %%%%%%%%%%%%%%%
% if Image2PlotResults == 1
%     if ~isempty(CurrentImgMask)
%         for tempi = 1:size(coordinatesFEMWorldDef,1)
%             try
%             if CurrentImgMask( floor(coordinatesFEMWorldDef(tempi,1)/um2px), ...
%                                 (size(CurrentImgMask,2)+1-ceil(coordinatesFEMWorldDef(tempi,2)/um2px)) ) == 0 
%                 coordinatesFEMWorldDef(tempi,:) = [nan,nan];
%             end
%             catch
%                 coordinatesFEMWorldDef(tempi,:) = [nan,nan];
%             end
%         end
%     else
%         CurrentImgMask = imread(CurrentImg)';
%         for tempi = 1:size(coordinatesFEMWorldDef,1)
%             try
%                 if CurrentImgMask( floor(coordinatesFEMWorldDef(tempi,1)/um2px), ...
%                         (size(CurrentImgMask,2)+1-ceil(coordinatesFEMWorldDef(tempi,2)/um2px)) ) < 0
%                     coordinatesFEMWorldDef(tempi,:) = [nan,nan];
%                 end
%             catch
%                 coordinatesFEMWorldDef(tempi,:) = [nan,nan];
%             end
%         end
%         % [row1,~] = find(isnan(U(1:2:end))==1);
%         % F(4*row1-3) = nan; F(4*row1-2) = nan; F(4*row1-1) = nan; F(4*row1) = nan; 
%     end
% end
%%%%%%%%%%% JY!!!Mask END %%%%%%%%%%%%%%%


%% Compute strain components
% Method 1
% strain_exx = dudx; strain_eyy = dvdy; strain_ezz = dwdz;
% strain_exy = 0.5*(dvdx + dudy);
% strain_exz = 0.5*(dwdx + dudz);
% strain_eyz = 0.5*(dvdz + dwdy);

% Method 2: Strain_tensor = 0.5* (F'*F - I)
strain_exx = zeros(size(coefficients,1),1);
strain_exy = zeros(size(coefficients,1),1);
strain_eyy = zeros(size(coefficients,1),1);
strain_ezz = zeros(size(coefficients,1),1);
strain_exz = zeros(size(coefficients,1),1);
dwdx = zeros(size(coefficients,1),1);
dwdy = zeros(size(coefficients,1),1);

for i = 1:size(coefficients,1)
    u_x = coefficients{i,1}(1,1); u_y = coefficients{i,1}(2,1); u_z = coefficients{i,1}(3,1);
    v_x = coefficients{i,1}(1,2); v_y = coefficients{i,1}(2,2); v_z = coefficients{i,1}(3,2);
    w_x = coefficients{i,1}(1,3); w_y = coefficients{i,1}(2,3); w_z = coefficients{i,1}(3,3);
    F = [1+u_x, u_y, u_z; v_x, 1+v_y, v_z; w_x, w_y, 1+w_z];
    temp_Strain_tensor = 0.5*(F'*F-eye(3));
    strain_exx(i) = temp_Strain_tensor(1,1);
    strain_exy(i) = temp_Strain_tensor(1,2);
    strain_eyy(i) = temp_Strain_tensor(2,2);
    dwdx(i) = w_x;
    dwdy(i) = w_y;
end


%% Plot on deformed images or Not
if DICpara.Image2PlotResults == 1
    coordinatesFEMWorldDef = [FirstFEM.coordinatesFEM(:,1)+U_2D_L(:,1), FirstFEM.coordinatesFEM(:,2)+U_2D_L(:,2)];
    elementsFEM = FirstFEM.elementsFEM;
    Img = CurrentImg;
else
    coordinatesFEMWorldDef = [FirstFEM.coordinatesFEM(:,1), FirstFEM.coordinatesFEM(:,2)];
    elementsFEM = FirstFEM.elementsFEM;
    Img = FirstImg;
end

% strain_exx = strain_eyy;
% strain_exy = strain_eyy;

strain_maxshear = sqrt((0.5*(strain_exx-strain_eyy)).^2 + strain_exy.^2);
% Principal strain
strain_principal_max = 0.5*(strain_exx+strain_eyy) + strain_maxshear;
strain_principal_min = 0.5*(strain_exx+strain_eyy) - strain_maxshear;
% equivalent von Mises strain
strain_vonMises = sqrt(strain_principal_max.^2 + strain_principal_min.^2 - ...
             strain_principal_max.*strain_principal_min + 3*strain_maxshear.^2);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== 1) Strain exx ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig1=figure; ax1=axes; 
try h1=imshow( flipud(imread(Img)),'InitialMagnification','fit');
catch h1=surf( flipud( imread(Img) ),'EdgeColor','none','LineStyle','none');
end

axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2); set(gca,'ydir','normal');
hold on; ax2=axes; h2=show([],elementsFEM(:,1:4),coordinatesFEMWorldDef/um2px,strain_exx,'NoEdgeColor');
set(gca,'fontSize',18);set(gca,'ydir','reverse'); view(2); box on; axis equal;  axis tight;   
alpha(h2,OrigDICImgTransparency); colormap(turbo); caxis auto; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% TODO: manually modify colormap and caxis %%%%%%
% caxis([-0.02,0.02]);
% colormap(turbo); caxis([-1,1]) % D Sample 
% colormap(turbo);   caxis([-0.25,0.25]) % foam
 % colormap(turbo); caxis([-0.1 0.5]) % D shaped
%caxis([-0.04 0.025]) %bulge
%caxis([-0.05 0.1]); 
% colormap(turbo); caxis([-0.004,0]); % Sample 12
% colormap(black_rainbow_plus);    clim([-0.005 0.005]);
% colormap(black_rainbow); caxis([-0.004,0.004]);
% ax1.XTick = [100,200,300]; % Unit: px
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

linkaxes([ax1,ax2]);  % Link axes together
ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = []; % Hide the top axes
colormap(ax1,'gray'); % Give each one its own colormap
set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
ax1.Visible = 'on'; %ax1.TickLabelInterpreter = 'latex'; 
%%%%% convert pixel unit to the physical world unit %%%%%
xticklabels(ax1, num2cell(round(um2px*ax1.XTick*10)/10, length(ax1.XTick) )' );
yticklabels(ax1, num2cell(round(um2px*ax1.YTick*10)/10, length(ax1.YTick) )' );
% cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'latex';
cb2 = colorbar('Position',[.17+0.685+0.012 .11+.128 .03 .557 ]); %cb2.TickLabelInterpreter = 'latex';


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%====== 2) Strain exy ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig1=figure; ax1=axes; 
try h1=imshow( flipud(imread(Img) ),'InitialMagnification','fit');
catch h1=surf( flipud( imread(Img) ),'EdgeColor','none','LineStyle','none');
end

axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2); set(gca,'ydir','normal');
hold on; ax2=axes; h2=show([],elementsFEM(:,1:4),coordinatesFEMWorldDef/um2px,strain_exy,'NoEdgeColor');
set(gca,'fontSize',18);set(gca,'ydir','reverse'); view(2); box on; axis equal;  axis tight;   
alpha(h2,OrigDICImgTransparency); colormap(turbo); caxis auto;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% TODO: manually modify colormap and caxis %%%%%%
% caxis([-0.02,0.02]);
% colormap(turbo); caxis([-1,1]) % D Sample 
% colormap(turbo);  caxis([-0.25,0.25]) % foam
% colormap(turbo); caxis([-0.1 0.1]) % D shaped
%caxis([-0.017 0.017]) % bulge
%caxis([-0.06 0.07]) 
% colormap(turbo); caxis([-0.008,0.008]); % Sample 12 
 %colormap(black_rainbow_plus);    caxis([-0.15 0.15]);
%   colormap(black_rainbow); caxis([-0.0018,0.0018]);
% ax1.XTick = [100,200,300]; % Unit: px
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

linkaxes([ax1,ax2]);  % Link axes together
ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = []; % Hide the top axes
colormap(ax1,'gray'); % Give each one its own colormap
set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
ax1.Visible = 'on';% ax1.TickLabelInterpreter = 'latex'; 
%%%%% convert pixel unit to the physical world unit %%%%%
xticklabels(ax1, num2cell(round(um2px*ax1.XTick*10)/10, length(ax1.XTick) )' );
yticklabels(ax1, num2cell(round(um2px*ax1.YTick*10)/10, length(ax1.YTick) )' );
% cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'latex';
cb2 = colorbar('Position',[.17+0.685+0.012 .11+.128 .03 .557 ]); %cb2.TickLabelInterpreter = 'latex';



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== 3) Strain eyy ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig1=figure; ax1=axes; 
try h1=imshow( flipud(imread(Img)),'InitialMagnification','fit');
catch h1=surf( flipud( imread(Img) ),'EdgeColor','none','LineStyle','none');
end

axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2); set(gca,'ydir','normal');
hold on; ax2=axes; h2=show([],elementsFEM(:,1:4),coordinatesFEMWorldDef/um2px,strain_eyy,'NoEdgeColor');
set(gca,'fontSize',18); set(gca,'ydir','reverse'); view(2); box on; axis equal;  axis tight;   
alpha(h2,OrigDICImgTransparency); colormap(turbo); caxis auto;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% TODO: manually modify colormap and caxis %%%%%%
% caxis([-0.6,0.]);
% colormap(turbo); caxis([-0.15,0]) % D Sample 
% colormap(turbo);  caxis([-0.25,0.25])% foam
% caxis([-0.05 0.05]) 
% colormap(turbo); caxis([-0.2 -0.05]) % D shaped cb1
% colormap(turbo); caxis([-0.1  0.01]) % D shaped cb2
%caxis([-0.02 0.025]) % bulge
% colormap(turbo); caxis([-0.002,0.017]); % Sample 12 
%colormap(black_rainbow_plus);  caxis([-0.15 0.15]);
%  colormap(black_rainbow); caxis([-0.0021,0.0021]);
% ax1.XTick = [100,200,300]; % Unit: px
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

linkaxes([ax1,ax2]);  % Link axes together
ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = []; % Hide the top axes
colormap(ax1,'gray'); % Give each one its own colormap
set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
ax1.Visible = 'on';% ax1.TickLabelInterpreter = 'latex'; 
%%%%% convert pixel unit to the physical world unit %%%%%
xticklabels(ax1, num2cell(round(um2px*ax1.XTick*10)/10, length(ax1.XTick) )' );
yticklabels(ax1, num2cell(round(um2px*ax1.YTick*10)/10, length(ax1.YTick) )' );
% cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'latex';
cb2 = colorbar('Position',[.17+0.685+0.012 .11+.128 .03 .557 ]); %cb2.TickLabelInterpreter = 'latex';


% %% Zach's task
% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % ====== 3) Strain eyy_2   ======
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fig1=figure; ax1=axes; 
% try h1=imshow( flipud(imread(Img)),'InitialMagnification','fit');
% catch h1=surf( flipud( imread(Img) ),'EdgeColor','none','LineStyle','none');
% end
% 
% axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2); set(gca,'ydir','normal');
% hold on; ax2=axes; h2=show([],elementsFEM(:,1:4),coordinatesFEMWorldDef/um2px,strain_eyy,'NoEdgeColor');
% set(gca,'fontSize',18); set(gca,'ydir','reverse'); view(2); box on; axis equal;  axis tight;   
% alpha(h2,OrigDICImgTransparency); colormap(turbo); caxis auto;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%% TODO: manually modify colormap and caxis %%%%%%
% % caxis([-0.6,0.]);
% % colormap(turbo); caxis([-0.15,0]) % D Sample 
% % colormap(turbo);  caxis([-0.25,0.25])% foam
% 
% %colormap(turbo); caxis([-0.3 -0.1]) %  cb1
% colormap(turbo); caxis([-0.1  0.01])  % cb2
% 
% 
% 
% % colormap(turbo); caxis([-0.002,0.017]); % Sample 12 
% %colormap(black_rainbow_plus);  caxis([-0.15 0.15]);
% %  colormap(black_rainbow); caxis([-0.0021,0.0021]);
% % ax1.XTick = [100,200,300]; % Unit: px
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% linkaxes([ax1,ax2]);  % Link axes together
% ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = []; % Hide the top axes
% colormap(ax1,'gray'); % Give each one its own colormap
% set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
% ax1.Visible = 'on';% ax1.TickLabelInterpreter = 'latex'; 
% %%%%% convert pixel unit to the physical world unit %%%%%
% xticklabels(ax1, num2cell(round(um2px*ax1.XTick*10)/10, length(ax1.XTick) )' );
% yticklabels(ax1, num2cell(round(um2px*ax1.YTick*10)/10, length(ax1.YTick) )' );
% % cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'latex';
% cb2 = colorbar('Position',[.17+0.685+0.012 .11+.128 .03 .557 ]); %cb2.TickLabelInterpreter = 'latex';

  
% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % ====== 4) Strain e_principal_max ======
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fig1=figure; ax1=axes; 
% try h1=imshow( flipud(imread(Img)),'InitialMagnification','fit');
% catch h1=surf( flipud( imread(Img) ),'EdgeColor','none','LineStyle','none');
% end
% 
% axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2); set(gca,'ydir','normal');
% hold on; ax2=axes; h2=show([],elementsFEM(:,1:4),coordinatesFEMWorldDef/um2px,strain_principal_max,'NoEdgeColor');
% set(gca,'fontSize',18); set(gca,'ydir','reverse'); view(2); box on; axis equal;  axis tight;   
% alpha(h2,OrigDICImgTransparency); colormap(turbo);; caxis auto;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%% TODO: manually modify colormap and caxis %%%%%%
% caxis([0,0.45]);
% % colormap(turbo); caxis([-0.15,0]) % D Sample 
% % colormap(turbo);  caxis([-0.25,0.25])% foam
% % colormap(turbo); caxis([-0.002,0.017]); % Sample 12 
% %colormap(black_rainbow_plus);  caxis([-0.15 0.15]);
% %  colormap(black_rainbow); caxis([-0.0021,0.0021]);
% % ax1.XTick = [100,200,300]; % Unit: px
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% linkaxes([ax1,ax2]);  % Link axes together
% ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = []; % Hide the top axes
% colormap(ax1,'gray'); % Give each one its own colormap
% set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
% ax1.Visible = 'on';% ax1.TickLabelInterpreter = 'latex'; 
% %%%%% convert pixel unit to the physical world unit %%%%%
% xticklabels(ax1, num2cell(round(um2px*ax1.XTick*10)/10, length(ax1.XTick) )' );
% yticklabels(ax1, num2cell(round(um2px*ax1.YTick*10)/10, length(ax1.YTick) )' );
% % cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'latex';
% cb2 = colorbar('Position',[.17+0.685+0.012 .11+.128 .03 .557 ]); %cb2.TickLabelInterpreter = 'latex';
% % 




%   
% % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % ====== 5) Strain e_principal_min ======
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fig1=figure; ax1=axes; 
% try h1=imshow( flipud(imread(Img)),'InitialMagnification','fit');
% catch h1=surf(  flipud( imread(Img) ),'EdgeColor','none','LineStyle','none');
% end
% 
% axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2); set(gca,'ydir','normal');
% hold on; ax2=axes; h2=show([],elementsFEM(:,1:4),coordinatesFEMWorldDef/um2px,strain_principal_min,'NoEdgeColor');
% set(gca,'fontSize',18); set(gca,'ydir','reverse');view(2); box on; axis equal;  axis tight;   
% alpha(h2,OrigDICImgTransparency); colormap(turbo); caxis auto;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%% TODO: manually modify colormap and caxis %%%%%%
% % colormap(turbo);  caxis auto; % D Sample 
% colormap(turbo); caxis auto % foam
% % colormap(turbo); caxis([-0.008,0]); % Sample 12 
% % colormap(turbo); caxis([-0.2,0.2])
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% linkaxes([ax1,ax2]);  % Link axes together
% ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = []; % Hide the top axes
% colormap(ax1,'gray'); % Give each one its own colormap
% set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
% ax1.Visible = 'on'; ax1.TickLabelInterpreter = 'latex'; 
% %%%%% convert pixel unit to the physical world unit %%%%%
% xticklabels(ax1, num2cell(round(um2px*ax1.XTick*10)/10, length(ax1.XTick) )' );
% yticklabels(ax1, num2cell(round(um2px*ax1.YTick*10)/10, length(ax1.YTick) )' );
% cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'latex';
% 
% 
% 
%  
% % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % ====== 6) Strain e_max_shear ======
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fig1=figure; ax1=axes; 
% try h1=imshow( flipud(imread(Img)),'InitialMagnification','fit');
% catch h1=surf(  flipud( imread(Img) ),'EdgeColor','none','LineStyle','none');
% end
% 
% axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2); set(gca,'ydir','normal');
% hold on; ax2=axes; h2=show([],elementsFEM(:,1:4),coordinatesFEMWorldDef/um2px,strain_maxshear,'NoEdgeColor');
% set(gca,'fontSize',18); set(gca,'ydir','reverse');view(2); box on; axis equal;  axis tight;   
% alpha(h2,OrigDICImgTransparency); colormap(turbo); caxis auto;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%% TODO: manually modify colormap and caxis %%%%%%
% % colormap(turbo);  caxis auto; % D Sample 
% colormap(turbo); caxis auto % foam
% % colormap(turbo); caxis([0,0.011]); % Sample 12 
% % colormap(turbo); caxis([0,0.2])
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% linkaxes([ax1,ax2]);  % Link axes together
% ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = []; % Hide the top axes
% colormap(ax1,'gray'); % Give each one its own colormap
% set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
% ax1.Visible = 'on'; ax1.TickLabelInterpreter = 'latex'; 
% %%%%% convert pixel unit to the physical world unit %%%%%
% xticklabels(ax1, num2cell(round(um2px*ax1.XTick*10)/10, length(ax1.XTick) )' );
% yticklabels(ax1, num2cell(round(um2px*ax1.YTick*10)/10, length(ax1.YTick) )' );
% cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'latex';
% 
%  
 
% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % ====== 7) von Mises equivalent strain ======
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fig1=figure; ax1=axes; 
% try h1=imshow( flipud(imread(Img)),'InitialMagnification','fit');
% catch h1=surf( flipud( imread(Img) ),'EdgeColor','none','LineStyle','none');
% end
% 
% axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2); set(gca,'ydir','normal');
% hold on; ax2=axes; h2=show([],elementsFEM(:,1:4),coordinatesFEMWorldDef/um2px,strain_vonMises,'NoEdgeColor');
% set(gca,'fontSize',18); set(gca,'ydir','reverse'); view(2); box on; axis equal;  axis tight;   
% alpha(h2,OrigDICImgTransparency); colormap(turbo);; caxis auto;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%% TODO: manually modify colormap and caxis %%%%%%
% caxis([0,1]);
% % colormap(turbo); caxis([-0.15,0]) % D Sample 
% % colormap(turbo);  caxis([-0.25,0.25])% foam
% % colormap(turbo); caxis([-0.002,0.017]); % Sample 12 
% %colormap(black_rainbow_plus);  caxis([-0.15 0.15]);
% %  colormap(black_rainbow); caxis([-0.0021,0.0021]);
% % ax1.XTick = [100,200,300]; % Unit: px
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% linkaxes([ax1,ax2]);  % Link axes together
% ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = []; % Hide the top axes
% colormap(ax1,'gray'); % Give each one its own colormap
% set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
% ax1.Visible = 'on';% ax1.TickLabelInterpreter = 'latex'; 
% %%%%% convert pixel unit to the physical world unit %%%%%
% xticklabels(ax1, num2cell(round(um2px*ax1.XTick*10)/10, length(ax1.XTick) )' );
% yticklabels(ax1, num2cell(round(um2px*ax1.YTick*10)/10, length(ax1.YTick) )' );
% % cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'latex';
% cb2 = colorbar('Position',[.17+0.685+0.012 .11+.128 .03 .557 ]); %cb2.TickLabelInterpreter = 'latex';

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== 8) dwdx ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig1=figure; ax1=axes; 
try h1=imshow( flipud(imread(Img)),'InitialMagnification','fit');
catch h1=surf( flipud( imread(Img) ),'EdgeColor','none','LineStyle','none');
end

axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2); set(gca,'ydir','normal');
hold on; ax2=axes; h2=show([],elementsFEM(:,1:4),coordinatesFEMWorldDef/um2px,dwdx,'NoEdgeColor');
set(gca,'fontSize',18); set(gca,'ydir','reverse'); view(2); box on; axis equal;  axis tight;   
alpha(h2,OrigDICImgTransparency); colormap(turbo); caxis auto;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% TODO: manually modify colormap and caxis %%%%%%
%caxis([0,1]);
% colormap(turbo); caxis([-0.15,0]) % D Sample 
% colormap(turbo);  caxis([-0.25,0.25])% foam
 % colormap(turbo);  caxis([-0.02 0.02]) % D shapoed
% caxis([-0.13 0.12])  % bulge
%caxis([-0.2 0.35]) 
% colormap(turbo); caxis([-0.002,0.017]); % Sample 12 
%colormap(black_rainbow_plus);  caxis([-0.15 0.15]);
%  colormap(black_rainbow); caxis([-0.0021,0.0021]);
% ax1.XTick = [100,200,300]; % Unit: px
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

linkaxes([ax1,ax2]);  % Link axes together
ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = []; % Hide the top axes
colormap(ax1,'gray'); % Give each one its own colormap
set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
ax1.Visible = 'on';% ax1.TickLabelInterpreter = 'latex'; 
%%%%% convert pixel unit to the physical world unit %%%%%
xticklabels(ax1, num2cell(round(um2px*ax1.XTick*10)/10, length(ax1.XTick) )' );
yticklabels(ax1, num2cell(round(um2px*ax1.YTick*10)/10, length(ax1.YTick) )' );
% cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'latex';
cb2 = colorbar('Position',[.17+0.685+0.012 .11+.128 .03 .557 ]); %cb2.TickLabelInterpreter = 'latex';


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== 9) dwdy ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig1=figure; ax1=axes; 
try h1=imshow( flipud(imread(Img)),'InitialMagnification','fit');
catch h1=surf( flipud( imread(Img) ),'EdgeColor','none','LineStyle','none');
end

axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2); set(gca,'ydir','normal');
hold on; ax2=axes; h2=show([],elementsFEM(:,1:4),coordinatesFEMWorldDef/um2px,dwdy,'NoEdgeColor');
set(gca,'fontSize',18); set(gca,'ydir','reverse'); view(2); box on; axis equal;  axis tight;   
alpha(h2,OrigDICImgTransparency); colormap(turbo);  caxis auto;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% TODO: manually modify colormap and caxis %%%%%%
% colormap(turbo); caxis([-0.15,0]) % D Sample 
% colormap(turbo);  caxis([-0.25,0.25])% foam
%  colormap(turbo);  caxis([-0.02 0.02]) % D shaped
%caxis([-0.15 0.15]) %bulge
%caxis([-0.2 0.22]) 
% colormap(turbo); caxis([-0.002,0.017]); % Sample 12 
%colormap(black_rainbow_plus);  caxis([-0.15 0.15]);
%  colormap(black_rainbow); caxis([-0.0021,0.0021]);
% ax1.XTick = [100,200,300]; % Unit: px
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

linkaxes([ax1,ax2]);  % Link axes together
ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = []; % Hide the top axes
colormap(ax1,'gray'); % Give each one its own colormap
set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
ax1.Visible = 'on';% ax1.TickLabelInterpreter = 'latex'; 
%%%%% convert pixel unit to the physical world unit %%%%%
xticklabels(ax1, num2cell(round(um2px*ax1.XTick*10)/10, length(ax1.XTick) )' );
yticklabels(ax1, num2cell(round(um2px*ax1.YTick*10)/10, length(ax1.YTick) )' );
% cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'latex';
cb2 = colorbar('Position',[.17+0.685+0.012 .11+.128 .03 .557 ]); %cb2.TickLabelInterpreter = 'latex';




end
 
 
