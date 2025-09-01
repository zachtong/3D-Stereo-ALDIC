function PlotWarpage(Warpage,U_2D_L,FirstFEM,FirstImg,CurrentImg,DICpara,voidIndex)

%% Initialization
warning off; load('./plotFiles/colormap_RdYlBu.mat','cMap');
% run('./plotFiles/Black_rainbow.m');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% convert pixel unit to the physical world unit %%%%%
try um2px = DICpara.um2px; 
catch um2px = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OrigDICImgTransparency = DICpara.OrigDICImgTransparency; % Original raw DIC image transparency

Image2PlotResults = DICpara.Image2PlotResults; % Choose image to plot over (first only, second and next images)

%% Plot on deformed images or Not
if Image2PlotResults == 1
    coordinatesFEMWorldDef = [FirstFEM.coordinatesFEM(:,1) + U_2D_L(:,1), FirstFEM.coordinatesFEM(:,2) + U_2D_L(:,2)];
    elementsFEM = FirstFEM.elementsFEM;
    Img = CurrentImg;
   
else
    coordinatesFEMWorldDef = [FirstFEM.coordinatesFEM(:,1), FirstFEM.coordinatesFEM(:,2)];
    elementsFEM = FirstFEM.elementsFEM;
    Img = FirstImg;
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== 1) warpage ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig1=figure; ax1=axes;  
try h1=imshow( flipud(imread(Img)) ,'InitialMagnification','fit');
catch h1=surf(  flipud( imread(Img) ),'EdgeColor','none','LineStyle','none');
end

axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2); set(gca,'ydir','normal');
hold on; ax2=axes; h2=show([], elementsFEM(:,1:4), coordinatesFEMWorldDef/um2px, Warpage{1}, 'NoEdgeColor');
set(gca,'fontSize',18); set(gca,'ydir','reverse');view(2); box on; axis equal;  axis tight;   
alpha(h2,OrigDICImgTransparency);  colormap("turbo"); caxis auto;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% TODO: manually modify colormap and caxis %%%%%%
% colormap(cMap); caxis([-40,40]); % colormap(jet);  
% caxis([-35,35]); % caxis([-0.025,0.025]); 
% caxis([-10,10]); D shaped
%caxis([-0.08,0.03]);  % Bulge
%caxis([-2.5 0]);  % Pig heart

 %caxis([-1 0.5]);
% colormap(black_rainbow);  
%  colormap(jet); caxis([-20 20]);
% ax1.XTick = [100,200,300]; % Unit: px
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

linkaxes([ax1,ax2]);  % Link axes together
ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = []; % Hide the top axes
colormap(ax1,'gray'); % Give each one its own colormap
set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
ax1.Visible = 'on'; % ax1.TickLabelInterpreter = 'latex'; 

%%%%% convert pixel unit to the physical world unit %%%%%
xticklabels(ax1, num2cell(round(um2px*ax1.XTick*10)/10, length(ax1.XTick) )' );
yticklabels(ax1, num2cell(round(um2px*ax1.YTick*10)/10, length(ax1.YTick))');
cb2 = colorbar('Position',[.17+0.685+0.012 .11+.128 .03 .557 ]); % cb2.TickLabelInterpreter = 'latex';
% cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'latex';
 