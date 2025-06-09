function [x2,y2,disp_x,disp_y,...
    dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz,...
    strain_exx, strain_eyy, strain_ezz, strain_exy, strain_eyz, strain_exz, ...
    strain_principal_max,strain_principal_min,strain_maxshear,strain_vonMises] = Plotstrain( ...
    U_2D,coefficients,Rad,x0,y0,sizeOfImg,CurrentImg,DICpara)
% PLOTSTRAIN: to compute and plot DIC solved strain fields on the original DIC images
%
%   INPUT:
%
%   OUTPUT: x2,y2                   x- and y-coordinates of points whose strain values are computed
%           disp_u,disp_v           Interpolated dispu and dispv at points {x2,y2}
%
%   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
warning off; load('./plotFiles/colormap_RdYlBu.mat','cMap');

OrigDICImgTransparency = DICpara.OrigDICImgTransparency; % Original raw DIC image transparency
Image2PlotResults = DICpara.Image2PlotResults; % Choose image to plot over (first only, second and next images)

%% For task1, I commented it
% if DICpara.DICIncOrNot == 0
%     temp(:,1) = U_2D.U(1:2:end,:);
%     temp(:,2) = U_2D.U(2:2:end,:);
%     U_2D = temp;
% end


%% Compute strain components
% initialization
x = x0(1+Rad:end-Rad,1+Rad:end-Rad);
y = y0(1+Rad:end-Rad,1+Rad:end-Rad);

M = size(x,1); N = size(x,2);


u_x = zeros(size(coefficients,1),1); u_y = zeros(size(coefficients,1),1); u_z = zeros(size(coefficients,1),1);
v_x = zeros(size(coefficients,1),1); v_y = zeros(size(coefficients,1),1); v_z = zeros(size(coefficients,1),1);
w_x = zeros(size(coefficients,1),1); w_y = zeros(size(coefficients,1),1); w_z = zeros(size(coefficients,1),1);
for i = 1:size(coefficients,1)
    u_x(i) = coefficients{i,1}(1,1); u_y(i) = coefficients{i,1}(2,1); u_z(i) = coefficients{i,1}(3,1);
    v_x(i) = coefficients{i,1}(1,2); v_y(i) = coefficients{i,1}(2,2); v_z(i) = coefficients{i,1}(3,2);
    w_x(i) = coefficients{i,1}(1,3); w_y(i) = coefficients{i,1}(2,3); w_z(i) = coefficients{i,1}(3,3);
end

if DICpara.winstepsize > 4 % When step size is smaller than 4, no need to refine 
    % imagesc([x(1,1) x(end,1)], [y(1,1) y(1,end)], flipud(g)); hold on;
    if M < 9, x2 = x(:,1)'; else x2 = linspace(x(1,1),x(end,1),4*(length(x(:,1))-1)+1); x2=x2(:)'; end
    if N < 9, y2 = y(1,:); else y2 = linspace(y(1,1),y(1,end),4*(length(y(1,:))-1)+1); y2=y2(:)'; end

    %% Compute displacement components to manipulate the reference image
    x_reshaped = reshape (U_2D(:,1),M+2*Rad,N+2*Rad);
    x_withStrainData =  x_reshaped(Rad+1:end-Rad, Rad+1:end-Rad);
    y_reshaped = reshape (U_2D(:,2),M+2*Rad,N+2*Rad);
    y_withStrainData =  y_reshaped(Rad+1:end-Rad, Rad+1:end-Rad);

    disp_x = gridfit(reshape(x,M*N,1),reshape(y,M*N,1),reshape(x_withStrainData,M,N,1),x2,y2);
    disp_y = gridfit(reshape(x,M*N,1),reshape(y,M*N,1),reshape(y_withStrainData,M,N,1),x2,y2);

    %% Compute strain components
    dudx = gridfit(reshape(x,M*N,1),reshape(y,M*N,1),reshape(u_x,M*N,1),x2,y2);
    dudy = gridfit(reshape(x,M*N,1),reshape(y,M*N,1),reshape(u_y,M*N,1),x2,y2);
    dudz = gridfit(reshape(x,M*N,1),reshape(y,M*N,1),reshape(u_z,M*N,1),x2,y2);
    dvdx = gridfit(reshape(x,M*N,1),reshape(y,M*N,1),reshape(v_x,M*N,1),x2,y2);
    dvdy = gridfit(reshape(x,M*N,1),reshape(y,M*N,1),reshape(v_y,M*N,1),x2,y2);
    dvdz = gridfit(reshape(x,M*N,1),reshape(y,M*N,1),reshape(v_z,M*N,1),x2,y2);
    dwdx = gridfit(reshape(x,M*N,1),reshape(y,M*N,1),reshape(w_x,M*N,1),x2,y2);
    dwdy = gridfit(reshape(x,M*N,1),reshape(y,M*N,1),reshape(w_y,M*N,1),x2,y2);
    dwdz = gridfit(reshape(x,M*N,1),reshape(y,M*N,1),reshape(w_z,M*N,1),x2,y2);


    % Please don't delete this line, to deal with the image and physical world coordinates
    [x2,y2]=ndgrid(x2,y2); x2=x2'; y2=y2';

else
    x2 = x';
    y2 = y';

    x_reshaped = reshape (U_2D(:,1),M+2*Rad,N+2*Rad);
    x_withStrainData =  x_reshaped(Rad+1:end-Rad, Rad+1:end-Rad);
    %x_withStrainData =  x_reshaped(Rad+5:end-Rad, Rad+5:end-Rad);
    y_reshaped = reshape (U_2D(:,2),M+2*Rad,N+2*Rad);
    y_withStrainData =  y_reshaped(Rad+1:end-Rad, Rad+1:end-Rad);
    %y_withStrainData =  y_reshaped(Rad+5:end-Rad, Rad+5:end-Rad);
    disp_x = x_withStrainData';
    disp_y = y_withStrainData';

    dudx = reshape(u_x,[M,N])';
    dudy = reshape(u_y,[M,N])';
    dudz = reshape(u_z,[M,N])';
    dvdx = reshape(v_x,[M,N])';
    dvdy = reshape(v_y,[M,N])';
    dvdz = reshape(v_z,[M,N])';
    dwdx = reshape(w_x,[M,N])';
    dwdy = reshape(w_y,[M,N])';
    dwdz = reshape(w_z,[M,N])';
end
% Method 1
% strain_exx = dudx; strain_eyy = dvdy; strain_ezz = dwdz;
% strain_exy = 0.5*(dvdx + dudy);
% strain_exz = 0.5*(dwdx + dudz);
% strain_eyz = 0.5*(dvdz + dwdy);
%clear u_x u_y u_z v_x v_y v_z w_x w_y w_z

% Method 2: Strain_tensor = 0.5* (F'*F - I)
strain_exx = zeros(size(dudx,1),size(dudx,2));
strain_exy = zeros(size(dudx,1),size(dudx,2));
strain_eyy = zeros(size(dudx,1),size(dudx,2));
strain_ezz = zeros(size(dudx,1),size(dudx,2));
strain_exz = zeros(size(dudx,1),size(dudx,2));
strain_eyz = zeros(size(dudx,1),size(dudx,2));

for i = 1:size(dudx,1)
    for j = 1:size(dudx,2)
        F = [1+dudx(i,j), dudy(i,j), dudz(i,j); dvdx(i,j), 1+dvdy(i,j), dvdz(i,j); dwdx(i,j), dwdy(i,j), 1+dwdz(i,j)];
        temp_Strain_tensor = 0.5*(F'*F-eye(3));
        strain_exx(i,j) = temp_Strain_tensor(1,1);
        strain_exy(i,j) = temp_Strain_tensor(1,2);
        strain_eyy(i,j) = temp_Strain_tensor(2,2);
    end
end

% Shear strain
strain_maxshear = sqrt((0.5*(strain_exx-strain_eyy)).^2 + strain_exy.^2);
% Principal strain
strain_principal_max = 0.5*(strain_exx+strain_eyy) + strain_maxshear;
strain_principal_min = 0.5*(strain_exx+strain_eyy) - strain_maxshear;
% equivalent von Mises strain
strain_vonMises = sqrt(strain_principal_max.^2 + strain_principal_min.^2 - ...
    strain_principal_max.*strain_principal_min + 3*strain_maxshear.^2);

% Debug
% x2 = x2(5:end,5:end);
% y2 = y2(5:end,5:end);
% strain_exx = strain_exx(5:end,5:end);
% strain_eyy = strain_eyy(5:end,5:end);





%% Crop option:
% h2=surf( (x2(:,21:89)+disp_x(:,21:89)), ...
%     sizeOfImg(2)+1-(y2(:,21:89)+disp_y(:,21:89)), strain_exx(:,21:89), 'EdgeColor','none','LineStyle','none');



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== 1) Strain exx ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig1=figure; ax1=axes;
try h1=imshow( flipud(imread(CurrentImg) ),'InitialMagnification','fit');
catch h1=surf( flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
end

axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2); set(gca,'ydir','normal');
hold on; ax2=axes; h2=surf( (x2+Image2PlotResults*disp_x), ...
    sizeOfImg(2)+1-(y2+Image2PlotResults*disp_y), strain_exx, 'EdgeColor','none','LineStyle','none');
set(gca,'fontSize',18); view(2); box on;  % set(gca,'ydir','normal');
alpha(h2,OrigDICImgTransparency);  axis equal;  axis tight;  colormap(cMap);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% TODO: manually modify colormap and caxis %%%%%%
% colormap(jet); %caxis([-0.025,0.025]);
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


%% Zach Debug

crop = 8;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== 1) Strain exx ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig1=figure; ax1=axes;
try h1=imshow( flipud(imread(CurrentImg) ),'InitialMagnification','fit');
catch h1=surf( flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
end

axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2); set(gca,'ydir','normal');
hold on; ax2=axes; h2=surf( (x2(crop+1:end-crop,5*crop+1:end-5*crop)+Image2PlotResults*disp_x(crop+1:end-crop,5*crop+1:end-5*crop)), ...
    sizeOfImg(2)+1-(y2(crop+1:end-crop,5*crop+1:end-5*crop)+Image2PlotResults*disp_y(crop+1:end-crop,5*crop+1:end-5*crop)), strain_exx(crop+1:end-crop,5*crop+1:end-5*crop), 'EdgeColor','none','LineStyle','none');
set(gca,'fontSize',18); view(2); box on;  % set(gca,'ydir','normal');
alpha(h2,OrigDICImgTransparency);  axis equal;  axis tight;  colormap(cMap);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% TODO: manually modify colormap and caxis %%%%%%
% colormap(jet); %caxis([-0.025,0.025]);
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

crop_strain_exx = strain_exx(crop+1:end-crop,crop+1:end-crop);
mean1 = mean(crop_strain_exx(:));
std1 = std(crop_strain_exx(:));
disp(['mean_exx = ', num2str(mean1) ]);
disp(['std_exx = ', num2str(std1) ]);



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== 2) Strain exy ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig1=figure; ax1=axes;
try h1=imshow( flipud(imread(CurrentImg)),'InitialMagnification','fit');
catch h1=surf( flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
end

axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2); set(gca,'ydir','normal');
hold on; ax2=axes; h2=surf( (x2+Image2PlotResults*disp_x), ...
    sizeOfImg(2)+1-(y2+Image2PlotResults*disp_y), strain_exy,'EdgeColor','none','LineStyle','none');
set(gca,'fontSize',18); view(2); box on;  clim auto; % set(gca,'ydir','normal');
alpha(h2,OrigDICImgTransparency);  axis equal;  axis tight; colormap(cMap);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% TODO: manually modify colormap and caxis %%%%%%
% colormap(jet); caxis([-0.008,0.008]);
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


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== 3) Strain eyy ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig1=figure; ax1=axes;
try h1=imshow( flipud(imread(CurrentImg)),'InitialMagnification','fit');
catch h1=surf( flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
end

axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2); set(gca,'ydir','normal');
hold on; ax2=axes; h2=surf( (x2+Image2PlotResults*disp_x), ...
    sizeOfImg(2)+1-(y2+Image2PlotResults*disp_y), strain_eyy,'EdgeColor','none','LineStyle','none');
set(gca,'fontSize',18); view(2); box on;  clim auto; % set(gca,'ydir','normal');
alpha(h2,OrigDICImgTransparency);  axis equal;  axis tight;  colormap(cMap);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% TODO: manually modify colormap and caxis %%%%%%
% colormap(jet); caxis([-0.002,0.014]);
% caxis([-0.55,0]);
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

if 0
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ====== 4) Strain e_principal_max ======
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fig1=figure; ax1=axes;
    try h1=imshow( flipud(imread(CurrentImg)),'InitialMagnification','fit');
    catch h1=surf(  flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
    end

    axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2); set(gca,'ydir','normal');
    hold on; ax2=axes; h2=surf( (x2+Image2PlotResults*disp_x), ...
        sizeOfImg(2)+1-(y2+Image2PlotResults*disp_y), strain_principal_max,'EdgeColor','none','LineStyle','none');
    set(gca,'fontSize',18); view(2); box on;  clim auto; % set(gca,'ydir','normal');
    alpha(h2,OrigDICImgTransparency);  axis equal;  axis tight;  colormap(cMap);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% TODO: manually modify colormap and caxis %%%%%%
    % colormap(jet); % caxis([-0.002,0.014]) % caxis([-0.025,0.025]);
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



    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ====== 5) Strain e_principal_min ======
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fig1=figure; ax1=axes;
    try h1=imshow( flipud(imread(CurrentImg)),'InitialMagnification','fit');
    catch h1=surf(  flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
    end

    axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2); set(gca,'ydir','normal');
    hold on; ax2=axes; h2=surf( (x2+Image2PlotResults*disp_x), ...
        sizeOfImg(2)+1-(y2+Image2PlotResults*disp_y), strain_principal_min,'EdgeColor','none','LineStyle','none');
    set(gca,'fontSize',18); view(2); box on;  clim auto; % set(gca,'ydir','normal');
    alpha(h2,OrigDICImgTransparency);  axis equal;  axis tight; colormap(cMap);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% TODO: manually modify colormap and caxis %%%%%%
    % colormap(jet); % caxis([-0.025,0.025]);
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


    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ====== 6) Strain e_max_shear ======
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fig1=figure; ax1=axes;
    try h1=imshow( flipud(imread(CurrentImg)),'InitialMagnification','fit');
    catch h1=surf(  flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
    end

    axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2); set(gca,'ydir','normal');
    hold on; ax2=axes; h2=surf( (x2+Image2PlotResults*disp_x), ...
        sizeOfImg(2)+1-(y2+Image2PlotResults*disp_y), strain_maxshear,'EdgeColor','none','LineStyle','none');
    set(gca,'fontSize',18); view(2); box on;  clim auto; % set(gca,'ydir','normal');
    alpha(h2,OrigDICImgTransparency);  axis equal;  axis tight; colormap(cMap);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% TODO: manually modify colormap and caxis %%%%%%
    % colormap(jet); caxis([-0.0,0.01]); % caxis([-0.025,0.025]);
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


    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ====== 7) von Mises equivalent strain ======
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fig1=figure; ax1=axes;
    try h1=imshow( flipud(imread(CurrentImg)),'InitialMagnification','fit');
    catch h1=surf(  flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
    end

    axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2);  set(gca,'ydir','normal');
    hold on; ax2=axes; h2=surf( (x2+Image2PlotResults*disp_x), ...
        sizeOfImg(2)+1-(y2+Image2PlotResults*disp_y), strain_vonMises,'EdgeColor','none','LineStyle','none');
    set(gca,'fontSize',18); view(2); box on;  clim auto;
    alpha(h2,OrigDICImgTransparency);  axis equal;  axis tight;  colormap(cMap);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% TODO: manually modify colormap and caxis %%%%%%
    % colormap(jet); caxis([-0.0,0.022]) % caxis([-0.025,0.025]);
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
end
end


