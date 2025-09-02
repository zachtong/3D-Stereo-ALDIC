function [strain_exx,strain_exy,strain_eyy,strain_principal_max,strain_principal_min, ...
    strain_maxshear,strain_vonMises] = PlotstrainQuadtreeMasks3D_inc_ST1(U_2D_L_inc,U_2D_L,coefficients,voidIndex,FirstFEM,CurrentFEM,FirstImg,CurrentImg,CurrentImgMask,DICpara)
%PLOTSTRAINQUADTREE: to plot DIC solved strain fields on a quadtree mesh 

%% Initialization

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% convert pixel unit to the physical world unit %%%%%
try um2px = DICpara.um2px; 
catch um2px = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OrigDICImgTransparency = DICpara.OrigDICImgTransparency; % Original raw DIC image transparency

%% Compute strain components
% Strain_tensor = 0.5* (F'*F - I)
strain_exx = zeros(size(coefficients,1),1);
strain_exy = zeros(size(coefficients,1),1);
strain_eyy = zeros(size(coefficients,1),1);
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
    coordinatesFEMWorldDef = [CurrentFEM.coordinatesFEM(:,1)+U_2D_L_inc(1:2:end), CurrentFEM.coordinatesFEM(:,2)+U_2D_L_inc(2:2:end)];
    elementsFEM = CurrentFEM.elementsFEM;
    Img = CurrentImg;

    % Exclude void region whose strain_e.. are zero! They are no need to be
    % used to create interpolant.
    Coor_temp = FirstFEM.coordinatesFEM'+U_2D_L';
    UsedCoor_temp = Coor_temp(:,~voidIndex);

    op3 = rbfcreate(UsedCoor_temp,strain_exx(~voidIndex)','RBFFunction','thinplate');
    strain_exx = rbfinterp(CurrentFEM.coordinatesFEM',op3)';
    op4 = rbfcreate(UsedCoor_temp,strain_exy(~voidIndex)','RBFFunction','thinplate');
    strain_exy = rbfinterp(CurrentFEM.coordinatesFEM',op4)';
    op5 = rbfcreate(UsedCoor_temp,strain_eyy(~voidIndex)','RBFFunction','thinplate');
    strain_eyy = rbfinterp(CurrentFEM.coordinatesFEM',op5)';
    op6 = rbfcreate(UsedCoor_temp,dwdx(~voidIndex)','RBFFunction','thinplate');
    dwdx = rbfinterp(CurrentFEM.coordinatesFEM',op6)';
    op7 = rbfcreate(UsedCoor_temp,dwdy(~voidIndex)','RBFFunction','thinplate');
    dwdy = rbfinterp(CurrentFEM.coordinatesFEM',op7)';
else
    coordinatesFEMWorldDef = [FirstFEM.coordinatesFEM(:,1), FirstFEM.coordinatesFEM(:,2)];
    elementsFEM = FirstFEM.elementsFEM;
    Img = FirstImg;
end

strain_maxshear = sqrt((0.5*(strain_exx-strain_eyy)).^2 + strain_exy.^2);
strain_principal_max = 0.5*(strain_exx+strain_eyy) + strain_maxshear;
strain_principal_min = 0.5*(strain_exx+strain_eyy) - strain_maxshear;
strain_vonMises = sqrt(strain_principal_max.^2 + strain_principal_min.^2 - ...
             strain_principal_max.*strain_principal_min + 3*strain_maxshear.^2);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== Selective Plotting ======
% This section is now a clean, configurable loop.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check which plots the user wants to generate
if isfield(DICpara, 'plots_strain_to_generate')
    plots_to_run = DICpara.plots_strain_to_generate;
else
    % If user doesn't specify, run a default set
    plots_to_run = {'exx', 'eyy', 'exy'};
    disp("Info: 'DICpara.plots_to_generate' is not set.");
    disp("Plotting default set: {'exx', 'eyy', 'exy'}.");
end

% Loop through the user's choices and plot them
for i = 1:length(plots_to_run)
    
    current_plot = plots_to_run{i};
    data_to_plot = [];
    caxis_field_name = '';
    colormap_field_name = '';

    % Use a switch to select the correct data and parameter names
    % This is a clean way to handle different cases without ugly if-else chains.
    switch current_plot
        case 'exx'
            data_to_plot = strain_exx;
            caxis_field_name = 'caxis_exx';
            colormap_field_name = 'colormap_exx';
        case 'exy'
            data_to_plot = strain_exy;
            caxis_field_name = 'caxis_exy';
            colormap_field_name = 'colormap_exy';
        case 'eyy'
            data_to_plot = strain_eyy;
            caxis_field_name = 'caxis_eyy';
            colormap_field_name = 'colormap_eyy';
        case 'principal_max'
            data_to_plot = strain_principal_max;
            caxis_field_name = 'caxis_principal_max';
            colormap_field_name = 'colormap_principal_max';
        case 'principal_min'
            data_to_plot = strain_principal_min;
            caxis_field_name = 'caxis_principal_min';
            colormap_field_name = 'colormap_principal_min';
        case 'max_shear'
            data_to_plot = strain_maxshear;
            caxis_field_name = 'caxis_max_shear';
            colormap_field_name = 'colormap_max_shear';
        case 'vonMises'
            data_to_plot = strain_vonMises;
            caxis_field_name = 'caxis_vonMises';
            colormap_field_name = 'colormap_vonMises';
        case 'dwdx'
            data_to_plot = dwdx;
            caxis_field_name = 'caxis_dwdx';
            colormap_field_name = 'colormap_dwdx';
        case 'dwdy'
            data_to_plot = dwdy;
            caxis_field_name = 'caxis_dwdy';
            colormap_field_name = 'colormap_dwdy';
        otherwise
            warning(['Plot type ''', current_plot, ''' is not recognized. Skipping.']);
            continue; % Skip to the next iteration
    end
    
    % Assemble parameters for the helper function
    plotParams.caxis = 'auto';
    plotParams.colormap = 'turbo';
    if isfield(DICpara, caxis_field_name)
        plotParams.caxis = DICpara.(caxis_field_name);
    end
    if isfield(DICpara, colormap_field_name)
        plotParams.colormap = DICpara.(colormap_field_name);
    end
    
    % Call the one, true plotting function
    plot_displacement_component(Img, elementsFEM, coordinatesFEMWorldDef, data_to_plot, um2px, OrigDICImgTransparency, plotParams);
    
end
end
 
 
