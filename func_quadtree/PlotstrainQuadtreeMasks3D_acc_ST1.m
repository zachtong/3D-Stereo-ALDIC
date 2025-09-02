function [strain_exx,strain_exy,strain_eyy,strain_principal_max,strain_principal_min, ...
    strain_maxshear,strain_vonMises,dwdx,dwdy] = PlotstrainQuadtreeMasks3D_acc_ST1(U_2D_L,coefficients,voidIndex,FirstFEM,FirstImg,CurrentImg,DICpara)
%PLOTSTRAINQUADTREE: to plot user-selected DIC solved strain fields.

%% --- Initialization, Unit Conversion, and Calculations ---
% This section remains the same as before.
try um2px = DICpara.um2px; 
catch um2px = 1;
end
OrigDICImgTransparency = DICpara.OrigDICImgTransparency;

% (Strain calculation code is exactly the same as before...)
strain_exx = zeros(size(coefficients,1),1);
strain_exy = zeros(size(coefficients,1),1);
strain_eyy = zeros(size(coefficients,1),1);
dwdx = zeros(size(coefficients,1),1);
dwdy = zeros(size(coefficients,1),1);
for i = 1:size(coefficients,1)
    u_x = coefficients{i,1}(1,1); u_y = coefficients{i,1}(2,1);
    v_x = coefficients{i,1}(1,2); v_y = coefficients{i,1}(2,2);
    w_x = coefficients{i,1}(1,3); w_y = coefficients{i,1}(2,3);
    F = [1+u_x, u_y, 0; v_x, 1+v_y, 0; w_x, w_y, 1];
    temp_Strain_tensor = 0.5*(F'*F-eye(3));
    strain_exx(i) = temp_Strain_tensor(1,1);
    strain_exy(i) = temp_Strain_tensor(1,2);
    strain_eyy(i) = temp_Strain_tensor(2,2);
    dwdx(i) = w_x;
    dwdy(i) = w_y;
end

% (Image and coordinate setup is the same...)
if DICpara.Image2PlotResults == 1
    coordinatesFEMWorldDef = [FirstFEM.coordinatesFEM(:,1)+U_2D_L(:,1), FirstFEM.coordinatesFEM(:,2)+U_2D_L(:,2)];
    elementsFEM = FirstFEM.elementsFEM;
    Img = CurrentImg;
else
    coordinatesFEMWorldDef = [FirstFEM.coordinatesFEM(:,1), FirstFEM.coordinatesFEM(:,2)];
    elementsFEM = FirstFEM.elementsFEM;
    Img = FirstImg;
end

% (Derived strain calculation is the same...)
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