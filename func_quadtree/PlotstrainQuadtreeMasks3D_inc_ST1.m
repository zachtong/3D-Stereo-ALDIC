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

%% Compute strain components (extracted to computeStrain3D helper)
[strain_exx, strain_eyy, strain_exy, strain_principal_max, strain_principal_min, ...
 strain_maxshear, strain_vonMises, dwdx, dwdy] = computeStrain3D(coefficients);

%% Plot on deformed images or Not
if DICpara.Image2PlotResults == 1
    % Q6=b: CurrentFEM.coordinatesFEM is already in current-frame image space
    % (the mesh was generated on this frame's mask), so it's the correct
    % plot location. No need to add U_2D_L_inc again.
    coordinatesFEMWorldDef = CurrentFEM.coordinatesFEM;
    elementsFEM = CurrentFEM.elementsFEM;
    Img = CurrentImg;

    % Interpolate strains from FirstFEM node "deformed" positions onto
    % CurrentFEM nodes (both are in current-frame image space).
    % C3.3: scatteredInterpolant('natural') replaces the old RBF thinplate
    % which was O(N^2) memory and 10-100x slower.
    src  = FirstFEM.coordinatesFEM + U_2D_L;  % deformed frame-1 node positions
    keep = ~voidIndex;
    Fi = @(v) scatteredInterpolant(src(keep,1), src(keep,2), v(keep), 'natural', 'none');
    Fi_exx = Fi(strain_exx); Fi_eyy = Fi(strain_eyy); Fi_exy = Fi(strain_exy);
    Fi_dwx = Fi(dwdx);       Fi_dwy = Fi(dwdy);

    q_x = CurrentFEM.coordinatesFEM(:,1);
    q_y = CurrentFEM.coordinatesFEM(:,2);
    strain_exx = Fi_exx(q_x, q_y);
    strain_eyy = Fi_eyy(q_x, q_y);
    strain_exy = Fi_exy(q_x, q_y);
    dwdx       = Fi_dwx(q_x, q_y);
    dwdy       = Fi_dwy(q_x, q_y);

    % Re-derive principal / shear / vonMises from the interpolated base
    % strains (the ones from computeStrain3D were on FirstFEM nodes).
    strain_maxshear      = sqrt((0.5*(strain_exx - strain_eyy)).^2 + strain_exy.^2);
    strain_principal_max = 0.5*(strain_exx + strain_eyy) + strain_maxshear;
    strain_principal_min = 0.5*(strain_exx + strain_eyy) - strain_maxshear;
    strain_vonMises      = sqrt(strain_principal_max.^2 + strain_principal_min.^2 - ...
        strain_principal_max.*strain_principal_min + 3*strain_maxshear.^2);
else
    coordinatesFEMWorldDef = [FirstFEM.coordinatesFEM(:,1), FirstFEM.coordinatesFEM(:,2)];
    elementsFEM = FirstFEM.elementsFEM;
    Img = FirstImg;
end

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
 
 
