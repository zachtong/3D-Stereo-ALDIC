function PlotdispQuadtreeMasks3D_inc_ST1(U_3D,U_2D_L_inc,U_2D_L,CurrentFEM,FirstFEM,FirstImg,CurrentImg,CurrentImgMask,DICpara,voidIndex)


%% Initialization
warning off; 

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
    coordinatesFEMWorldDef = [CurrentFEM.coordinatesFEM(:,1)+U_2D_L_inc(1:2:end), CurrentFEM.coordinatesFEM(:,2)+U_2D_L_inc(2:2:end)];
    elementsFEM = CurrentFEM.elementsFEM;
    Img = CurrentImg;

    % Exclude void region whose strain_e.. are zero! They are no need to be
    % used to create interpolant.
    Coor_temp = FirstFEM.coordinatesFEM'+U_2D_L';
    UsedCoor_temp = Coor_temp(:,~voidIndex);

    op3 = rbfcreate(UsedCoor_temp,U_3D{1}(~voidIndex)','RBFFunction','thinplate');
    disp_u = rbfinterp(CurrentFEM.coordinatesFEM',op3)';
    op4 = rbfcreate(UsedCoor_temp,U_3D{2}(~voidIndex)','RBFFunction','thinplate');
    disp_v = rbfinterp(CurrentFEM.coordinatesFEM',op4)';
    op5 = rbfcreate(UsedCoor_temp,U_3D{3}(~voidIndex)','RBFFunction','thinplate');
    disp_w = rbfinterp(CurrentFEM.coordinatesFEM',op5)';
else
    coordinatesFEMWorldDef = [FirstFEM.coordinatesFEM(:,1), FirstFEM.coordinatesFEM(:,2)];
    elementsFEM = FirstFEM.elementsFEM;
    Img = FirstImg;
    disp_u = U_3D{1}; disp_v = U_3D{2}; disp_w = U_3D{3};
end

%%%%%%%%%%% JY!!!Mask START %%%%%%%%%%%%%%%
% if Image2PlotResults == 1
%     if ~isempty(CurrentImgMask)
%         for tempi = 1:size(coordinatesFEMWorldDef,1)
%             try
%                 if CurrentImgMask( round(coordinatesFEMWorldDef(tempi,1)/um2px), ...
%                                     (size(CurrentImgMask,2)+1-round(coordinatesFEMWorldDef(tempi,2)/um2px)) ) == 0 
%                     coordinatesFEMWorldDef(tempi,:) = [nan,nan];
%                 end
%             catch
%                 coordinatesFEMWorldDef(tempi,:) = [nan,nan];
%             end
% 
%         end
%     else
%         CurrentImgMask = imread(CurrentImg)';
%         for tempi = 1:size(coordinatesFEMWorldDef,1)
%             try
%                 if CurrentImgMask( round(coordinatesFEMWorldDef(tempi,1)/um2px), ...
%                         (size(CurrentImgMask,2)+1-round(coordinatesFEMWorldDef(tempi,2)/um2px)) ) < 0
%                     coordinatesFEMWorldDef(tempi,:) = [nan,nan];
%                 end
%             catch
%                 coordinatesFEMWorldDef(tempi,:) = [nan,nan];
%             end
%         end
%     end
% end
%%%%%%%%%%% JY!!!Mask END %%%%%%%%%%%%%%%


% Check which plots the user wants to generate
if isfield(DICpara, 'plots_disp_to_generate')
    plots_to_run = DICpara.plots_disp_to_generate;
else
    % If user doesn't specify, run a default set
    plots_to_run = {'u', 'v', 'w'};
    disp("Info: 'DICpara.plots_to_generate' is not set.");
    disp("Plotting default set: {'u', 'v', 'w'}.");
end

% Loop through the user's choices and plot them
for i = 1:length(plots_to_run)
    
    current_plot = plots_to_run{i};
    data_to_plot = [];
    caxis_field_name = '';
    colormap_field_name = '';

    % Use a switch to select the correct data and parameter names
    switch current_plot
        case 'u'
            data_to_plot = disp_u;
            caxis_field_name = 'caxis_u';
            colormap_field_name = 'colormap_u';
        case 'v'
            data_to_plot = disp_v;
            caxis_field_name = 'caxis_v';
            colormap_field_name = 'colormap_v';
        case 'w'
            data_to_plot = disp_w;
            caxis_field_name = 'caxis_w';
            colormap_field_name = 'colormap_w';
        case 'magnitude'
            data_to_plot = sqrt(disp_u.^2 + disp_v.^2 + disp_w.^2);
            caxis_field_name = 'caxis_mag';
            colormap_field_name = 'colormap_mag';
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


